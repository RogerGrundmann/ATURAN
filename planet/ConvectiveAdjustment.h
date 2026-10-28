/*
 * SHARED PHYSICS — dry convective adjustment, one implementation for every planet.
 *
 * THIS FILE MUST BE BYTE-IDENTICAL IN EVERY MODEL THAT USES IT. It knows nothing about which
 * planet it is running on: everything planet-specific arrives through the Model template
 * parameter, and the only two things it asks of a model beyond the usual fields are
 * Model::planet_tag() for the log prefix and the environment-variable prefix built from it.
 *
 * WHY IT IS SHARED. ConvectiveAdjustmentJup.h/.cpp and ConvectiveAdjustmentSat.h/.cpp were 126
 * and 124 lines of which, after normalising the names away, the CODE was identical — every one
 * of the 44 differing lines was a comment. Two copies of an algorithm is two places to fix a
 * bug and one place to forget. ATNEPT and ATURAN have no convective adjustment at all yet;
 * without this file, adding it would make four copies.
 *
 * NOT a preprocessor planet switch. #ifdef JUPITER / #elif SATURN would compile three of the
 * four branches OUT of every build, so a mistake in the Uranus branch would sit undetected until
 * someone next built Uranus. A template is instantiated by every planet that uses it, so a
 * mistake breaks all four builds at once — which is the whole point.
 *
 * Reference algorithm:
 *   Manabe, S. and Strickler, R. F.: "Thermal Equilibrium of the Atmosphere with a Convective
 *   Adjustment", J. Atmos. Sci. 21, 361-385, 1964.
 *
 * WHY A MODEL NEEDS THIS. Measured on ATJUP from its own restart files against its own cp_mix:
 * it develops a thin superadiabatic layer and deepens it as a run proceeds — at iteration 100
 * two levels around 21-24 km exceed the dry adiabat by 0.016 K/km, at iteration 200 by 0.035, by
 * iteration 500 five levels from 21 to 35 km by up to 0.125. It can accumulate that because the
 * momentum equation never felt the buoyancy: with the body forces in the wrong unit system
 * nothing responded to the instability, so nothing relieved it.
 *
 * The buoyancy term cannot fix this by itself, and it is worth being clear about why. It is
 * written as an anomaly about the horizontal mean of each level, so a purely one-dimensional
 * superadiabatic column produces exactly zero force. The term responds to horizontal density
 * contrasts; the unstable stratification is what makes those contrasts grow rather than
 * oscillate. Nothing else in these models restores a column to its adiabat.
 *
 * WHAT IT DOES. Each column is swept from the bottom up. Wherever a layer pair is steeper than
 * the dry adiabat, the whole unstable SEGMENT is mixed to exactly the adiabatic lapse rate,
 * conserving the mass-weighted enthalpy of the segment. Sweeps repeat until the column is
 * stable. Its two properties are the ones that matter: it removes the instability completely
 * rather than damping it, and it does not create or destroy energy.
 *
 * WHAT IT DELIBERATELY DOES NOT DO. It mixes temperature only, not composition. Real convection
 * carries the species with it, and a moist adjustment would use the saturated adiabat where
 * cloud is present rather than the dry one. Both are reasonable extensions; neither is done
 * in this class, because each is a modelling decision with consequences for the microphysics that
 * already runs in the saturation adjustment and the precipitation scheme. The moist version is
 * the second class in this file, MoistConvectiveAdjustment, and has its own note there.
 */

#pragma once

#include "ATPhys.h"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Shared helpers (environment knobs, saturation vapour pressure) live in ATPhys.h.


template<class Model>
class ConvectiveAdjustment {
public:
    explicit ConvectiveAdjustment(Model& model) : m(model) {}

    void run();

private:
    Model& m;
};


template<class Model>
inline void ConvectiveAdjustment<Model>::run(){
    const char* TAG = Model::planet_tag();
    std::cout << std::endl << "      " << TAG
              << ": ConvectiveAdjustment (dry, Manabe-Strickler)" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ---- The critical lapse rate, as a temperature drop per grid layer ----
    //
    // dT/dz = -g/cp is the dry adiabat; over one layer of thickness dz that is a drop of
    // (g/cp)*dz kelvin, and the model's temperature is stored as T/t_ref, so the number compared
    // against below is (g/cp)*dz/t_ref. cp is the mixture value the model computes for itself,
    // which on ATJUP is cp_mix = 10655.4 J/(kg K), making g/cp = 2.33 K/km: with dz = 3.5 km a
    // stable column falls by at most 8.15 K per layer and anything steeper is adjusted. Do not
    // substitute a textbook cp here — the 12000 J/(kg K) that looks right for H2/He gives
    // 2.07 K/km, and a 0.26 K/km error in the criterion is larger than the whole superadiabatic
    // excess ATJUP develops.
    //
    // <PLANET>_CONV_ADJ_LAPSE scales it, for asking what a different critical lapse rate would
    // do — 0 gives an isothermal criterion, values above 1 make the scheme stricter than the dry
    // adiabat, which is one crude way to stand in for a moist adiabat in a condensing region.
    const double dz_m      = (m.L_atm * 1.0e3) / double(m.im - 1);
    const double lapse_fac = ATPhys::env_double(TAG, "CONV_ADJ_LAPSE", 1.0);
    const double dT_ad_nd  = lapse_fac * (m.g / m.cp_mix) * dz_m / m.t_ref;

    // A column is left alone unless it is superadiabatic by more than this, so that round-off
    // does not make the scheme fire on a column that is already neutral.
    constexpr double tol_nd = 1.0e-12;

    // Sweeps are repeated until the column is stable. A deep unstable layer needs one sweep per
    // layer in the worst case; the cap only exists so a pathological column cannot spin here.
    const int max_pass = [TAG](){
        const int v = ATPhys::env_int(TAG, "CONV_ADJ_PASSES", 64);
        return v > 0 ? v : 64; }();

    // ---- Diagnostics, reduced across the whole grid ----
    long long n_columns_adjusted = 0;    // columns that needed at least one sweep
    long long n_layers_mixed     = 0;    // individual layers put back on the adiabat
    int       max_passes_used    = 0;    // worst column, to show whether the cap was reached
    double    max_dT_K           = 0.0;  // largest temperature change applied, in kelvin
    double    max_rel_drift      = 0.0;  // worst enthalpy non-conservation, relative

    std::vector<double> w(m.im), t_col(m.im);

    #pragma omp parallel for collapse(2) schedule(dynamic, 8) firstprivate(w, t_col) \
        reduction(+:n_columns_adjusted, n_layers_mixed) \
        reduction(max:max_passes_used, max_dT_K, max_rel_drift)
    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){

            // ---- The fluid part of this column ----
            // Solid cells hold boundary values, not a fluid state, so the column starts above the
            // topography and stops at the first solid cell above it (there should be none, but a
            // column that is walled in higher up must not be mixed across the wall). On a model
            // with no obstacle every SeaMount entry is 0 and this reduces to i0 = 0, i1 = im-1.
            int i0 = 0;
            while(i0 < m.im && m.SeaMount.x[i0][j][k] == 1.0) i0++;
            int i1 = i0;
            while(i1 + 1 < m.im && m.SeaMount.x[i1 + 1][j][k] != 1.0) i1++;
            if(i1 - i0 < 1) continue;                       // nothing to mix

            // ---- Mass weights ----
            // The enthalpy per unit area of a layer is (cp/g) * dp, so the conserved quantity is
            // the dp-weighted temperature. dp is taken from the hydrostatic p_stat with faces at
            // the midpoints, and the two end layers get their half-cell. Using dp rather than a
            // density avoids a trap: rho = p/(R*T) makes rho*T identically p/R, so a density
            // weight would conserve nothing at all.
            bool weights_ok = true;
            for(int i = i0; i <= i1; i++){
                const double p_lo = (i > i0) ? 0.5 * (m.p_stat.x[i-1][j][k] + m.p_stat.x[i][j][k])
                                             : m.p_stat.x[i0][j][k];
                const double p_hi = (i < i1) ? 0.5 * (m.p_stat.x[i][j][k] + m.p_stat.x[i+1][j][k])
                                             : m.p_stat.x[i1][j][k];
                w[i] = p_lo - p_hi;
                if(!(w[i] > 0.0) || !std::isfinite(w[i])) weights_ok = false;
                t_col[i] = m.t.x[i][j][k];
                if(!std::isfinite(t_col[i]) || t_col[i] <= 0.0) weights_ok = false;
            }
            // A column with a non-monotonic p_stat or a bad temperature is left untouched: this
            // routine is not the place to repair either, and mixing across a NaN would spread it
            // through the whole column.
            if(!weights_ok) continue;

            double sum_before = 0.0;
            for(int i = i0; i <= i1; i++) sum_before += w[i] * t_col[i];

            // ---- Sweeps ----
            //
            // Whole unstable SEGMENTS are mixed at once, not adjacent pairs. Both converge to the
            // same profile, but pairwise mixing moves heat one layer per sweep, so a deep unstable
            // block needs as many sweeps as it has layers: a test that made the criterion
            // artificially strict (CONV_ADJ_LAPSE=0.5, so nearly every column qualifies) ran 137
            // million pair adjustments and still hit the 64-sweep cap. Mixing the segment settles
            // it in one step.
            //
            // A segment [a..b] is put on the adiabat, T_q = C - dT_ad*(q-a), with C chosen so the
            // dp-weighted temperature of the segment is unchanged:
            //     C = [ sum w_q T_q + dT_ad * sum w_q (q-a) ] / sum w_q
            // Mixing can destabilise the joint with the layer below or above, so the segment is
            // grown in whichever direction is still too steep and re-mixed. It can only grow, and
            // only within the column, so the inner loop terminates.
            int passes = 0;
            long long layers_mixed = 0;
            bool changed = true;
            while(changed && passes < max_pass){
                changed = false;
                passes++;
                int i = i0;
                while(i < i1){
                    if(t_col[i] - t_col[i+1] <= dT_ad_nd + tol_nd){ i++; continue; }

                    int a = i, b = i + 1;
                    for(;;){
                        double sw = 0.0, swt = 0.0, swk = 0.0;
                        for(int q = a; q <= b; q++){
                            sw  += w[q];
                            swt += w[q] * t_col[q];
                            swk += w[q] * double(q - a);
                        }
                        const double C = (swt + dT_ad_nd * swk) / sw;
                        for(int q = a; q <= b; q++){
                            const double t_new = C - dT_ad_nd * double(q - a);
                            const double dK = std::fabs(t_new - t_col[q]) * m.t_ref;
                            if(dK > max_dT_K) max_dT_K = dK;
                            t_col[q] = t_new;
                        }

                        bool extended = false;
                        if(a > i0 && t_col[a-1] - t_col[a] > dT_ad_nd + tol_nd){ a--; extended = true; }
                        if(b < i1 && t_col[b] - t_col[b+1] > dT_ad_nd + tol_nd){ b++; extended = true; }
                        if(!extended) break;
                    }

                    layers_mixed += (b - a + 1);
                    changed = true;
                    i = b;                       // carry on above the block just mixed
                }
            }

            if(layers_mixed == 0) continue;

            double sum_after = 0.0;
            for(int i = i0; i <= i1; i++) sum_after += w[i] * t_col[i];
            const double drift = (sum_before != 0.0)
                               ? std::fabs(sum_after - sum_before) / std::fabs(sum_before) : 0.0;
            if(drift > max_rel_drift) max_rel_drift = drift;

            for(int i = i0; i <= i1; i++) m.t.x[i][j][k] = t_col[i];

            // COUNT ONLY COLUMNS THAT WERE ACTUALLY MIXED. This was unconditional, so it
            // counted every column that reached here — i.e. every column whose weights passed the
            // finiteness guard above, adjusted or not. The variable is named
            // n_columns_adjusted and the report prints it as "N of M columns (P %)", so the
            // percentage read as "how much of the model is convectively unstable" when it meant
            // "how much of the model has usable weights".
            //
            // On ATNEPT that printed 100.00 %, which survived a neutral initial profile, the
            // physics block not running on odd iterations, and the removal of a temperature clamp
            // — because it was jm*km and never varied. On ATSAT it prints 97.79 %, the 2.21 %
            // shortfall being columns that fail the weights guard rather than stable ones.
            //
            // layers_mixed and the max dT are the honest signals and were always right.
            if(layers_mixed > 0) n_columns_adjusted++;
            n_layers_mixed += layers_mixed;
            if(passes > max_passes_used) max_passes_used = passes;
        }
    }

    const double frac = 100.0 * double(n_columns_adjusted) / double(m.jm * m.km);
    printf("      %s: convective adjustment — critical drop %.3f K per %.1f km layer;"
           " %lld of %d columns (%.2f %%), %lld layers, worst column %d sweeps of %d,"
           " max dT %.3f K, enthalpy drift %.2e\n",
           TAG, dT_ad_nd * m.t_ref, dz_m * 1.0e-3,
           n_columns_adjusted, m.jm * m.km, frac, n_layers_mixed,
           max_passes_used, max_pass, max_dT_K, max_rel_drift);
    if(max_passes_used >= max_pass)
        std::cout << "      " << TAG << ": WARNING - the sweep cap was reached, a column may "
                     "still be superadiabatic (raise " << TAG << "_CONV_ADJ_PASSES)" << std::endl;

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for ConvectiveAdjustment\n", elapsed.count() * 1e-9);
    std::cout << "      " << TAG << ": ConvectiveAdjustment ended" << std::endl;
}


/*
 * ===== MOIST CONVECTIVE ADJUSTMENT, LEDOUX CRITERION (2026-10-26) =====
 *
 * <TAG>_MCA=<N> at the model's call site (0 = off, the default; N = every N iterations).
 *
 * WHY. In an H2-He atmosphere the condensable vapours are HEAVIER than the air they sit in (NH3 17,
 * H2O 18, against a mean molar mass of ~2.3), so a vapour content that falls with height makes the
 * column denser below and stabilises it -- the opposite of Earth, where vapour is light. Whether a
 * column convects is therefore decided by DENSITY, not by temperature alone (the Ledoux criterion),
 * and a saturated parcel rises along its moist adiabat, not the dry one. The dry adjustment above
 * knows neither. ATOM_Precipitation has no moist-stability law for such an atmosphere (its Tiedtke
 * scheme works on theta_e with a dry kappa and an Earth parcel); this is written for the giants.
 *
 * WHAT IT DOES. Per column, bottom up, every adjacent pair (a, a+1) is tested with a parcel: the
 * cell a is lifted to the level of a+1 conserving its liquid-ice static energy
 *     h_L = cp T + g z - sum_gas (lv q_liquid + ls q_ice) - sum_tracer L_tr q_tr
 * and the total mass fraction q_t = q_v + q_liquid + q_ice of every gas, and condenses (or
 * evaporates) to equilibrium there. The pair is unstable when the parcel is LIGHTER than the air at
 * a+1, at the same pressure, i.e. when its virtual temperature is higher:
 *     T_v = T (1 + sum_gas q_v (1/ep - 1) - q_condensate)
 * (ep = mu_gas / mu_air, the saturation adjustment's own ep; condensate carries mass, no pressure).
 * An unstable pair is MIXED: h_L, every gas's q_t and every passive tracer are made uniform over the
 * segment, mass-weighted, and every cell of the segment is put back to equilibrium at its own p and
 * z. Uniform h_L and q_t is exactly neutral (a parcel from any cell of it lands on the state of the
 * cell above), dry adiabatic where the segment is subsaturated and moist adiabatic where it is
 * saturated. The segment then grows while its edges are unstable, and sweeps repeat, as in the dry
 * class. A subsaturated column with no composition gradient reduces EXACTLY to the dry criterion
 * (drop g dz / cp per layer).
 *
 * Equilibrium means: below the gas's t_0 the condensate is ice and the vapour is at ice saturation,
 * above it liquid at liquid saturation -- the saturation adjustment's BERGERON=2 phase rule, ATJUP's
 * default. With another rule the next saturation adjustment re-partitions the mixed cells; the
 * adjustment itself conserves h_L and mass whatever the phase. Nothing precipitates here.
 *
 * CONSERVED, per column: sum w h_L and sum w q_t of every gas and tracer, w = rho dz (the mass of
 * the cell per unit area: rho the model's rho_at, dz the face-centred thickness from the layer
 * heights, half cells at the ends) -- so the column mass of every species, sum c dz, is unchanged.
 * The worst relative drift is printed.
 *
 * TRACER ENTHALPY (2026-10-26). A tracer that is the product of a reaction releasing heat L_tr [J/kg]
 * on formation (ATJUP: NH4SH from NH3 + H2S, with its equilibrium adjustment on) enters h_L like a
 * condensate, - L_tr q_tr. The tracer itself is only mixed here, not re-equilibrated; the model's own
 * adjustment does that afterwards with the same L_tr, so the energy of the mixing and of the
 * re-equilibration add up. L_tr = 0 leaves a tracer passive in the energy, as before.
 *
 * WHAT IT DOES NOT DO. Momentum, k and omega are not mixed. The test is LOCAL (adjacent pairs): a
 * parcel from deep below that would have to cross a stable layer first (CAPE behind CIN) does not
 * fire -- that is penetrative convection, which a mass-flux scheme would need. The time scale is
 * instantaneous unless <TAG>_MCA_RELAX < 1.
 *
 * Knobs (read here, under the model's prefix):
 *   <TAG>_MCA_LEDOUX   1  composition in the density test; 0 = temperature only (Schwarzschild,
 *                         still with the moist parcel)
 *   <TAG>_MCA_LOADING  1  condensate weight in the density test (only with LEDOUX)
 *   <TAG>_MCA_RELAX    1  fraction of the adjustment applied per call, on the conserved variables
 *   <TAG>_MCA_PASSES  64  sweep cap per column
 * The report line counts the pairs that are unstable in temperature but held by their composition
 * ("held by composition") -- the number that says what the Ledoux term is worth.
 */

class Array;

template<class Model>
class MoistConvectiveAdjustment {
public:
    // One condensable gas: vapour, liquid and ice as densities [kg/m3], and the coefficient sets the
    // model hands its saturation adjustment for the same gas.
    struct Gas {
        const char* name;
        Array* v; Array* c; Array* ice;
        double ep;                                  // mu_gas / mu_air: q_sat = ep e_sat / p
        double lv, ls;                              // J/kg
        double t_0;                                 // K: liquid above, ice below
        double C, L0, R, del_alf, del_bet;          // liquid set, ATPhys::saturation_vapour_pressure
        double C_i, L0_i, del_alf_i, del_bet_i;     // ice set
    };
    // A field mixed without a phase change, as a density [kg/m3]. ep > 0: a gas, in the virtual
    // temperature with 1/ep - 1; condensed: a solid or liquid, in the loading; neither: passive.
    // L: heat released per kg when the tracer is FORMED [J/kg], carried in h_L (0 = none).
    struct Tracer {
        const char* name;
        Array* f;
        double ep;
        bool condensed;
        double L;
    };

    MoistConvectiveAdjustment(Model& model, std::vector<Gas> gases, std::vector<Tracer> tracers)
        : m(model), G(std::move(gases)), TR(std::move(tracers)) {}

    void run();

private:
    Model& m;
    std::vector<Gas>    G;
    std::vector<Tracer> TR;
};


template<class Model>
inline void MoistConvectiveAdjustment<Model>::run(){
    const char* TAG = Model::planet_tag();
    std::cout << std::endl << "      " << TAG
              << ": MoistConvectiveAdjustment (Ledoux, h_L and q_t mixed)" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    static const bool   ledoux   = ATPhys::env_int(TAG, "MCA_LEDOUX", 1) != 0;
    static const bool   loading  = ATPhys::env_int(TAG, "MCA_LOADING", 1) != 0;
    static const double relax    = std::min(1.0, std::max(0.0, ATPhys::env_double(TAG, "MCA_RELAX", 1.0)));
    static const int    max_pass = [TAG](){ const int v = ATPhys::env_int(TAG, "MCA_PASSES", 64); return v > 0 ? v : 64; }();

    const int    ng = (int)G.size(), nt = (int)TR.size(), im = m.im;
    const double cp = m.cp_mix, grav = m.g;
    constexpr double tol_K = 1.0e-7;                // a pair is unstable when the parcel is warmer than this
    constexpr int    n_iter   = 60;             // cap of the equilibrium solve
    constexpr double tol_T    = 1.0e-10;        // K, its bracket width

    // Mass fraction at saturation over the phase the temperature selects, and its latent heat.
    auto qsat = [&](int g, double T, double p) -> double {
        const Gas& s = G[g];
        const double E = (T >= s.t_0)
            ? ATPhys::saturation_vapour_pressure(T, s.C,   s.L0,   s.R, s.del_alf,   s.del_bet)
            : ATPhys::saturation_vapour_pressure(T, s.C_i, s.L0_i, s.R, s.del_alf_i, s.del_bet_i);
        const double q = s.ep * E / p;
        return std::isfinite(q) ? q : 0.0;
    };
    auto lat = [&](int g, double T){ return T >= G[g].t_0 ? G[g].lv : G[g].ls; };

    // Equilibrium temperature of (h_L, q_t) at (p, z): cp T + g z - sum L(T) max(0, q_t - q_sat(T)) = h_L.
    // The left side rises monotonically with T (also across t_0, where ls > lv), so the root is
    // bracketed by the all-vapour value T_lo (f <= 0) and T_lo + sum L q_t / cp (f >= 0); it is found
    // by the Illinois variant of regula falsi inside that bracket (6-10 steps where a bisection took
    // 60, and the bracket keeps it safe at the kinks of t_0 and of q_t = q_sat). cond[g] receives the
    // condensate. Returns false if the state has no positive temperature.
    auto equil = [&](double hL, const double* qt, double p, double z, double& T, double* cond) -> bool {
        const double T_lo = (hL - grav * z) / cp;
        if(!(T_lo > 0.0) || !(p > 0.0)) return false;
        bool any = false;
        for(int g = 0; g < ng; g++) if(qt[g] > qsat(g, T_lo, p)){ any = true; break; }
        if(!any){ T = T_lo; for(int g = 0; g < ng; g++) cond[g] = 0.0; return true; }
        double span = 0.0;
        for(int g = 0; g < ng; g++) span += std::max(G[g].lv, G[g].ls) * std::max(0.0, qt[g]);
        auto fof = [&](double TT){
            double f = cp * TT + grav * z - hL;
            for(int g = 0; g < ng; g++) f -= lat(g, TT) * std::max(0.0, qt[g] - qsat(g, TT, p));
            return f; };
        double a = T_lo, b = T_lo + span / cp + 1.0e-9;
        double fa = fof(a), fb = fof(b);
        int side = 0;
        for(int it = 0; it < n_iter && b - a > tol_T; it++){
            double c = (fb - fa != 0.0) ? b - fb * (b - a) / (fb - fa) : 0.5 * (a + b);
            if(!(c > a && c < b)) c = 0.5 * (a + b);
            const double fc = fof(c);
            if(fc < 0.0){ a = c; fa = fc; if(side == -1) fb *= 0.5; side = -1; }
            else        { b = c; fb = fc; if(side == +1) fa *= 0.5; side = +1; }
            if(std::fabs(fc) <= cp * tol_T){ a = b = c; break; }
        }
        T = b;                                      // the side where the condensate is not overcounted
        for(int g = 0; g < ng; g++) cond[g] = std::max(0.0, qt[g] - qsat(g, T, p));
        return true;
    };

    static int dbg_j = -1, dbg_k = -1;
    static bool dbg_done = false;
    if(const char* e = ATPhys::env_for(TAG, "MCA_DEBUG")){ if(dbg_j < 0) sscanf(e, "%d,%d", &dbg_j, &dbg_k); }

    long long n_columns = 0, n_layers = 0, n_moist_trig = 0, n_dry_trig = 0, n_held = 0, n_unstable_T = 0, n_skipped = 0;
    int    max_passes_used = 0;
    double max_dT_K = 0.0, max_drift_h = 0.0, max_drift_q = 0.0;

    #pragma omp parallel for collapse(2) schedule(dynamic, 8) \
        reduction(+:n_columns, n_layers, n_moist_trig, n_dry_trig, n_held, n_unstable_T, n_skipped) \
        reduction(max:max_passes_used, max_dT_K, max_drift_h, max_drift_q)
    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){
            int i0 = 0;
            while(i0 < im && m.SeaMount.x[i0][j][k] == 1.0) i0++;
            int i1 = i0;
            while(i1 + 1 < im && m.SeaMount.x[i1 + 1][j][k] != 1.0) i1++;
            // <TAG>_MCA_INTERIOR=1 (2026-10-28, ATJUP row GIANTS-MCA-DEFAULT; default Model::mca_default_interior(),
            // true on ATURAN only, false = the column as before): leave the
            // radial boundary planes i = 0 and i = im-1 out. bcRadius writes them from the interior and
            // the Runge-Kutta step does not integrate them, so mixing them is undone every iteration.
            // MEASURED on ATURAN (satchk/mcagiants_2026-10-28/utb, ATURAN_TB_DIAG, 8 it): the MCA put
            // +7.44 K into i0 and took it from i1-3, the BC took -8.76 K back out of i0, and the pair
            // was unstable again: 100 % of columns on every call. A column above a SeaMount keeps its
            // first fluid cell, which is an integrated cell.
            static const bool interior = ATPhys::env_int(TAG, "MCA_INTERIOR", Model::mca_default_interior() ? 1 : 0) != 0;
            if(interior){ if(i0 == 0) i0 = 1; if(i1 == im - 1) i1 = im - 2; }
            if(i1 - i0 < 1) continue;

            // ---- Column state (mass fractions, K, bar, m) ----
            std::vector<double> T(im), p(im), z(im), rho(im), w(im), hL(im), T0(im);
            std::vector<double> qt(im * ng), qv(im * ng), qc(im * ng), qi(im * ng), qr(im * std::max(nt, 1));
            std::vector<char>   touched(im, 0);
            bool ok = true;
            for(int i = i0; i <= i1 && ok; i++){
                T[i] = T0[i] = m.t.x[i][j][k] * m.t_ref;
                p[i] = m.p_stat.x[i][j][k];
                z[i] = double(m.get_layer_height(i)) * 1.0e3;
                rho[i] = m.rho_at(i, j, k);
                if(!(T[i] > 0.0) || !(p[i] > 0.0) || !(rho[i] > 0.0) || !std::isfinite(T[i] + p[i] + rho[i])) ok = false;
            }
            if(!ok){ n_skipped++; continue; }
            for(int i = i0; i <= i1 && ok; i++){
                const double lo = (i > i0) ? 0.5 * (z[i-1] + z[i]) : z[i0];
                const double hi = (i < i1) ? 0.5 * (z[i] + z[i+1]) : z[i1];
                w[i] = rho[i] * (hi - lo);
                if(!(w[i] > 0.0)) ok = false;
                double h = cp * T[i] + grav * z[i];
                for(int g = 0; g < ng; g++){
                    const int x = i * ng + g;
                    qv[x] = std::max(0.0, G[g].v->x[i][j][k])   / rho[i];
                    qc[x] = std::max(0.0, G[g].c->x[i][j][k])   / rho[i];
                    qi[x] = std::max(0.0, G[g].ice->x[i][j][k]) / rho[i];
                    qt[x] = qv[x] + qc[x] + qi[x];
                    h -= G[g].lv * qc[x] + G[g].ls * qi[x];
                }
                for(int r = 0; r < nt; r++){
                    qr[i * nt + r] = TR[r].f->x[i][j][k] / rho[i];
                    h -= TR[r].L * qr[i * nt + r];
                }
                hL[i] = h;
            }
            if(!ok){ n_skipped++; continue; }

            // The tracers' share of h_L in a cell: equil() solves for T with the gases' latent heat only,
            // so it is handed h_L + this (the tracers are not changed by the solve).
            auto trE = [&](int i_tr) -> double {
                double e = 0.0;
                for(int r = 0; r < nt; r++) e += TR[r].L * qr[i_tr * nt + r];
                return e; };

            // Virtual-temperature factor of a cell state; composition only with LEDOUX.
            auto tv = [&](double TT, const double* v_, const double* cnd, int i_tr) -> double {
                if(!ledoux) return TT;
                double f = 1.0;
                for(int g = 0; g < ng; g++){ f += v_[g] * (1.0 / G[g].ep - 1.0); if(loading) f -= cnd[g]; }
                for(int r = 0; r < nt; r++){
                    const double q = qr[i_tr * nt + r];
                    if(TR[r].ep > 0.0) f += q * (1.0 / TR[r].ep - 1.0);
                    else if(TR[r].condensed && loading) f -= q;
                }
                return TT * f;
            };
            std::vector<double> pc(ng), pv(ng), ec(ng), ev(ng);
            // Each cell's own equilibrium (T, condensate), solved once and kept until the cell is mixed.
            std::vector<double> Teq(im), ceq(im * ng);
            std::vector<char>   eq_ok(im, 0);
            // Parcel from a lifted to b; unstable if lighter. held: unstable in T alone but not in T_v.
            // The air at b is taken in ITS equilibrium too (same h_L and q_t, condensed or evaporated
            // in place), so vapour that is supersaturated only because the saturation adjustment has
            // not run since the last step cannot fire the scheme by its latent heat alone.
            auto unstable = [&](int a, int b, bool* held, bool* moist) -> bool {
                double Tp;
                if(!equil(hL[a] + trE(a), &qt[a * ng], p[b], z[b], Tp, pc.data())) return false;
                if(!eq_ok[b]){
                    if(!equil(hL[b] + trE(b), &qt[b * ng], p[b], z[b], Teq[b], &ceq[b * ng])) return false;
                    eq_ok[b] = 1;
                }
                const double Tb = Teq[b];
                for(int g = 0; g < ng; g++) ec[g] = ceq[b * ng + g];
                bool wet = false;
                for(int g = 0; g < ng; g++){ pv[g] = qt[a * ng + g] - pc[g]; if(pc[g] > 0.0) wet = true; ev[g] = qt[b * ng + g] - ec[g]; }
                const double dTv = tv(Tp, pv.data(), pc.data(), a) - tv(Tb, ev.data(), ec.data(), b);
                if(held) *held = (Tp - Tb > tol_K) && !(dTv > tol_K);
                if(moist) *moist = wet;
                return dTv > tol_K;
            };
            auto settle = [&](int q) -> bool {             // cell q to equilibrium from hL / qt
                double TT;
                if(!equil(hL[q] + trE(q), &qt[q * ng], p[q], z[q], TT, pc.data())) return false;
                T[q] = TT;
                Teq[q] = TT; eq_ok[q] = 1;
                for(int g = 0; g < ng; g++){
                    const int x = q * ng + g;
                    ceq[x] = pc[g];
                    qv[x] = qt[x] - pc[g];
                    if(TT >= G[g].t_0){ qc[x] = pc[g]; qi[x] = 0.0; } else { qi[x] = pc[g]; qc[x] = 0.0; }
                }
                touched[q] = 1;
                return true;
            };

            // <TAG>_MCA_DEBUG=j,k (print-only): the parcel test of one column, first call only.
            if(dbg_j == j && dbg_k == k && !dbg_done){
                dbg_done = true;
                printf("      %s: [MCA-DEBUG] column j %d k %d, cp %.2f g %.3f\n", TAG, j, k, cp, grav);
                for(int i = i0; i < i1; i++){
                    double Tp = 0.0, Tb = 0.0;
                    equil(hL[i] + trE(i), &qt[i * ng], p[i+1], z[i+1], Tp, pc.data());
                    equil(hL[i+1] + trE(i+1), &qt[(i+1) * ng], p[i+1], z[i+1], Tb, ec.data());
                    for(int g = 0; g < ng; g++){ pv[g] = qt[i * ng + g] - pc[g]; ev[g] = qt[(i+1) * ng + g] - ec[g]; }
                    const double tvp = tv(Tp, pv.data(), pc.data(), i), tve = tv(Tb, ev.data(), ec.data(), i + 1);
                    printf("      %s: [MCA-DEBUG] i %2d T %.3f z %.0f p %.4f rho %.4e w %.4e qt0 %.3e qt1 %.3e | parcel %.3f env %.3f dT %+.4f | Tv %.3f %.3f dTv %+.4f\n",
                           TAG, i, T[i], z[i], p[i], rho[i], w[i], qt[i * ng], ng > 1 ? qt[i * ng + 1] : 0.0, Tp, Tb, Tp - Tb, tvp, tve, tvp - tve);
                }
            }

            // Census of the unadjusted column: pairs unstable in temperature, and those held by composition.
            std::vector<char> cen_u(im, 0), cen_m(im, 0);   // kept for the first sweep
            for(int i = i0; i < i1; i++){
                bool held = false, moist = false;
                const bool u = unstable(i, i + 1, &held, &moist);
                cen_u[i] = u; cen_m[i] = moist;
                if(u || held) n_unstable_T++;
                if(held) n_held++;
            }

            std::vector<double> hL0(hL.begin(), hL.end()), qt0(qt.begin(), qt.end()), qr0(qr.begin(), qr.end());
            double sum_h0 = 0.0; std::vector<double> sum_q0(ng + nt, 0.0);
            for(int i = i0; i <= i1; i++){
                sum_h0 += w[i] * hL[i];
                for(int g = 0; g < ng; g++) sum_q0[g] += w[i] * qt[i * ng + g];
                for(int r = 0; r < nt; r++) sum_q0[ng + r] += w[i] * qr[i * nt + r];
            }

            int passes = 0; long long layers = 0; bool changed = true, fail = false;
            while(changed && passes < max_pass && !fail){
                changed = false; passes++;
                int i = i0;
                while(i < i1 && !fail){
                    bool moist = false, uns;
                    if(passes == 1 && !touched[i] && !touched[i + 1]){ uns = cen_u[i]; moist = cen_m[i]; }
                    else uns = unstable(i, i + 1, nullptr, &moist);
                    if(!uns){ i++; continue; }
                    if(passes == 1){ if(moist) n_moist_trig++; else n_dry_trig++; }
                    int a = i, b = i + 1;
                    for(;;){
                        double sw = 0.0, sh = 0.0;
                        std::vector<double> sq(ng, 0.0), sr(nt, 0.0);
                        for(int q = a; q <= b; q++){
                            sw += w[q]; sh += w[q] * hL[q];
                            for(int g = 0; g < ng; g++) sq[g] += w[q] * qt[q * ng + g];
                            for(int r = 0; r < nt; r++) sr[r] += w[q] * qr[q * nt + r];
                        }
                        for(int q = a; q <= b; q++){
                            hL[q] = sh / sw;
                            for(int g = 0; g < ng; g++) qt[q * ng + g] = sq[g] / sw;
                            for(int r = 0; r < nt; r++) qr[q * nt + r] = sr[r] / sw;
                            if(!settle(q)){ fail = true; break; }
                        }
                        if(fail) break;
                        bool ext = false;
                        if(a > i0 && unstable(a - 1, a, nullptr, nullptr)){ a--; ext = true; }
                        if(b < i1 && unstable(b, b + 1, nullptr, nullptr)){ b++; ext = true; }
                        if(!ext) break;
                    }
                    layers += (b - a + 1);
                    changed = true;
                    i = b;
                }
            }
            if(fail){ n_skipped++; continue; }              // nothing written: the column keeps its state
            if(layers == 0) continue;

            // <TAG>_MCA_RELAX: only a fraction of the adjustment, on the conserved variables.
            if(relax < 1.0){
                for(int q = i0; q <= i1; q++){
                    if(!touched[q]) continue;
                    hL[q] = hL0[q] + relax * (hL[q] - hL0[q]);
                    for(int g = 0; g < ng; g++){ const int x = q * ng + g; qt[x] = qt0[x] + relax * (qt[x] - qt0[x]); }
                    for(int r = 0; r < nt; r++){ const int x = q * nt + r; qr[x] = qr0[x] + relax * (qr[x] - qr0[x]); }
                    if(!settle(q)){ fail = true; break; }
                }
                if(fail){ n_skipped++; continue; }
            }

            double sum_h1 = 0.0; std::vector<double> sum_q1(ng + nt, 0.0);
            for(int i = i0; i <= i1; i++){
                sum_h1 += w[i] * hL[i];
                for(int g = 0; g < ng; g++) sum_q1[g] += w[i] * qt[i * ng + g];
                for(int r = 0; r < nt; r++) sum_q1[ng + r] += w[i] * qr[i * nt + r];
            }
            if(sum_h0 != 0.0) max_drift_h = std::max(max_drift_h, std::fabs(sum_h1 - sum_h0) / std::fabs(sum_h0));
            for(int x = 0; x < ng + nt; x++)
                if(sum_q0[x] > 0.0) max_drift_q = std::max(max_drift_q, std::fabs(sum_q1[x] - sum_q0[x]) / sum_q0[x]);

            for(int q = i0; q <= i1; q++){
                if(!touched[q]) continue;
                max_dT_K = std::max(max_dT_K, std::fabs(T[q] - T0[q]));
                m.t.x[q][j][k] = T[q] / m.t_ref;
                for(int g = 0; g < ng; g++){
                    const int x = q * ng + g;
                    G[g].v->x[q][j][k]   = qv[x] * rho[q];
                    G[g].c->x[q][j][k]   = qc[x] * rho[q];
                    G[g].ice->x[q][j][k] = qi[x] * rho[q];
                }
                for(int r = 0; r < nt; r++) TR[r].f->x[q][j][k] = qr[q * nt + r] * rho[q];
            }
            n_columns++;
            n_layers += layers;
            if(passes > max_passes_used) max_passes_used = passes;
        }
    }

    printf("      %s: moist convective adjustment (LEDOUX %d, LOADING %d, RELAX %.3f) -- %lld of %d columns"
           " (%.2f %%), %lld layers, triggers moist %lld / dry %lld, pairs unstable in T %lld of which held by"
           " composition %lld, worst column %d sweeps of %d, max dT %.4f K, drift h_L %.2e q %.2e, skipped %lld\n",
           TAG, ledoux ? 1 : 0, loading ? 1 : 0, relax, n_columns, m.jm * m.km,
           100.0 * double(n_columns) / double(m.jm * m.km), n_layers, n_moist_trig, n_dry_trig,
           n_unstable_T, n_held, max_passes_used, max_pass, max_dT_K, max_drift_h, max_drift_q, n_skipped);
    if(max_passes_used >= max_pass)
        std::cout << "      " << TAG << ": WARNING - the MCA sweep cap was reached (raise " << TAG << "_MCA_PASSES)" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    printf(" time measured: %.3f seconds for MoistConvectiveAdjustment\n",
           std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() * 1e-9);
    std::cout << "      " << TAG << ": MoistConvectiveAdjustment ended" << std::endl;
}

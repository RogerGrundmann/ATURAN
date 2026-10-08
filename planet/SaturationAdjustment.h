/*
 * SHARED PHYSICS — saturation adjustment and mixed-phase partitioning. MUST BE BYTE-IDENTICAL IN
 * EVERY MODEL THAT USES IT; `make check-shared` verifies that against planet/SHARED.md5.
 *
 * Reference algorithm:
 *   Tao, W.-K., Simpson, J., and McCumber, M.:
 *   "An Ice-Water Saturation Adjustment", AMS Notes and Correspondence, 1988.
 *
 * ===== WHY THIS IS SHARED =====
 *
 * SaturationAdjustmentJup.cpp and SaturationAdjustmentSat.cpp were the same routine written
 * twice: the same triple loop, the same entry guards, the same Tao mixed-phase iteration to the
 * same budget and tolerance, the same clamps, the same write-back and the same report. Comparing
 * them line by line left exactly two real differences, and both are model FACTS rather than
 * algorithm variants, so both are asked of the model instead of being written into two copies:
 *
 *   m.is_solid(i,j,k)                cells with no thermodynamic state. ATJUP's obstacle interior
 *                                    is SeaMount.x == 1.0, where bcSolidGround sets t = p_stat = 0;
 *                                    ATSAT has no solid body and returns false, so the test folds
 *                                    away. This accessor already existed for exactly this purpose —
 *                                    it is what let Turbulence and PressureSolver be shared.
 *
 *   Planet::satadj_updates_pstat()   whether the hydrostatic pressure is rewritten from the
 *                                    adjusted temperature at the end of each cell. ATSAT does,
 *                                    ATJUP does not. This is a real modelling disagreement about
 *                                    where p_stat is allowed to respond to latent heating, it was
 *                                    already flagged-not-settled in ATSAT's header, and sharing
 *                                    the file is NOT the moment to settle it. The hook keeps both
 *                                    behaviours exactly as they were.
 *
 * ===== WHAT IS NOT A PARAMETER ANY MORE =====
 *
 * ATJUP's run() took coeff_A, coeff_B, coeff_A_i, coeff_B_i, cp_gas, r_gas and m_mol and used
 * none of them — seven parameters that appeared once each, in the signature. They are gone. The
 * two-coefficient Clausius-Clapeyron pair they belong to is used by Thermo_Jup.cpp, which calls
 * ATPhys::clausius_clapeyron directly.
 *
 * ===== THE ICE COEFFICIENTS ARE STILL THE CALLER'S PROBLEM =====
 *
 * The routine takes a full ice quadruple (C_i, L0_i, del_alf_i, del_bet_i) and uses it for the
 * deposition/sublimation branch, as it always did. ATJUP passes real ice values; ATSAT passes its
 * LIQUID pair, because its parameter set has no ice pair and inventing numbers for Saturn's H2O,
 * NH3 and CH4 ices is a physics decision rather than a port. That substitution lives at ATSAT's
 * call site, not here, so supplying real values later changes one place.
 *
 * ===== <TAG>_SATADJ_DIAG AND <TAG>_SATADJ_CONSERVE (2026-10-08), BOTH DEFAULT 0 = OFF =====
 *
 * Mirrored from ATOM_Precipitation's ATM_SATADJ_DIAG / ATM_SATADJ_PHASE, where the same routine
 * was charged +5.5e5 mm/a of column water by an independent budget and ONE clip turned out to be
 * 100 % of it. The iteration conserves total water in its RATES (d_q_c + d_q_i = -d_q_v); what
 * is WRITTEN passes through five places that do not:
 *
 *   entry ice delete    if(t_u > t_0) ice = 0            ice above the melting point is deleted,
 *                                                        not melted: no liquid, no latent heat
 *   entry negatives     c / cloud / ice < 0 -> 0         a source of whatever was negative
 *   loop clips          q_v_b / q_c_b / q_i_b < 0 -> 0   the increment is split by TEMPERATURE
 *                                                        (CND / DEP), so an evaporation step takes
 *                                                        from a phase that may be empty; the vapour
 *                                                        gains it in full and nothing loses it
 *   loop ice delete     if(T >= t_0) q_i_b = 0           the same deletion, inside the iteration
 *   final clamp         <= 0 -> 0 over the whole field
 *
 * <TAG>_SATADJ_DIAG=1 is PRINT-ONLY: per call (per gas) it reports the column water
 * (vapour + cloud + ice, sin(colatitude)- and layer-thickness-weighted mean, g/m2) before and
 * after, and charges the change to those five buckets. after - before - sum(buckets) is printed
 * as `unattributed` and is an IDENTITY (round-off) — if it is not, a write was missed. The top
 * level carries no layer thickness and is not counted.
 *
 * <TAG>_SATADJ_CONSERVE=1 is the repair of the two that are defects rather than safety nets:
 *   - the loop increment takes from the phase that is there: what the temperature split would
 *     take from an empty phase is taken from the other one, and what neither holds is left in
 *     the vapour (d_q_v is reduced); condensation cannot take more vapour than the cell holds.
 *     The latent-heat line consumes the corrected pair, so the energy follows the mass.
 *   - ice above the melting point MELTS into the liquid and cools the cell by the latent heat of
 *     fusion (ls - lv), limited to the amount that brings the cell to t_0 — at entry (every
 *     fluid cell, the skipped subsaturated ones too) and inside the iteration. With ls <= lv
 *     (a planet whose ice pair is its liquid pair) it moves the mass and no heat.
 * The entry negatives and the final clamp are left as they are; the instrument sizes them.
 * With both knobs unset every written field is byte-identical to the routine before them.
 */

#pragma once

#include "ATPhys.h"   // saturation_vapour_pressure

#include <cmath>
#include <chrono>
#include <cstdio>
#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class Array;


template<class Planet>
class SaturationAdjustment {
public:
    explicit SaturationAdjustment(Planet& model) : m(model) {}

    void run(const std::string& gas,
             double t_0,       double t_00,
             double ep,        double lv,   double ls,
             double C,         double L0,   double R,
             double del_alf,   double del_bet,
             double C_i,       double L0_i,
             double del_alf_i, double del_bet_i,
             Array& c,         Array& cloud, Array& ice);

private:
    Planet& m;

    // ATJUP's budget and tolerance. ATSAT's inherited routine used 30 passes to 1e-4 and the
    // mirror adopted these when it was written, so the two agree and this is one constant.
    static constexpr int    iter_prec_end = 15;
    static constexpr double q_diff_min    = 1.0e-3;
};


template<class Planet>
void SaturationAdjustment<Planet>::run(
        const std::string& gas,
        double t_0,       double t_00,
        double ep,        double lv,   double ls,
        double C,         double L0,   double R,
        double del_alf,   double del_bet,
        double C_i,       double L0_i,
        double del_alf_i, double del_bet_i,
        Array& c,         Array& cloud, Array& ice)
{
    std::cout << std::endl << "      SaturationAdjustment of " << gas << " begin" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double exp_pressure = m.g / (m.gam * m.R_ref);
    const double t_range_inv  = 1.0 / (t_0 - t_00);

    // <TAG>_SATADJ_DIAG / <TAG>_SATADJ_CONSERVE, see the note at the top. Read once per planet.
    static const bool diag = [](){
        const char* e = ATPhys::env_for(Planet::planet_tag(), "SATADJ_DIAG"); return e && atoi(e) != 0; }();
    static const bool conserve = [](){
        const char* e = ATPhys::env_for(Planet::planet_tag(), "SATADJ_CONSERVE"); return e && atoi(e) != 0; }();
    const double Lf = ls - lv;                       // latent heat of fusion, same unit as lv / ls
    // Budget sums: weight sin(colatitude) * layer thickness [m] on densities [kg/m3] -> kg/m2.
    double w_before = 0.0, w_after = 0.0, w_norm = 0.0;
    double b_entry_ice = 0.0, b_entry_neg = 0.0, b_loop_clip = 0.0, b_loop_ice = 0.0, b_final = 0.0;
    long   n_entry_ice = 0, n_loop_clip = 0, n_loop_ice = 0, n_noconv = 0, n_adjusted = 0;
    if(diag){
        for(int j = 0; j < m.jm; j++) w_norm += std::sin(m.the.z[j]) * m.km;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:w_before)
        for(int k = 0; k < m.km; k++){
            for(int j = 0; j < m.jm; j++){
                const double sj = std::sin(m.the.z[j]);
                for(int i = 0; i < m.im; i++){
                    const double dz = m.layer_thickness_m(i);
                    if(dz > 0.0) w_before += sj * dz * (c.x[i][j][k] + cloud.x[i][j][k] + ice.x[i][j][k]);
                }
            }
        }
    }

    // Shared diagnostic state — written under omp critical, so which cell gets reported still
    // depends on thread arrival, but the writes cannot tear. Nothing here feeds the solution.
    bool   sat_found       = false;
    int    iter_prec_found = 0;
    int    i_sat = 0, j_sat = 0, k_sat = 0;
    double height_sat = 0.0;
    double t_latent = 0.0, p_latent = 0.0;
    double t_sat    = 0.0, p_sat    = 0.0;
    double t_u_sat  = 0.0, p_u_sat  = 0.0;
    double q_v_b_sat = 0.0, q_c_b_sat = 0.0, q_i_b_sat = 0.0;
    double saturation = 0.0;
    // Which cell gets reported is decided by POSITION and not by thread arrival order; -1 is
    // "no cell found yet", and every real key is >= 0. See the critical block below.
    long long sat_key = -1;

    // -----------------------------------------------------------------------
    // Main saturation-adjustment loop — fully independent per cell
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:b_entry_ice,b_entry_neg,b_loop_clip,b_loop_ice,n_entry_ice,n_loop_clip,n_loop_ice,n_noconv,n_adjusted)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){

                double t_u = m.t.x[i][j][k] * m.t_ref;         // [K]; changed only by the CONSERVE melt
                const double p_u = m.p_stat.x[i][j][k];        // [bar]

                // Cells with no usable thermodynamic state are skipped: a solid interior, where
                // t and p_stat are zero, and any cell that has lost a positive temperature or
                // pressure. Both formulas below are singular there —
                // ATPhys::saturation_vapour_pressure forms -L0/T (-inf at T=0) and
                // del_alf*log(T) (0*-inf = NaN), and q_Rain_0 divides by p_u — and the entry test
                // further down cannot catch a NaN, because every comparison against one is false.
                // Written as !(x > 0.0) so a NaN that arrived from elsewhere is skipped rather
                // than propagated; on ATJUP this was one of the two seeds of a domain-wide NaN.
                if(!(t_u > 0.0) || !(p_u > 0.0) || m.is_solid(i, j, k)) continue;

                // Budget weight of this cell (0 where the level has no thickness, i.e. the top).
                double wgt = 0.0;
                if(diag){ const double dz = m.layer_thickness_m(i); if(dz > 0.0) wgt = std::sin(m.the.z[j]) * dz; }

                // Enforce physical bounds
                if(t_u > t_0){
                    if(conserve){
                        // Melt, do not delete: into the liquid, cooled by the latent heat of fusion,
                        // and only as much as brings the cell down to t_0.
                        if(ice.x[i][j][k] > 0.0){
                            double dm = ice.x[i][j][k];
                            if(Lf > 0.0){
                                const double hc  = m.cp_mix * m.rho_at(i, j, k);
                                const double cap = hc * (t_u - t_0) / Lf;
                                if(dm > cap) dm = cap;
                                t_u -= Lf * dm / hc;
                                m.t.x[i][j][k] = t_u / m.t_ref;
                            }
                            ice.x[i][j][k]   -= dm;
                            cloud.x[i][j][k] += dm;
                        }
                    } else {
                        if(ice.x[i][j][k] != 0.0){ b_entry_ice -= wgt * ice.x[i][j][k]; n_entry_ice++; }
                        ice.x[i][j][k] = 0.0;
                    }
                }
                if(c.x[i][j][k]     < 0.0){ b_entry_neg -= wgt * c.x[i][j][k];     c.x[i][j][k]     = 0.0; }
                if(cloud.x[i][j][k] < 0.0){ b_entry_neg -= wgt * cloud.x[i][j][k]; cloud.x[i][j][k] = 0.0; }
                if(ice.x[i][j][k]   < 0.0){ b_entry_neg -= wgt * ice.x[i][j][k];   ice.x[i][j][k]   = 0.0; }

                // Density for both the saturation quantity and the latent-heat divisor. rho_at()
                // is each model's single gated accessor — the constant reference density r_mix
                // unless <TAG>_LOCAL_RHO is set — so the microphysics and the adjustment saturate
                // on the same density by construction.
                //
                // Why the gate defaults to the reference density: on ATSAT, switching it to the
                // local one diverged the run, T to 4.65e5 degC in eight iterations. Both halves of
                // the feedback push the same way — the divisor cp_mix*rho is small in a thin cell
                // so the heating per unit condensate is large, and the target rho*ep*E/p is scaled
                // by the same small rho so more vapour reads as supersaturated.
                const double rho_c = m.rho_at(i, j, k);

                const double E_Rain_0 = ATPhys::saturation_vapour_pressure(t_u, C, L0, R,
                                                                           del_alf, del_bet);
                const double q_Rain_0 = rho_c * ep * E_Rain_0 / p_u;

                // Skip: subsaturated, already saturated, or the SVP underflowed to zero in a
                // very cold cell.
                if(c.x[i][j][k] <= q_Rain_0 || q_Rain_0 <= 0.0) continue;

                // ---- Mixed-phase iteration (Tao et al. 1988) ----
                double q_v_b   = c.x[i][j][k];
                double q_c_b   = cloud.x[i][j][k];
                double q_i_b   = ice.x[i][j][k];
                double T       = t_u;
                double q_v_hyp = q_v_b;

                bool cell_found = false;
                int  cell_iter  = 0;
                double cell_clip = 0.0, cell_ice = 0.0;      // water the clips created / the delete removed

                for(int itr = 1; itr <= iter_prec_end; itr++){
                    double CND = (T - t_00) * t_range_inv;
                    double DEP = (t_0 - T)  * t_range_inv;
                    if(T <= t_00){ CND = 0.0; DEP = 1.0; }
                    if(T >= t_0) { CND = 1.0; DEP = 0.0; }

                    double d_q_v = q_v_hyp - q_v_b;
                    double d_q_c = -d_q_v * CND;
                    double d_q_i = -d_q_v * DEP;
                    if(conserve){
                        // Condensation cannot take more vapour than the cell holds ...
                        if(d_q_v < -q_v_b){ d_q_v = -q_v_b; d_q_c = -d_q_v * CND; d_q_i = -d_q_v * DEP; }
                        // ... and evaporation takes from the phase that is there: the part the
                        // temperature split asks of an empty phase goes to the other one, and what
                        // neither holds stays in the vapour.
                        if(q_c_b + d_q_c < 0.0){ const double ex = -(q_c_b + d_q_c); d_q_c = -q_c_b; d_q_i -= ex; }
                        if(q_i_b + d_q_i < 0.0){ const double ex = -(q_i_b + d_q_i); d_q_i = -q_i_b; d_q_c -= ex;
                                                 if(q_c_b + d_q_c < 0.0) d_q_c = -q_c_b; }
                        d_q_v = -(d_q_c + d_q_i);
                    }

                    T     += (lv * d_q_c + ls * d_q_i) / (m.cp_mix * rho_c);
                    q_v_b += d_q_v;
                    q_c_b += d_q_c;
                    q_i_b += d_q_i;

                    if(q_v_b < 0.0){ cell_clip -= q_v_b; q_v_b = 0.0; }
                    if(q_c_b < 0.0){ cell_clip -= q_c_b; q_c_b = 0.0; }
                    if(q_i_b < 0.0){ cell_clip -= q_i_b; q_i_b = 0.0; }

                    // AT THE UPDATED TEMPERATURE T, not at the entry temperature. This is the
                    // feedback the scheme exists to resolve: latent heat raises T, which raises
                    // the saturation vapour pressure, which limits further condensation.
                    const double E_Rain = ATPhys::saturation_vapour_pressure(T, C,   L0,   R,
                                                                             del_alf,   del_bet);
                    const double E_Ice  = ATPhys::saturation_vapour_pressure(T, C_i, L0_i, R,
                                                                             del_alf_i, del_bet_i);
                    const double q_Rain = rho_c * ep * E_Rain / p_u;
                    const double q_Ice  = rho_c * ep * E_Ice  / p_u;

                    if(q_c_b > 0.0 && q_i_b > 0.0)
                        q_v_hyp = (q_c_b * q_Rain + q_i_b * q_Ice) / (q_c_b + q_i_b);
                    else if(q_i_b == 0.0) q_v_hyp = q_Rain;
                    else                  q_v_hyp = q_Ice;

                    if(T >= t_0){
                        if(conserve){
                            if(q_i_b > 0.0){                 // melt, limited to what brings T to t_0
                                double dm = q_i_b;
                                if(Lf > 0.0){
                                    const double hc  = m.cp_mix * rho_c;
                                    const double cap = hc * (T - t_0) / Lf;
                                    if(dm > cap) dm = cap;
                                    T -= Lf * dm / hc;
                                }
                                q_i_b -= dm;
                                q_c_b += dm;
                            }
                        } else {
                            cell_ice += q_i_b;
                            q_i_b = 0.0;
                        }
                    }

                    const double q_diff = std::fabs(q_v_b - q_v_hyp) / (q_v_hyp + 1e-20);
                    if(q_diff <= q_diff_min){
                        cell_found = true;
                        cell_iter  = itr;
                        break;
                    }
                    q_v_hyp = 0.5 * (q_v_hyp + q_v_b);   // has smoothing effect
                }

                n_adjusted++;
                if(!cell_found) n_noconv++;
                if(cell_clip != 0.0){ b_loop_clip += wgt * cell_clip; n_loop_clip++; }
                if(cell_ice  != 0.0){ b_loop_ice  -= wgt * cell_ice;  n_loop_ice++;  }

                // Write converged cell values back
                c.x[i][j][k]     = q_v_b;
                cloud.x[i][j][k] = q_c_b;
                ice.x[i][j][k]   = q_i_b;
                m.t.x[i][j][k]   = T / m.t_ref;

                // Hydrostatic pressure rebuilt from the adjusted temperature — ATSAT only. See
                // the note at the top: this is a flagged modelling disagreement, not a variant of
                // the algorithm, and the hook preserves each model's existing behaviour.
                if(Planet::satadj_updates_pstat()){
                    m.p_stat.x[i][j][k] = m.p_ref
                        * std::pow(m.t.x[i][j][k], exp_pressure);   // [bar]
                }

                if(cell_found){
                    // ===== THE WINNER IS CHOSEN BY POSITION, NOT BY ARRIVAL ORDER =====
                    //
                    // This block used to assign unconditionally, so the cell reported was
                    // whichever thread happened to enter the critical section LAST. That is not
                    // a race — the section is properly synchronised and writes reporting
                    // variables only, and every output file is bit-identical across it — but it
                    // made the LOG non-comparable above one thread: two runs of the same binary
                    // at the same thread count named different cells.
                    //
                    // The loop nest is k, then j, then i, so serial traversal visits keys in
                    // increasing order and "the last cell found wins" is exactly "the largest
                    // key wins". Taking the maximum therefore reproduces the single-threaded
                    // answer at ANY thread count: this is a determinism fix, not a change of
                    // which cell is meant. The key is built in the nest's own order for that
                    // reason — (i,j,k) would be just as deterministic and would pick a DIFFERENT
                    // cell, silently changing every reported diagnostic.
                    //
                    // ATURAN and ATNEPT carry the same block in their own inherited routines
                    // (SaturationAdjustment{Uran,Nept}.cpp, the default there) and took this
                    // same fix first; this is the shared copy that ATSAT and ATJUP run.
                    const long long key = ((long long)k * m.jm + j) * m.im + i;
                    #pragma omp critical
                    if(key > sat_key)
                    {
                        sat_key         = key;
                        sat_found       = true;
                        iter_prec_found = cell_iter;
                        i_sat = i; j_sat = j; k_sat = k;
                        height_sat = m.get_layer_height(i_sat);
                        t_u_sat    = t_u;
                        p_u_sat    = p_u;
                        t_sat      = T;
                        p_sat      = m.p_ref * std::pow(t_sat / m.t_ref, exp_pressure);
                        t_latent   = t_sat - t_u_sat;
                        p_latent   = p_sat - p_u_sat;
                        q_v_b_sat  = q_v_b;
                        q_c_b_sat  = q_c_b;
                        q_i_b_sat  = q_i_b;
                        saturation = q_v_b - q_Rain_0;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Clamp negative values left by the iteration
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:b_final,w_after)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){
                double wgt = 0.0;
                if(diag){ const double dz = m.layer_thickness_m(i); if(dz > 0.0) wgt = std::sin(m.the.z[j]) * dz; }
                if(c.x[i][j][k]     <= 0.0){ b_final -= wgt * c.x[i][j][k];     c.x[i][j][k]     = 0.0; }
                if(cloud.x[i][j][k] <= 0.0){ b_final -= wgt * cloud.x[i][j][k]; cloud.x[i][j][k] = 0.0; }
                if(ice.x[i][j][k]   <= 0.0){ b_final -= wgt * ice.x[i][j][k];   ice.x[i][j][k]   = 0.0; }
                w_after += wgt * (c.x[i][j][k] + cloud.x[i][j][k] + ice.x[i][j][k]);
            }
        }
    }

    if(diag && w_norm > 0.0){
        const double f = 1.0e3 / w_norm;             // kg/m2 summed -> g/m2 mean
        const double sum_b = b_entry_ice + b_entry_neg + b_loop_clip + b_loop_ice + b_final;
        printf("      %s: [SATADJ-DIAG] %s  column water (vapour+cloud+ice) [g/m2]: before %.9e  after %.9e  change %+.6e   (CONSERVE=%d)\n",
               Planet::planet_tag(), gas.c_str(), w_before * f, w_after * f, (w_after - w_before) * f, conserve ? 1 : 0);
        printf("      %s: [SATADJ-DIAG] %s    entry ice deleted above t_0  %+.6e  (%ld cells)\n", Planet::planet_tag(), gas.c_str(), b_entry_ice * f, n_entry_ice);
        printf("      %s: [SATADJ-DIAG] %s    entry negatives clipped      %+.6e\n", Planet::planet_tag(), gas.c_str(), b_entry_neg * f);
        printf("      %s: [SATADJ-DIAG] %s    loop clips (phase split)     %+.6e  (%ld cells)\n", Planet::planet_tag(), gas.c_str(), b_loop_clip * f, n_loop_clip);
        printf("      %s: [SATADJ-DIAG] %s    loop ice deleted at T >= t_0 %+.6e  (%ld cells)\n", Planet::planet_tag(), gas.c_str(), b_loop_ice * f, n_loop_ice);
        printf("      %s: [SATADJ-DIAG] %s    final clamp                  %+.6e\n", Planet::planet_tag(), gas.c_str(), b_final * f);
        printf("      %s: [SATADJ-DIAG] %s    unattributed                 %+.6e   cells adjusted %ld, not converged in %d passes %ld\n",
               Planet::planet_tag(), gas.c_str(), (w_after - w_before - sum_b) * f, n_adjusted, iter_prec_end, n_noconv);
    }

    // -----------------------------------------------------------------------
    // Timing and report
    // -----------------------------------------------------------------------
    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

    if(!sat_found){
        std::cout << "      NO saturation in SaturationAdjustment of " << gas << " found"
                  << std::endl
                  << "      iter_prec_end = " << iter_prec_end << std::endl;
    } else {
        std::cout << "      saturation in SaturationAdjustment of " << gas << " found" << std::endl
             << "      iter_prec_found = " << iter_prec_found
             << "   iter_prec_end = "      << iter_prec_end   << std::endl
             << "      i_sat = "     << i_sat
             << "   j_sat = "        << j_sat
             << "   k_sat = "        << k_sat
             << "   height_sat[km] = " << height_sat  << std::endl
             << "      p_stat[bar] = "  << p_sat
             << "   p_u[bar] = "        << p_u_sat
             << "   p_latent[bar] = "   << p_latent   << std::endl
             << "      T[°C] = "        << t_sat    - m.t_ref
             << "   t_u[°C] = "         << t_u_sat  - m.t_ref
             << "   t_latent[°C] = "    << t_latent  << std::endl
             << "      saturation[g/m³] = " << saturation * 1e3 << std::endl
             << "      " << gas << " humid[g/m³] = "  << q_v_b_sat * 1e3
             << "   cloud[g/m³] = "  << q_c_b_sat * 1e3
             << "   ice[g/m³] = "    << q_i_b_sat * 1e3 << std::endl;
    }

    std::cout << "      SaturationAdjustment of " << gas << " ended" << std::endl;
}

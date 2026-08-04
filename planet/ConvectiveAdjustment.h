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
 * here, because each is a modelling decision with consequences for the microphysics that already
 * runs in the saturation adjustment and the precipitation scheme.
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

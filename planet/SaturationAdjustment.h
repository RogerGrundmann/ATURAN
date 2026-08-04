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

    // -----------------------------------------------------------------------
    // Main saturation-adjustment loop — fully independent per cell
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){

                const double t_u = m.t.x[i][j][k] * m.t_ref;   // [K]
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

                // Enforce physical bounds
                if(t_u > t_0)              ice.x[i][j][k]   = 0.0;
                if(c.x[i][j][k]     < 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] < 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   < 0.0) ice.x[i][j][k]   = 0.0;

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

                for(int itr = 1; itr <= iter_prec_end; itr++){
                    double CND = (T - t_00) * t_range_inv;
                    double DEP = (t_0 - T)  * t_range_inv;
                    if(T <= t_00){ CND = 0.0; DEP = 1.0; }
                    if(T >= t_0) { CND = 1.0; DEP = 0.0; }

                    const double d_q_v = q_v_hyp - q_v_b;
                    const double d_q_c = -d_q_v * CND;
                    const double d_q_i = -d_q_v * DEP;

                    T     += (lv * d_q_c + ls * d_q_i) / (m.cp_mix * rho_c);
                    q_v_b += d_q_v;
                    q_c_b += d_q_c;
                    q_i_b += d_q_i;

                    if(q_v_b < 0.0) q_v_b = 0.0;
                    if(q_c_b < 0.0) q_c_b = 0.0;
                    if(q_i_b < 0.0) q_i_b = 0.0;

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

                    if(T >= t_0) q_i_b = 0.0;

                    const double q_diff = std::fabs(q_v_b - q_v_hyp) / (q_v_hyp + 1e-20);
                    if(q_diff <= q_diff_min){
                        cell_found = true;
                        cell_iter  = itr;
                        break;
                    }
                    q_v_hyp = 0.5 * (q_v_hyp + q_v_b);   // has smoothing effect
                }

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
                    #pragma omp critical
                    {
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
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){
                if(c.x[i][j][k]     <= 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] <= 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   <= 0.0) ice.x[i][j][k]   = 0.0;
            }
        }
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

#include "SaturationAdjustmentUran.h"
#include "cUranusModel.h"
#include "Utils.h"
#include "Array.h"

#include <iostream>
#include <cmath>
#include <cstdio>

using namespace std;
using namespace AtomUtils;

void SaturationAdjustmentUran::run_legacy(const std::string& gas,
    double coeff_A,   double coeff_B,
    double coeff_A_i, double coeff_B_i,
    double t_0,       double t_00,
    double ep,        double lv,  double ls,
    double cp,        double r,
    double C,         double L0,  double R,
    double del_alf,   double del_bet,   double m_mol,
    Array& c,         Array& cloud,   Array& ice)
{
    const int im = m.im, jm = m.jm, km = m.km;

    cout << endl << "      SaturationAdjustment of " << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    cout.precision(9);

    const double exp_pressure = m.g / (m.gam * m.R_ref);
    const double t_range_inv  = 1.0 / (t_0 - t_00);

    // Shared diagnostic state — protected by omp critical where written
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

    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){

                const double t_u = m.t.x[i][j][k] * m.t_ref;
                const double p_u = m.p_stat.x[i][j][k];

                if(t_u > t_0)              ice.x[i][j][k]   = 0.0;
                if(c.x[i][j][k]     < 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] < 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   < 0.0) ice.x[i][j][k]   = 0.0;

                const double E_Rain_0 = saturation_vapour_pressure(t_u, C, L0, R, del_alf, del_bet);
                const double q_Rain_0 = m.r_mix * ep * E_Rain_0 / p_u;

                // skip subsaturated cells and cells where SVP underflowed to 0
                if(c.x[i][j][k] <= q_Rain_0 || q_Rain_0 <= 0.0) continue;

                double q_v_b = c.x[i][j][k];
                double q_c_b = cloud.x[i][j][k];
                double q_i_b = ice.x[i][j][k];
                double T     = t_u;

                bool cell_found = false;
                int  cell_iter  = 0;

                for(int itr = 1; itr <= iter_prec_end; itr++){

                    double CND = (T - t_00) * t_range_inv;
                    double DEP = (t_0 - T)  * t_range_inv;
                    if(T <= t_00){ CND = 0.0; DEP = 1.0; }
                    if(T >= t_0) { CND = 1.0; DEP = 0.0; }

                    const double E_Rain = saturation_vapour_pressure(T, C, L0, R, del_alf, del_bet);
                    const double q_Rain = m.r_mix * ep * E_Rain / p_u;

                    // convergence check before applying correction
                    const double q_diff = std::fabs(q_v_b - q_Rain) / (q_Rain + 1.0e-20);
                    if(q_diff <= q_diff_min){
                        cell_found = true;
                        cell_iter  = itr;
                        break;
                    }

                    // implicit Tao correction — denominator accounts for saturation feedback
                    // dq_sat/dT = q_Rain * f'(T), f'(T) = (L0/T² + del_alf/T + del_bet)/(1e-3*R)
                    const double L_mix     = lv * CND + ls * DEP;
                    const double dqsat_dT  = q_Rain * (L0/(T*T) + del_alf/T + del_bet) / (1.0e-3 * R);
                    const double denom     = 1.0 + L_mix * dqsat_dT / (m.cp_mix * m.r_mix);
                    const double d_q_total = (q_v_b - q_Rain) / denom;

                    T     += L_mix * d_q_total / (m.cp_mix * m.r_mix);
                    q_v_b -= d_q_total;
                    q_c_b += CND * d_q_total;
                    q_i_b += DEP * d_q_total;

                    if(q_v_b < 0.0) q_v_b = 0.0;
                    if(q_c_b < 0.0) q_c_b = 0.0;
                    if(q_i_b < 0.0) q_i_b = 0.0;
                    if(T >= t_0)    q_i_b = 0.0;

                } // itr

                // write back best approximation regardless of convergence
                c.x[i][j][k]        = q_v_b;
                cloud.x[i][j][k]    = q_c_b;
                ice.x[i][j][k]      = q_i_b;
                m.t.x[i][j][k]      = T / m.t_ref;

                if(cell_found){
                    // ===== THE WINNER IS CHOSEN BY POSITION, NOT BY ARRIVAL ORDER =====
                    //
                    // This block used to assign unconditionally, so the cell reported was
                    // whichever thread happened to enter the critical section LAST. That is not
                    // a race — the section is properly synchronised and writes reporting
                    // variables only, and every output file was bit-identical across it — but it
                    // made the LOG non-comparable above one thread: two 16-thread runs of the
                    // same binary named different cells, and those were the only log lines that
                    // differed once the pressure-solver race was fixed.
                    //
                    // The loop nest is k, then j, then i, so serial traversal visits keys in
                    // increasing order and "the last cell found wins" is exactly "the largest
                    // key wins". Taking the maximum therefore reproduces the single-threaded
                    // answer at ANY thread count, which is the whole point: this is a
                    // determinism fix, not a change of which cell is meant. The key is built in
                    // the nest's own order for that reason — (i,j,k) would pick a different cell
                    // and change the reported answer.
                    const long long key = ((long long)k * jm + j) * im + i;
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

            } // end i
        } // end j
    } // end k

    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k]     <= 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] <= 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   <= 0.0) ice.x[i][j][k]   = 0.0;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

    if(!sat_found)
        cout << "      NO saturation found in SaturationAdjustment of " << gas << endl
             << "      iter_prec_end = " << iter_prec_end << endl;
    else
        cout << "      saturation of water vapour in SaturationAdjustment of " << gas << " found" << endl
             << "      iter_prec_found = " << iter_prec_found
             << "   iter_prec_end = "      << iter_prec_end   << endl
             << "      i_sat = "     << i_sat
             << "   j_sat = "        << j_sat
             << "   k_sat = "        << k_sat
             << "   height_sat[km] = " << height_sat  << endl
             << "      p_stat[bar] = "  << p_sat
             << "   p_u[bar] = "        << p_u_sat
             << "   p_latent[bar] = "   << p_latent   << endl
             << "      T[°C] = "        << t_sat    - m.t_ref
             << "   t_u[°C] = "         << t_u_sat  - m.t_ref
             << "   t_latent[°C] = "    << t_latent  << endl
             << "      saturation[g/m³] = " << saturation * 1e3 << endl
             << "      " << gas << " humid[g/m³] = "  << q_v_b_sat * 1e3
             << "   cloud[g/m³] = "  << q_c_b_sat * 1e3
             << "   ice[g/m³] = "    << q_i_b_sat * 1e3 << endl;

    cout << "      SaturationAdjustment of " << gas << " ended" << endl;
    return;
}

/*
 * Dispatch between the inherited routine above and the SHARED SaturationAdjustment<Planet>.
 * Default is the inherited one; ATURAN_SATADJ=1 selects the shared algorithm. See the header for
 * what is known about the difference (measured on ATSAT, not on Uranus) and for why the ice
 * quadruple below is the liquid one.
 */
void SaturationAdjustmentUran::run(const std::string& gas,
    double coeff_A,   double coeff_B,
    double coeff_A_i, double coeff_B_i,
    double t_0,       double t_00,
    double ep,        double lv,  double ls,
    double cp,        double r,
    double C,         double L0,  double R,
    double del_alf,   double del_bet,   double m_mol,
    Array& c,         Array& cloud,   Array& ice)
{
    if(mirrored_enabled() == 0){
        run_legacy(gas, coeff_A, coeff_B, coeff_A_i, coeff_B_i, t_0, t_00, ep, lv, ls,
                   cp, r, C, L0, R, del_alf, del_bet, m_mol, c, cloud, ice);
        return;
    }

    // THE ICE QUADRUPLE IS THE LIQUID ONE — see the header. This is the single place to change
    // when real Uranus ice coefficients exist; the shared algorithm already treats the liquid
    // and ice pairs separately, so nothing else has to move.
    SaturationAdjustment<cUranusModel>(m).run(gas, t_0, t_00, ep, lv, ls,
                                              C, L0, R, del_alf, del_bet,
                                              C, L0,    del_alf, del_bet,
                                              c, cloud, ice);
}

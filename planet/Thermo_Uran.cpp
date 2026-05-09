/*
 * Atmosphere General Circulation Modell(ATURAN) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cUranusModel.h"
#include "SaturationAdjustmentUran.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


void cUranusModel::TropopauseLocation(){
//    cout << endl << "      ATURAN: TropopauseLocation" << endl;

    im_tropopause = std::vector<int>(jm, 0);
    const int    j_half      = (jm - 1) / 2;
    const double d_j_half    = (double)j_half;
    const double trop_u2_eff = (double)(i_beg - i_max);

    #pragma omp parallel for schedule(static)
    for(int j = 0; j < jm; j++){
        const double d_j = (double)j;
        im_tropopause[j] = (int)((trop_u2_eff
            * (d_j*d_j / (d_j_half*d_j_half) - 2.0*d_j / d_j_half))
            + (double)i_beg);
    }

//    cout << "      ATURAN: TropopauseLocation ended" << endl;
}
/*
*
*/
void cUranusModel::Latent_Heat(){
    cout << endl << "      ATURAN: Latent_Heat" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(static) collapse(2)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Q_Latent.x[0][j][k]   = 0.0;
            Q_Sensible.x[0][j][k] = 0.0;
        }
    }

    const double inv_2dr        = 1.0 / (2.0 * dr);
    const double inv_L_atm2     = 1.0 / (L_atm * L_atm);

    #pragma omp parallel for schedule(static)
    for(int j = 1; j < jm-1; j++){
        const double sinthe = sin(the.z[j]);

        for(int k = 1; k < km-1; k++){
            for(int i = im-2; i >= 1; i--){
                const double rm       = rad.z[i];
                const double rmsinthe = rm * sinthe;
                const double inv_2rm_dthe      = 1.0 / (2.0 * rm * dthe);
                const double inv_2rmsinthe_dphi = 1.0 / (2.0 * rmsinthe * dphi);

                const double t_u = t.x[i][j][k] * t_ref;
                const double p_u = p_stat.x[i][j][k];

                const double E_Rain     = 1e3 * SaturationAdjustmentUran::clausius_clapeyron(t_u, coeff_h2o_A,   coeff_h2o_B);
                const double E_Ice      = 1e3 * SaturationAdjustmentUran::clausius_clapeyron(t_u, coeff_h2o_A_i, coeff_h2o_B_i);
                const double q_Rain     = ep_h2o * E_Rain     / (p_u - E_Rain);
                const double q_Ice      = ep_h2o * E_Ice      / (p_u - E_Ice);

                const double E_Rain_h2s = 1e3 * SaturationAdjustmentUran::clausius_clapeyron(t_u, coeff_h2s_A,   coeff_h2s_B);
                const double E_Ice_h2s  = 1e3 * SaturationAdjustmentUran::clausius_clapeyron(t_u, coeff_h2s_A_i, coeff_h2s_B_i);
                const double q_Rain_h2s = ep_h2s * E_Rain_h2s / (p_u - E_Rain_h2s);
                const double q_Ice_h2s  = ep_h2s * E_Ice_h2s  / (p_u - E_Ice_h2s);

                const double E_Rain_nh3 = 1e3 * SaturationAdjustmentUran::clausius_clapeyron(t_u, coeff_nh3_A,   coeff_nh3_B);
                const double E_Ice_nh3  = 1e3 * SaturationAdjustmentUran::clausius_clapeyron(t_u, coeff_nh3_A_i, coeff_nh3_B_i);
                const double q_Rain_nh3 = ep_nh3 * E_Rain_nh3 / (p_u - E_Rain_nh3);
                const double q_Ice_nh3  = ep_nh3 * E_Ice_nh3  / (p_u - E_Ice_nh3);

                const double u_av = 0.5 * (u.x[i+1][j][k] + u.x[i-1][j][k]);
                const double v_av = 0.5 * (v.x[i+1][j][k] + v.x[i-1][j][k]);
                const double w_av = 0.5 * (w.x[i+1][j][k] + w.x[i-1][j][k]);
                const double velocity_av = -sqrt((u_av*u_av + v_av*v_av + w_av*w_av) / 3.0);

                const double dtdr   = (t.x[i+1][j][k]   - t.x[i-1][j][k])   * inv_2dr;
                const double dtdthe = (t.x[i][j+1][k]   - t.x[i][j-1][k])   * inv_2rm_dthe;
                const double dtdphi = (t.x[i][j][k+1]   - t.x[i][j][k-1])   * inv_2rmsinthe_dphi;

                const double dh2odr   = (h2o.x[i+1][j][k] - h2o.x[i-1][j][k]) * inv_2dr;
                const double dh2odthe = (h2o.x[i][j+1][k] - h2o.x[i][j-1][k]) * inv_2rm_dthe;
                const double dh2odphi = (h2o.x[i][j][k+1] - h2o.x[i][j][k-1]) * inv_2rmsinthe_dphi;

                const double dh2sdr   = (h2s.x[i+1][j][k] - h2s.x[i-1][j][k]) * inv_2dr;
                const double dh2sdthe = (h2s.x[i][j+1][k] - h2s.x[i][j-1][k]) * inv_2rm_dthe;
                const double dh2sdphi = (h2s.x[i][j][k+1] - h2s.x[i][j][k-1]) * inv_2rmsinthe_dphi;

                const double dnh3dr   = (nh3.x[i+1][j][k] - nh3.x[i-1][j][k]) * inv_2dr;
                const double dnh3dthe = (nh3.x[i][j+1][k] - nh3.x[i][j-1][k]) * inv_2rm_dthe;
                const double dnh3dphi = (nh3.x[i][j][k+1] - nh3.x[i][j][k-1]) * inv_2rmsinthe_dphi;

                const double dtemp = dtdr   + dtdthe   + dtdphi;
                const double dh2o  = dh2odr + dh2odthe + dh2odphi;
                const double dh2s  = dh2sdr + dh2sdthe + dh2sdphi;
                const double dnh3  = dnh3dr + dnh3dthe + dnh3dphi;

                // H2O liquid latent heat
                double Q_lat = 0.0;
                if(h2o.x[i][j][k] >= q_Rain)
                    Q_lat = lv_h2o * velocity_av * dh2o * inv_L_atm2;

                // Ice latent heat — declared per-cell to avoid data race under OpenMP
                double Latency_Ice = 0.0;
                if(h2o.x[i][j][k] >= q_Ice)
                    Latency_Ice = ls_h2o * velocity_av * dnh3 * inv_L_atm2;

                // H2S liquid latent heat
                if(h2s.x[i][j][k] >= q_Rain_h2s)
                    Q_lat += lv_h2s * velocity_av * dh2s * inv_L_atm2;

                if(h2s.x[i][j][k] >= q_Ice_h2s)
                    Latency_Ice += ls_h2s * velocity_av * dh2s * inv_L_atm2;

                // NH3 liquid latent heat
                if(nh3.x[i][j][k] >= q_Rain_nh3)
                    Q_lat += lv_nh3 * velocity_av * dnh3 * inv_L_atm2;

                if(nh3.x[i][j][k] >= q_Ice_nh3)
                    Latency_Ice += ls_nh3 * velocity_av * dnh3 * inv_L_atm2;

                Q_Latent.x[i][j][k]   = Q_lat + Latency_Ice;
                Q_Sensible.x[i][j][k] = rho_mix.x[i][j][k] * cp_mix
                    * velocity_av * dtemp * t_ref * inv_L_atm2;

                if(SeaMount.x[i][j][k] == 1.0){
                    Q_Latent.x[i][j][k]   = 0.0;
                    Q_Sensible.x[i][j][k] = 0.0;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Latent_Heat\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: Latent_Heat ended" << endl;
}
/*
*
*/

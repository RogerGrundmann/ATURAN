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


void cUranusModel::init_temperature(){
    cout << endl << "      ATURAN: init_Temperature" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    int j_half = (jm-1)/2;
    double d_j_half = (double)j_half;
    double t_eff = t_pole - t_equator; // non-dimensional
    double d_j = 0.0;
    double height = 0.0;

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            d_j = (double)(j);
            // t_equator and t_pole are the TOP (i=im-1) temperatures [K].
            // T_bottom = T_top + gam * L_atm; gam [K/km], L_atm [km] via get_layer_height.

            const double T_top    = t_eff * (d_j * d_j/(d_j_half * d_j_half)
                - 2.0 * d_j/d_j_half) + t_pole;                         // K at i=im-1
            const double T_bottom = T_top + gam * L_atm;  // K at i=0

            for(int i = 0; i < im; i++){
                height = get_layer_height(i);                           // km
                t.x[i][j][k] = (T_bottom - gam * height) / t_ref;

/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     height[km] = " << height << endl
        << "     t_equator[°C] = " << t_equator - t_ref 
        << "     t_pole[°C] = " << t_pole - t_ref << endl
        << "     gam[K/km] = " << gam
        << "     gam * height[K] = " << gam * height << endl
        << "     t_u[°C] = " << t.x[i][j][k] * t_ref - t_ref
        << "     t_u[K] = " << t.x[i][j][k] * t_ref << endl;
*/

            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_Temperature\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: init_Temperature ended" << endl;
}
/*
*
*/
void cUranusModel::init_PressureDynamic(){
    cout << endl << "      ATNEPT: init_PressureStatic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double coeff = u_0 * u_0 / 6.0 * 1e-5;

    #pragma omp parallel for schedule(static) collapse(2)
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                const double ui = u.x[i][j][k];
                const double vi = v.x[i][j][k];
                const double wi = w.x[i][j][k];
                p_dyn.x[i][j][k] = rho_mix.x[i][j][k] * coeff * (ui*ui + vi*vi + wi*wi);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATNEPT: init_PressureStatic ended" << endl;
}
/*
*
*/
void cUranusModel::init_PressureStatic(){
    cout << endl << "      ATURAN: init_PressureStatic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // P(i,j,k) = p_bottom * (T(i)/T(0))^n,  n = g/(gam*R_ref).
    const double exp_pressure = g / (gam * R_ref);

    #pragma omp parallel for schedule(static) collapse(2)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            const double p_bottom = r_mix * R_mix * t.x[0][j][k] * t_ref * 1e-5;  // highest static pressure
            for(int i = 0; i < im; i++){
                p_stat.x[i][j][k] = p_bottom
                    * pow(t.x[i][j][k] / t.x[0][j][k], exp_pressure);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: init_PressureStatic ended" << endl;
}
/*
*
*/
void cUranusModel::init_vapour_cloud_ice(std::string gas,
    double &c_tropopause, double &coeff_A, double &coeff_B,
    double &coeff_A_i, double &coeff_B_i,
    double &t_0, double &t_00,
    double &ep, double &r, double &m,
    double &C, double &L0, double &R,
    double &del_alf, double &del_bet,
    Array &c, Array &cloud, Array &ice, Array &cloudiness){

    cout << endl << "      ATURAN: init_vapour_cloud_ice of " << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    double r_max_equator     = 0.0;
    double r_max_pole        = 0.0;
    double r_max_add_equator = 0.0;
    double r_max_add_pole    = 0.0;
    double t_add_equator     = 0.0;
    double t_add_pole        = 0.0;
    double magnus            = 1e3;

    if(gas == "CH4"){
        r_max_equator = r_ch4;
        r_max_pole    = 0.7 * r_ch4;
    }
    if(gas == "H2S"){
        r_max_equator = r_h2s;
        r_max_pole    = 0.7 * r_h2s;
    }
    if(gas == "H2O"){
        r_max_equator = r_h2o;
        r_max_pole    = 0.7 * r_h2o;
    }
    if(gas == "NH3"){
        r_max_equator    = r_nh3;
        r_max_pole       = 0.7 * r_nh3;
        r_max_add_equator = r_nh3_add;
        r_max_add_pole   = 0.7 * r_nh3_add;
        t_add_equator    = t_0_nh3;
        t_add_pole       = 0.7 * t_add_equator;
    }

    r_max     = std::vector<double>(jm, r_max_pole);
    r_max_add = std::vector<double>(jm, r_max_add_pole);
    t_add     = std::vector<double>(jm, t_add_pole);

    const double d_j_half     = (double)(jm - 1) / 2.0;
    const double r_max_eff     = r_max_pole     - r_max_equator;
    const double r_max_add_eff = r_max_add_pole - r_max_add_equator;
    const double t_add_eff     = t_add_pole     - t_add_equator;

    for(int j = 0; j < jm; j++){
        const double d_j = (double)j;
        const double par = AtomUtils::parabola(d_j / d_j_half);
        r_max[j]     = r_max_eff     * par + r_max_pole;
        r_max_add[j] = r_max_add_eff * par + r_max_add_pole;
        t_add[j]     = t_add_eff     * par + t_add_pole;
    }

    #pragma omp parallel for schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(is_land(SeaMount, i, j, k)){
                    c.x[i][j][k] = 0.0;
                    continue;
                }
                const double t_u   = t.x[i][j][k] * t_ref;
                const double p_u   = p_stat.x[i][j][k];
                const double E_Rain = SaturationAdjustmentUran::saturation_vapour_pressure(
                    t_u, C, L0, R, del_alf, del_bet);
                const double q_Rain = ep * E_Rain / p_u;

                double val = magnus * r_mix * q_Rain;

                if(val >= r_max[j])            val = r_max[j];
                if((gas == "NH3") && (p_u >= p_00_nh3)) val = r_max[j];
                if((gas == "H2S") && (p_u >= p_00_h2s)) val = r_max[j];
                if((gas == "H2O") && (p_u >= p_0_h2o))  val = r_max[j];
                if((gas == "CH4") && (p_u >= p_0_ch4))  val = r_max[j];

                c.x[i][j][k] = val;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_vapour_cloud_ice\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: init_vapour_cloud_ice of " << gas << endl;
}
/*
*
*/
void cUranusModel::Forces(){
    cout << endl << "      ATURAN: Forces" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);

    #pragma omp parallel for schedule(static) collapse(2)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            const double rm       = rad.z[i];
            const double sinthe   = sin(the.z[j]);
            const double costhe   = cos(the.z[j]);
            const double rmsinthe = rm * sinthe;

            for(int k = 1; k < km-1; k++){

                const double Coriolis_rad = -2.0 * omega * sinthe * w.x[i][j][k];
                const double Coriolis_the = +2.0 * omega * costhe * w.x[i][j][k];
                const double Coriolis_phi = +2.0 * omega * (-costhe * v.x[i][j][k]
                    + sinthe * u.x[i][j][k]);

                const double dpdr   = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr;
                const double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
                const double dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;

                // u,v,w are non-dimensional (divided by u_0), so restore
                // physical units: F_Cor [N/m³] = ρ · u_0 · |2Ω×v_ND|
                CoriolisForce.x[i][j][k] = Coriolis * r_mix * u_0
                    * sqrt((Coriolis_rad*Coriolis_rad
                          + Coriolis_the*Coriolis_the
                          + Coriolis_phi*Coriolis_phi) / 3.0);

                // rm is dimensionless; physical radius = rm * L_atm [m]
                // F_cen [N/m³] = ρ · ω² · r_phys · (1 + |sinθ|)
                CentrifugalForce.x[i][j][k] = centrifugal * r_mix
                    * omega * omega * rm * L_atm * (1.0 + fabs(sinthe));

                BuoyancyForce.x[i][j][k] = buoyancy
                    * r_mix * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
                    / (r_mix * R_mix * t.x[i][j][k] * t_ref) * 1e5;

                const double dpdthe_rm    = dpdthe / rm;
                const double dpdphi_rmsin = dpdphi / rmsinthe;
                PresGradForce.x[i][j][k] =
                    -sqrt((dpdr*dpdr
                         + dpdthe_rm*dpdthe_rm
                         + dpdphi_rmsin*dpdphi_rmsin) / 3.0) / L_atm * 1.0e5;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: Forces ended" << endl;
}
/*
*
*/

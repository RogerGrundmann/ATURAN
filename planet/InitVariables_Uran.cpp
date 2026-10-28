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
            // ATURAN_T_INIT_ANCHOR=1: the adiabat through t_1bar_K() at 1 bar, the same at
            // every latitude, isothermal at T_top above the level where it has cooled to T_top.
            const bool anchor = t_init_anchor() != 0;
            const double T_bottom = anchor ? t_bottom_anchor_K() : T_top + gam * L_atm;  // K at i=0

            for(int i = 0; i < im; i++){
                height = get_layer_height(i);                           // km
                t.x[i][j][k] = (T_bottom - gam * height) / t_ref;
                if(anchor) t.x[i][j][k] = std::max(T_bottom - gam * height, T_top) / t_ref;

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
/*
 * ATURAN_TW_BALANCE=<strength>: THE ZONAL WIND IN THERMAL-WIND BALANCE WITH THE TEMPERATURE.
 * Default 1 with ATURAN_GEO_BALANCE (the default), else 0. Ported from ATOM_Precipitation (VelocityInitializer.h, balance_thermal_wind),
 * whose note gives the reason it is an initial condition: geostrophic adjustment takes ~1/f,
 * tens of thousands of iterations, and cannot be waited for.
 *
 * theta is the colatitude and w the zonal wind (east positive), so
 *
 *     dw/dz = + g / (a f T) * dT/dtheta,     f = 2 omega cos(theta)
 *
 * The wind at the ANCHOR level is kept (ATURAN_TW_ANCHOR=<level>; default each column's jet
 * level, where the prescribed |w| is largest; ATOM anchors at the surface) and the shear the model's own
 * temperature implies is integrated away from it, trapezoid in height, at the given strength.
 * Tapered out towards the equator with sin2(lat) / (sin2(lat) + sin2(lat_min)), lat_min =
 * ATURAN_TW_LATMIN (15 deg), because f -> 0 there; capped at ATURAN_TW_WMAX (400 m/s). Only w
 * is written (and wn); v and u are left alone.
 */
void cUranusModel::balanceThermalWind(){
    static const double s = [](){ const char* e = getenv("ATURAN_TW_BALANCE"); return e ? atof(e) : (geo_balance() ? 1.0 : 0.0); }();
    if(s == 0.0) return;
    static const double lat_min_deg = [](){ const char* e = getenv("ATURAN_TW_LATMIN"); return e ? atof(e) : 15.0; }();
    static const double w_max = [](){ const char* e = getenv("ATURAN_TW_WMAX"); return e ? atof(e) : 400.0; }();
    static const int anchor_env = [](){ const char* e = getenv("ATURAN_TW_ANCHOR"); return e ? atoi(e) : -1; }();
    // the anchor: a fixed level when the knob names one, otherwise each column's JET LEVEL, the
    // level where the prescribed |w| is largest (the initialiser puts the measured wind there)
    const bool jet_anchor = !(anchor_env >= 0 && anchor_env < im);
    tw_anchor_i.assign(jm * km, im - 1);
    const double s2_min = std::pow(std::sin(lat_min_deg * M_PI / 180.0), 2);
    const double dz = L_atm * 1.0e3 * dr;                           // layer thickness [m]
    const double w_lim = w_max / u_0;

    long n_cols = 0, n_capped = 0;
    double worst_dw = 0.0;
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:n_cols,n_capped) reduction(max:worst_dw)
    for(int j = 1; j < jm-1; j++){
        for(int k = 0; k < km; k++){
            const double slat = std::cos(the.z[j]);                 // sin(latitude)
            const double f = 2.0 * omega * slat;
            if(f == 0.0) continue;
            const double wgt = (slat * slat) / (slat * slat + s2_min);
            // shear at a level [(m/s)/m]
            auto shear = [&](int i)->double{
                const double T = t.x[i][j][k] * t_ref;
                if(!(T > 0.0)) return 0.0;
                const double dTdthe = (t.x[i][j+1][k] - t.x[i][j-1][k]) * t_ref / (2.0 * dthe);
                const double a_m = metricRadius(rad.z[i]) * L_atm * 1.0e3;
                const double sh = g * dTdthe / (a_m * f * T);
                return std::isfinite(sh) ? sh : 0.0; };
            int ia = anchor_env;
            if(jet_anchor){
                ia = 0;
                for(int i = 1; i < im; i++) if(std::fabs(w.x[i][j][k]) > std::fabs(w.x[ia][j][k])) ia = i;
            }
            tw_anchor_i[j * km + k] = ia;
            const double w0 = w.x[ia][j][k];
            auto put = [&](int i, double acc){
                double dw = s * wgt * acc;                          // [m/s]
                if(!std::isfinite(dw)) dw = 0.0;
                double wnew = w0 + dw / u_0;
                if(std::fabs(wnew) > w_lim){ wnew = (wnew > 0.0 ? w_lim : -w_lim); n_capped++; }
                if(std::fabs(dw) > worst_dw) worst_dw = std::fabs(dw);
                w.x[i][j][k] = wnew; wn.x[i][j][k] = wnew; };
            double acc = 0.0;
            for(int i = ia + 1; i < im; i++){ acc += 0.5 * (shear(i-1) + shear(i)) * dz; put(i, acc); }
            acc = 0.0;
            for(int i = ia - 1; i >= 0; i--){ acc -= 0.5 * (shear(i+1) + shear(i)) * dz; put(i, acc); }
            n_cols++;
        }
    }
    double wmx = 0.0, wmn = 0.0;
    for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
        wmx = std::max(wmx, w.x[i][j][k]); wmn = std::min(wmn, w.x[i][j][k]); }
    printf("      ATURAN: thermal-wind balance, strength %.3f, anchor %s (level %d at the equator), lat_min %.1f deg: %ld columns, largest shear added %.2f m/s, w now %.2f .. %.2f m/s, capped cells %ld\n",
           s, jet_anchor ? "each column's jet level" : "fixed level", tw_anchor_i[((jm-1)/2) * km], lat_min_deg, n_cols, worst_dw, wmn * u_0, wmx * u_0, n_capped);
}

void cUranusModel::init_PressureDynamic(){
    cout << endl << "      ATURAN: init_PressureDynamic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ATURAN_PDYN_INIT (2026-10-12): what the dynamic pressure starts from. 1 is the kinetic energy
    // density of the prescribed wind, which this routine has always written; 0 starts it at zero
    // and leaves it to the pressure solver -- ATSAT's repair of the same day (ATSAT_PDYN_INIT,
    // InitVariables_Sat.cpp has the derivation): p_dyn enters the momentum equations through its
    // gradient, and the radial gradient of that kinetic-energy field is a force nothing balances.
    // ATJUP has always started p_dyn at zero.
    // DEFAULT 0 SINCE 2026-10-12 on the user's word (1 restores, byte-identical to 10e724f).
    // On Uranus the field is built with the local density and is small (max 0.027 bar), and the
    // force is correspondingly weaker than on Saturn or Neptune, but it is the same force.
    // MEASURED, 224 iterations, 4 threads, 1 against 0 (ATJUP/satchk/giants/UP1, UP0):
    //     continuity residuum at iteration 16 / 112 / 224   0.032 / 0.066 / 0.051 -> 0.009 / 0.003 / 0.003
    //     radial wind at 224                                -1.68 .. +0.66 -> -0.02 .. +0.18 m/s
    //     H2O / NH3 / H2S column                            -0.38 / -0.10 / -0.02 % -> +0.09 / +0.02 / +0.02 %
    // v, w and the temperature extrema agree within 0.5 %, the mean OLR is 0.92 W/m2 in both.
    static const int pdyn_init = [](){ const char* e = getenv("ATURAN_PDYN_INIT"); return e ? atoi(e) : 2; }();
    // ATURAN_PDYN_GEO=1: THE BAROTROPIC PRESSURE OF A GEOSTROPHIC START. With the radial forces
    // carried by p_hydro, what is left for p_dyn at the start is the meridional balance: the
    // Coriolis force and the metric term w2 cot(theta)/r on the zonal wind need a meridional
    // pressure gradient, the same at every height because the wind below the jet level is. As
    // rhs_v forms them,
    //     scale_p / rm * dp/dtheta = Coriolis * scale_Cor * 2 omega cos(theta) w + w2 cot(theta) / rm
    // integrated from the north pole with the anchor wind of each column (balanceThermalWind),
    // trapezoid in theta. The centrifugal term of rhs_v is 1e-4 of these and is left out.
    if(pdyn_geo() && metric_split()){
        const double scale_p = 1.0e5 / (r_mix * u_0 * u_0);
        // ATURAN_PDYN_GEO_EXACT (2026-10-14, ported from ATNEPT; default follows ATURAN_GEO_BALANCE,
        // 0 = the trapezoid): integrate so that the centred difference rhs_v takes returns the
        // integrand exactly at every integrated row, p[j+1] = p[j-1] + 2 dthe G[j], seeded with the
        // trapezoid's first step. The trapezoid's centred difference is the 1-2-1 mean of G, wrong
        // where the jets change over a few rows. Two interleaved chains: a 2-cell mode in p is
        // invisible to the centred difference that reads it.
        static const bool geo_exact = [](){ const char* e = getenv("ATURAN_PDYN_GEO_EXACT"); return e ? atoi(e) != 0 : geo_balance(); }();
        auto integrate = [&](const std::vector<double>& G, std::vector<double>& P){
            P[0] = 0.0;
            for(int j = 1; j < jm; j++)
                P[j] = (geo_exact && j >= 2) ? P[j-2] + 2.0 * dthe * G[j-1]
                                             : P[j-1] + 0.5 * (G[j-1] + G[j]) * dthe;
        };
        double pmax = 0.0, pmin = 0.0;
        #pragma omp parallel for schedule(static) reduction(max:pmax) reduction(min:pmin)
        for(int k = 0; k < km; k++){
            std::vector<double> G(jm, 0.0), P(jm, 0.0);
            for(int j = 0; j < jm; j++){
                const int ia = ((int)tw_anchor_i.size() == jm * km) ? tw_anchor_i[j * km + k] : im - 1;
                const double rm = metricRadius(rad.z[ia]);
                const double sth = std::max(sinthe_min(), std::abs(sin(the.z[j])));
                const double cth = cos(the.z[j]);
                const double wa = w.x[ia][j][k];
                G[j] = (j == 0 || j == jm-1 || sth < 1.0e-6) ? 0.0
                    : (Coriolis * coriolis_scale() * 2.0 * omega * cth * wa * rm + wa * wa * cth / sth) / scale_p;
            }
            integrate(G, P);
            for(int j = 0; j < jm; j++){
                for(int i = 0; i < im; i++) p_dyn.x[i][j][k] = P[j];
                pmax = std::max(pmax, P[j]); pmin = std::min(pmin, P[j]);
            }
        }
        printf("      ATURAN: init_PressureDynamic: p_dyn is the barotropic pressure that balances Coriolis and metric force of the zonal wind in the meridional equation; %.6f .. %.6f bar (ATURAN_PDYN_GEO=0 starts it at zero)\n", pmin, pmax);
        // ATURAN_PDYN_GEO=2: the rest of the balance, level by level (see pdyn_geo_mode()). rhs_v
        // takes -scale_p dp_dyn/dtheta / rm - dp_hydro/dtheta / rm, so p_baro = Q_i - scale_p p_dyn
        // - p_hydro(start), Q_i the same integral of the same integrand with the wind and metric
        // radius of level i. p_hydro(start) includes the buoyancy of the hydrostatic split.
        if(pdyn_geo_mode() == 2){
            p_baro.assign((size_t)im * jm * km, 0.0);
            computeBuoyancyRefLevel();
            computeHydrostaticPressure();                       // p_baro still zero
            double bmax = 0.0, bmin = 0.0;
            #pragma omp parallel for schedule(static) reduction(max:bmax) reduction(min:bmin)
            for(int k = 0; k < km; k++){
                std::vector<double> G(jm, 0.0), Q(jm, 0.0);
                for(int i = 0; i < im; i++){
                    const double rm = metricRadius(rad.z[i]);
                    for(int j = 0; j < jm; j++){
                        const double sth = std::max(sinthe_min(), std::abs(sin(the.z[j])));
                        const double cth = cos(the.z[j]);
                        const double wi = w.x[i][j][k];
                        G[j] = (j == 0 || j == jm-1 || sth < 1.0e-6) ? 0.0
                            : Coriolis * coriolis_scale() * 2.0 * omega * cth * wi * rm + wi * wi * cth / sth;
                    }
                    integrate(G, Q);
                    for(int j = 0; j < jm; j++){
                        const double b_ = Q[j] - scale_p * p_dyn.x[i][j][k] - p_hydro.x[i][j][k];
                        p_baro[((size_t)i * jm + j) * km + k] = b_;
                        bmax = std::max(bmax, b_); bmin = std::min(bmin, b_);
                    }
                }
            }
            computeHydrostaticPressure();                       // now with p_baro
            printf("      ATURAN: init_PressureDynamic: ATURAN_PDYN_GEO=2, p_baro (the per-level rest, added to p_hydro, frozen) %.6e .. %.6e\n", bmin, bmax);
        }
        return;
    }
    // With ATURAN_METRIC_SPLIT the metric term is balanced by p_hydro, so a p_dyn that balances it
    // too would be an unopposed force: it starts at zero then.
    if(metric_split() && pdyn_init == 2){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++) p_dyn.x[i][j][k] = 0.0;
        cout << "      ATURAN: init_PressureDynamic: p_dyn starts at zero, the metric term is carried by p_hydro (ATURAN_METRIC_SPLIT=0 restores the balanced p_dyn)" << endl;
        return;
    }
    // 2 (2026-10-12): p_dyn starts in balance with the metric term of the radial momentum
    // equation. rhs_u carries + (v2 + w2) / rm, the centrifugal acceleration of the horizontal
    // flow itself, and the prescribed jet makes it large (w = 400 m/s gives 6.5e-3 m/s2 on
    // Neptune). With p_dyn = 0 nothing opposes it: it pushes u upward in every column, and
    // against a rigid lid that is a divergence in the wall cells. The pressure solver would
    // build the balancing gradient itself (the term is in its source), but one sweep per call
    // does not get there. So it is integrated here, d p_dyn / dr = (v2 + w2) / (rm * scale_p),
    // from zero at the bottom, trapezoid in r, with the RHS's own scale_p = 1e5 / (r_mix u_0^2).
    // DEFAULT SINCE 2026-10-12 on the user's word ("repair ... the raised continuity residuum from
    // the rigid lid"). MEASURED with a per-term print of rhs_u (not kept): at iteration 16 the term
    // is the whole of the level-mean tendency of u, everything else 20 times smaller or more; with
    // 2 the pressure gradient cancels it in the interior, a remainder stays in the first and last
    // cell. 224 iterations, rigid lid on, 0 against 2 (ATJUP/satchk/giants/UL224, UB224):
    //     continuity residuum at 16 / 224   0.0127 / 0.0133 -> 0.0095 / 0.0043
    //     radial wind at 224                -0.01 .. +0.19 -> -0.04 .. +0.03 m/s
    //     H2O / NH3 / H2S column            +0.069 / +0.020 / +0.019 % -> +0.052 / -0.003 / -0.004 %
    // v, w, the temperature extrema, the OLR and the precipitation are unchanged.
    if(pdyn_init == 2){
        const double scale_p = 1.0e5 / (r_mix * u_0 * u_0);
        double pmax = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(max:pmax)
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                double f_prev = 0.0, p_ = 0.0;
                for(int i = 0; i < im; i++){
                    const double f = (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k])
                                   / (metricRadius(rad.z[i]) * scale_p);
                    if(i > 0) p_ += 0.5 * (f_prev + f) * dr;
                    p_dyn.x[i][j][k] = p_; f_prev = f;
                    if(p_ > pmax) pmax = p_;
                }
            }
        }
        printf("      ATURAN: init_PressureDynamic: p_dyn balances the metric term (v2 + w2)/rm of rhs_u; max %.6f bar (ATURAN_PDYN_INIT=0 starts it at zero)\n", pmax);
        return;
    }
    if(pdyn_init == 0){
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++) p_dyn.x[i][j][k] = 0.0;
        cout << "      ATURAN: init_PressureDynamic: p_dyn starts at zero (ATURAN_PDYN_INIT=1 restores the kinetic-energy field)" << endl;
        return;
    }

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
    printf(" time measured: %.3f seconds for init_PressureDynamic\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: init_PressureDynamic ended" << endl;
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
            // ATURAN_T_INIT_ANCHOR=1: above the level where the adiabat has cooled to the top
            // temperature the layer is isothermal, and its pressure is the hydrostatic one of an
            // isothermal layer (same gas constant as the exponent). r_mix has been scaled to the
            // anchored bottom state by then, so p_bottom IS p_bottom_bar().
            const bool anchor = t_init_anchor() != 0;
            const double T_bottom = t.x[0][j][k] * t_ref;
            const double T_top    = t.x[im-1][j][k] * t_ref;
            const double h_trop   = (T_bottom - T_top) / gam;           // km
            const double p_trop   = p_bottom * pow(T_top / T_bottom, exp_pressure);
            for(int i = 0; i < im; i++){
                p_stat.x[i][j][k] = p_bottom
                    * pow(t.x[i][j][k] / t.x[0][j][k], exp_pressure);
                if(anchor && get_layer_height(i) > h_trop)
                    p_stat.x[i][j][k] = p_trop * exp(-(get_layer_height(i) - h_trop) * g / (R_ref * T_top));
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
    // ATURAN_INIT_SUPERSAT=<f> (2026-10-28, ATJUP row GIANTS-MCA-DEFAULT; DEFAULT 1 since 2026-10-28, user: "(a)";
    // =1000 restores the 1e3 above): the vapour set-up
    // is magnus x saturation, capped at the deep abundance, so with 1e3 every deck starts up to 1000x
    // supersaturated; the set-up saturation adjustment condenses the excess and its latent heat
    // warms the decks +14.5 K (i1) / +18.7 K (i12) over the dry adiabat at the equator -- the steep
    // pairs the MCA fires on (satchk/mcagiants_2026-10-28/dbg). 1 = vapour at saturation.
    static const double init_supersat = [](){ const char* e = getenv("ATURAN_INIT_SUPERSAT"); return e ? atof(e) : 1.0; }();
    if(init_supersat > 0.0) magnus = init_supersat;

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

                const double fr = init_species_fr(i, j, k);   // <TAG>_INIT_SPECIES_LOCAL_RHO, 1.0 unset
                double val = magnus * fr * r_mix * q_Rain;

                if(val >= fr * r_max[j])            val = fr * r_max[j];
                if((gas == "NH3") && (p_u >= p_00_nh3)) val = fr * r_max[j];
                if((gas == "H2S") && (p_u >= p_00_h2s)) val = fr * r_max[j];
                if((gas == "H2O") && (p_u >= p_0_h2o))  val = fr * r_max[j];
                if((gas == "CH4") && (p_u >= p_0_ch4))  val = fr * r_max[j];

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

    // ATURAN_BC_PHI_PERIODIC (ATPhys.h): with 1 the seam meridian k = 0 is computed like any other.
    const bool phi_per = ATPhys::phi_periodic<cUranusModel>();
    const int  k_first = ATPhys::phi_first(phi_per);

    #pragma omp parallel for schedule(static) collapse(2)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            const double rm       = metricRadius(rad.z[i]);
            const double sinthe   = sin(the.z[j]);
            const double costhe   = cos(the.z[j]);
            const double rmsinthe = rm * sinthe;

            for(int k = k_first; k < km-1; k++){
                const int kw = ATPhys::phi_wrap(phi_per, k-1, km);
                const int ke = ATPhys::phi_wrap(phi_per, k+1, km);

                const double Coriolis_rad = -2.0 * omega * sinthe * w.x[i][j][k];
                const double Coriolis_the = +2.0 * omega * costhe * w.x[i][j][k];
                const double Coriolis_phi = +2.0 * omega * (-costhe * v.x[i][j][k]
                    + sinthe * u.x[i][j][k]);

                const double dpdr   = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr;
                const double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
                const double dpdphi = (p_dyn.x[i][j][ke] - p_dyn.x[i][j][kw]) * inv_2dphi;

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

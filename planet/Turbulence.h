/*
 * SHARED PHYSICS — turbulence closure (k-epsilon / k-omega / k-omega SST), one implementation
 * for every planet. MUST BE BYTE-IDENTICAL IN EVERY MODEL THAT USES IT; `make check-shared`
 * verifies that against planet/SHARED.md5.
 *
 * It knows nothing about which planet it runs on. Everything planet-specific arrives through the
 * Planet template parameter, and it asks only three things of a model beyond the usual fields:
 *
 *     Planet::planet_tag()        "ATJUP" / "ATSAT" — the log prefix, and the prefix its knobs
 *                                are read under (ATJUP_TURB_MODEL, ATSAT_NUE_MAX, ...)
 *     m.surface_index(j,k)       the first fluid level of a column
 *     m.is_solid(i,j,k)          whether a cell is inside the topography
 *
 * THE LAST TWO ARE THE ENTIRE DIFFERENCE between the two copies this replaces. TurbulenceJup.h
 * and TurbulenceSat.h were 837 and 841 lines differing in 42, of which 28 were comments — every
 * one of the six real differences was the surface of a column, written as i_topography[j][k] and
 * SeaMount.x == 1.0 in one and as 0 and false in the other. Behind two accessors, the same code
 * serves a model with an obstacle and a model without.
 */
/*
 * Jupiter Atmosphere Circulation Model (ATJUP) applied to turbulent flow
 * Standalone turbulence class: k-epsilon (Chien 1982),
 *                              k-omega   (Wilcox 1988),
 *                              k-omega SST (Menter 1994)
 *
 * A faithful mirror of ATOM_Precipitation atmosphere/TurbulenceAtm.h. All three
 * models are kept, selected by cJupiterModel::turb_model ("k_epsilon",
 * "k_omega", "k_omega_SST"), overridable at runtime with ATJUP_TURB_MODEL.
 * Turbulence arrays (tke, tken, dis, disn, nue, prod, tke_source, dis_source and
 * the 2D vel_star) are owned by cJupiterModel and reached through the reference m.
 *
 * The equations, coefficients, blending functions, limiters, caps and Neumann
 * land-face treatment are carried over unchanged. Only the four things that MUST
 * differ between the two codebases were adapted:
 *
 *  1. UNITS. ATOM's L_atm is in METRES (param.py: 400.0); ATJUP's is in
 *     KILOMETRES (140.0), and ATJUP's get_layer_height() likewise returns km
 *     where ATOM's returns m. Every place ATOM writes m.L_atm or
 *     get_layer_height() this file uses L_atm_m = m.L_atm*1e3 and
 *     height_m(i) = get_layer_height(i)*1e3. Copying the expressions verbatim
 *     would put every nondimensionalisation out by 1000x. The dimensionless
 *     convention itself is unchanged: k* = k/u_0^2, omega* = omega*L_atm/u_0,
 *     eps* = eps*L_atm/u_0^3, nue* = nue/(u_0*L_atm).
 *
 *  2. LAND MASK. ATOM carries an h array with AtomUtils::is_land(). ATJUP marks
 *     solid cells with SeaMount.x == 1.0 (BC_Jup.h bcSeaMount), so is_land() here
 *     tests that instead. Note what this means physically: on Earth the mask is
 *     the planetary surface, whereas ATJUP's only solid body is the SeaMount, the
 *     obstacle installed to make a GRS-like storm possible. So the wall-bounded
 *     branches of these models act on that obstacle, and i=0 is the deep interior
 *     of a gas giant rather than a ground surface. The i=0 wall BC is retained as
 *     ATOM has it, since this is a mirror; whether a deep-interior wall condition
 *     is the physics wanted is a separate question from whether the port is faithful.
 *
 *  3. C++ STANDARD. ATOM uses std::clamp (C++17); ATJUP builds with -std=c++11,
 *     so those become std::max(lo, std::min(hi, x)). Also, C++11 has no inline
 *     variables, so any constexpr member passed to std::min/std::max by const&
 *     would be odr-used and fail to link from a header — hence the local copies.
 *
 *  4. GAS PROPERTIES. nue_air defaults to the H2/He kinematic viscosity rather
 *     than air's 1.5e-5 m2/s.
 *
 * Gated by ATJUP_TURB (default 0 = off, bit-identical). Like ATOM, this fills
 * tke/dis/nue/prod/tke_source/dis_source. nue* feeds the momentum and scalar equations
 * (ATJUP_TURB_COUPLING), and — as in ATOM — the two turbulence transport equations are
 * assembled in RHS_Jup_Turb.cpp and integrated by RungeKutta_Jup_Turb.cpp, which re-derive
 * the P-Y balance from the current k* and dis* at every RK4 sub-stage. What stays here is the
 * eddy viscosity itself, the friction velocity, the wall/ABL conditioning and the masking
 * of the solid cells.
*/

#pragma once

#include "ConvectiveAdjustment.h"   // ATPhys::env_for and friends

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

class Array;

template<class Planet>
class Turbulence {
public:
    enum Model { k_epsilon, k_omega, k_omega_SST };

    explicit Turbulence(Planet& model,
                           double vel_star_ref = 0.4,        // reference friction velocity for re_turb [m/s]
                           double z_0          = 0.1,        // roughness length   [m]
                           double nue_gas      = 1.8e-5);    // kin. viscosity H2/He [m2/s]

    void init();
    void run();
    void apply_wall_bc();

    // Eddy-viscosity ceiling in PHYSICAL units [m2/s]. ATOM caps at 1000 m2/s, a terrestrial
    // boundary-layer value — but this closure's own equilibrium here is nu_t = k/omega ~ 5e3
    // m2/s (k ~ 67 m2/s2, omega ~ 0.014 1/s), so 1000 sat BELOW the physics and the cap, not the
    // turbulence model, was setting the eddy viscosity everywhere. Jovian estimates: vertical
    // eddy diffusivity ~1e3-1e4 m2/s in the troposphere, jet-scale horizontal eddy viscosity
    // ~1e5-1e6. The default therefore sits above the closure's own value so it acts as a
    // runaway guard rather than the operative limiter. Override with ATJUP_NUE_MAX.
    static double nue_max_phys(){
        static const double v = [](){ const char* e = ATPhys::env_for(Planet::planet_tag(), "NUE_MAX"); return e ? atof(e) : 1.0e5; }();
        return v;
    }
    // ATOM tapers nu_t to zero above abl_height because residual viscosity aloft homogenised the
    // shear its Hadley-Ferrel cells live on. Jupiter has no such boundary-layer ceiling —
    // turbulence here is generated by jet shear and by the SeaMount wake throughout the free
    // atmosphere — so the taper is OFF by default. ATJUP_ABL_TAPER=1 restores ATOM's behaviour.
    static int abl_taper(){
        static const int v = [](){ const char* e = ATPhys::env_for(Planet::planet_tag(), "ABL_TAPER"); return e ? atoi(e) : 0; }();
        return v;
    }

private:
    Planet& m;
    Model  turb_model;
    double vel_star_ref;
    double z_0;
    double nue_air;      // kept the ATOM name so the formulae read identically
    double re_turb;

    // ABL model constants (unchanged from ATOM)
    static constexpr double C_nue   = 0.028;
    static constexpr double Karman  = 0.42;
    static constexpr double zeta    = 3.715;      // coordinate-stretching factor
    static constexpr double dis_min = 1.0e-10;    // minimum dis* to avoid nue -> inf

    static Model parse_model(const std::string& s){
        if(s == "k_epsilon")   return k_epsilon;
        if(s == "k_omega")     return k_omega;
        if(s == "k_omega_SST") return k_omega_SST;
        return k_omega_SST;                        // default / "none" fallback
    }
    static double blend(double inner, double outer, double F1){
        return F1 * inner + (1.0 - F1) * outer;
    }
    // C++11 replacement for std::clamp
    static double clamp3(double x, double lo, double hi){
        return (x < lo) ? lo : ((x > hi) ? hi : x);
    }

    void print_model_name() const;
    bool is_land(int i, int j, int k) const;
    double height_m(int i) const;      // get_layer_height in METRES (the models store km)
    double L_atm_m() const;            // L_atm in METRES

    void compute_vel_star();
    void init_fields();
    void compute_sources();
    void zero_land_cells();
    void clamp_nue();

    void compute_k_epsilon(int i, int j, int k, double y_mount,
        double dudr, double dvdr, double dwdr,
        double dudthe, double dvdthe, double dwdthe,
        double dudphi, double dvdphi, double dwdphi,
        double Omega);
    void compute_k_omega(int i, int j, int k,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double Omega_mag, double W_unused);
    void compute_k_omega_SST(int i, int j, int k,
        double y_mount, double rm, double sinthe,
        double dtkedr, double ddisdr,
        double dtkedthe, double ddisdthe,
        double dtkedphi, double ddisdphi,
        double Omega, double W);
};

// ---------------------------------------------------------------------------

template<class Planet>
inline Turbulence<Planet>::Turbulence(Planet& model,
                                    double vel_star_ref_, double z_0_, double nue_gas_)
    : m(model),
      turb_model(k_omega_SST),
      vel_star_ref(vel_star_ref_),
      z_0(z_0_),
      nue_air(nue_gas_),
      re_turb(vel_star_ref_ * z_0_ / nue_gas_)
{
    // Model selection: the model member, overridable from the environment.
    const char* e = ATPhys::env_for(Planet::planet_tag(), "TURB_MODEL");
    turb_model = parse_model(e ? std::string(e) : m.turb_model);
}

template<class Planet>
inline bool Turbulence<Planet>::is_land(int i, int j, int k) const {
    return m.is_solid(i, j, k);
}
template<class Planet>
inline double Turbulence<Planet>::height_m(int i) const {
    return (double)m.get_layer_height(i) * 1.0e3;    // km -> m
}
template<class Planet>
inline double Turbulence<Planet>::L_atm_m() const {
    return m.L_atm * 1.0e3;                          // km -> m
}

template<class Planet>
inline void Turbulence<Planet>::print_model_name() const {
    using namespace std;
    switch(turb_model){
        case k_epsilon:   cout << "      k-epsilon turbulence model (Chien 1982)" << endl; break;
        case k_omega:     cout << "      k-omega turbulence model (Wilcox 1988)" << endl;  break;
        case k_omega_SST: cout << "      k-omega SST turbulence model (Menter 1994)" << endl; break;
    }
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::init(){
    using namespace std;
    cout << endl << "      " << Planet::planet_tag() << ": Turbulence<Planet>::init" << endl;
    print_model_name();
    auto begin = chrono::high_resolution_clock::now();

    m.re_turb = re_turb;
    compute_vel_star();
    init_fields();
    apply_wall_bc();
    compute_sources();   // prime nue, tke_source, dis_source from the ABL profile
    zero_land_cells();
    clamp_nue();

    auto end     = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Turbulence<Planet>::init\n", elapsed.count() * 1e-9);
    cout << "      " << Planet::planet_tag() << ": Turbulence<Planet>::init ended" << endl;
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::run(){
    using namespace std;
    cout << endl << "      " << Planet::planet_tag() << ": Turbulence<Planet>::run" << endl;
    print_model_name();
    auto begin = chrono::high_resolution_clock::now();

    compute_vel_star();
    compute_sources();
    // No time integration here: k* and dis* are prognostic variables of the RK4 system
    // (rhs_tke / rhs_dis in RHS_Jup_Turb.cpp), exactly as in ATOM. Integrating them here as
    // well would advance the same source terms twice per iteration.
    apply_wall_bc();
    zero_land_cells();
    clamp_nue();

    auto end     = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Turbulence<Planet>::run\n", elapsed.count() * 1e-9);
    cout << "      " << Planet::planet_tag() << ": Turbulence<Planet>::run ended" << endl;
}

// -----------------------------------------------------------------------
// Friction velocity u_tau per (j,k) column from the horizontal wind speed at the
// first fluid cell above the topographic surface.
//   U_horiz = u_0*sqrt(v^2 + w^2)  at i = i_mount + 1
//   u_tau   = max(vel_star_min, sqrt(C_D)*U_horiz)
template<class Planet>
inline void Turbulence<Planet>::compute_vel_star(){
    const double C_D          = 0.002;
    const double vel_star_min = 0.05;   // [m/s]
    // Ceiling on u_tau, from the strongest horizontal wind Jupiter is known to sustain
    // (~200 m/s in the equatorial jet, twice the model's velocity scale u_0):
    //   u_tau_max = sqrt(C_D) * U_max.
    // The reason this is needed at all: ATJUP's momentum field develops non-finite and
    // absurd (1e150 m/s) cells near the pole and at the first fluid layer from iteration 1,
    // with or without the closure — a pre-existing problem of the momentum solver. u_tau
    // feeds k_bg, which compute_sources applies as a FLOOR on k*, so one such cell used to
    // seed k* ~ 1e304 across the domain the moment k* became prognostic. Bounding u_tau
    // keeps the closure's response to a locally failed momentum field finite and local; it
    // is the counterpart of ATOM's velocity clamp in BC_Atm.h, done here so that ATJUP's
    // momentum fields themselves are left untouched. Override with ATJUP_VEL_STAR_MAX [m/s].
    static const double vel_star_max = [](){
        const char* e = ATPhys::env_for(Planet::planet_tag(), "VEL_STAR_MAX");
        return e ? atof(e) : std::sqrt(0.002) * 200.0; }();   // ~8.9 m/s

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){
            const int i_surf = std::min(m.surface_index(j, k) + 1, m.im - 1);
            const double v_s = m.v.x[i_surf][j][k];
            const double w_s = m.w.x[i_surf][j][k];
            // A non-finite velocity in the sampled cell must not become a non-finite u_tau:
            // k_bg in init_fields/compute_sources is proportional to u_tau^2 and is applied as
            // a FLOOR on k*, so one infinite column would seed k* = inf and — now that k* is
            // integrated by the RK4 — spread it. ATJUP does produce non-finite v,w at the first
            // fluid layer near the pole from iteration 1 (also with the closure off: it is a
            // pre-existing problem of the momentum solver, not of this model), so the guard is
            // not hypothetical. Such a column falls back to the calm-air floor.
            const double vw2 = v_s * v_s + w_s * w_s;
            if(!std::isfinite(vw2)){
                m.vel_star.y[j][k] = vel_star_min;
                continue;
            }
            const double U_horiz = m.u_0 * std::sqrt(vw2);
            m.vel_star.y[j][k] = clamp3(std::sqrt(C_D) * U_horiz, vel_star_min, vel_star_max);
        }
    }
}

// -----------------------------------------------------------------------
// Wall boundary condition for omega at i=0. Physical formulas give dimensional
// omega [1/s]; normalise to omega* = omega*L_atm/u_0.
template<class Planet>
inline void Turbulence<Planet>::apply_wall_bc(){
    // ATOM applies this at i=0, its ground surface. ATJUP's i=0 is the DEEP INTERIOR of a gas
    // giant: there is no wall there, and imposing a no-slip log-layer condition invented a
    // surface boundary layer that does not exist. The only solid body in ATJUP is the SeaMount,
    // so the wall law acts at the first fluid cell above it and nowhere else. Columns without
    // an obstacle get a zero-gradient (non-amplifying) condition at i=0 instead.
    const double nd_omega = L_atm_m() / m.u_0;                 // omega_phys -> omega*
    const double C_nue_l  = C_nue;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){
            const int i_w = m.surface_index(j, k);

            if(i_w > 0 && i_w < m.im - 1){
                // Wall face lies between the solid cell i_w-1 and the first fluid cell i_w.
                const double y1      = std::max(m.layer_thickness_m(std::max(i_w - 1, 0)), 1.0);
                const double vs      = m.vel_star.y[j][k];

                // tke: log-layer equilibrium k = u_tau^2/sqrt(C_mu) (all models)
                m.tke.x[i_w][j][k] = vs * vs / (std::sqrt(C_nue_l) * m.u_0 * m.u_0);
                // nue vanishes at the no-slip wall
                m.nue.x[i_w][j][k] = 0.0;

                if(turb_model == k_epsilon){
                    m.dis.x[i_w][j][k] = 0.0;                  // Chien: eps_wall = 0
                } else if(turb_model == k_omega){
                    m.dis.x[i_w][j][k] = 6.0 * nue_air / (0.075 * y1 * y1) * nd_omega;
                } else if(turb_model == k_omega_SST){
                    m.dis.x[i_w][j][k] = 60.0 * nue_air / (0.0333 * y1 * y1) * nd_omega;
                }
            } else {
                // No obstacle in this column: zero-gradient at the deep boundary.
                m.tke.x[0][j][k] = m.tke.x[1][j][k];
                m.dis.x[0][j][k] = m.dis.x[1][j][k];
                m.nue.x[0][j][k] = m.nue.x[1][j][k];
            }
        }
    }
}

// -----------------------------------------------------------------------
// Initialise tke and dis throughout the domain (ABL profile, Narjisse).
template<class Planet>
inline void Turbulence<Planet>::init_fields(){
    const double inv_u02  = 1.0 / (m.u_0 * m.u_0);
    const double L_m      = L_atm_m();
    const double nd_omega = L_m / m.u_0;                          // omega_phys -> omega*
    const double nd_eps   = L_m / (m.u_0 * m.u_0 * m.u_0);        // eps_phys   -> eps*
    const double C_nue_l  = C_nue;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 1; j < m.jm - 1; j++){
        for(int k = 1; k < m.km - 1; k++){
            // Per-column friction velocity and derived background floors.
            const double vs     = m.vel_star.y[j][k];
            const double k_bg   = 0.01 * vs * vs / std::sqrt(C_nue_l);                        // [m2/s2]
            const double om_bg  = std::pow(C_nue_l, -0.25) * std::sqrt(k_bg) / (Karman * L_m); // [1/s]
            const double eps_bg = std::pow(C_nue_l, 0.75) * std::pow(k_bg, 1.5) / (Karman * L_m);

            // ABL anchored to the LOCAL surface, so it rides on top of topography.
            const double y_mount = height_m(m.surface_index(j, k));   // surface height [m]

            for(int i = 1; i < m.im - 1; i++){
                const double z_i    = height_m(i) - y_mount;         // height above local surface [m]
                // Wall distance floored at ONE GRID LAYER, not ATOM's 1e-6 m. z_i is exactly 0
                // in the layer sitting on the local surface, and omega = C^-0.25 sqrt(k)/(kappa z)
                // then diverges. ATOM tolerates the 1e-6 floor only because its nd_omega =
                // L_atm/u_0 is 4 s (L_atm = 400 m); ATJUP's is 1400 s (L_atm = 140 km), 350x
                // larger, so the same singularity reached omega* ~ 1e10 and dis_source ~ -1e20.
                // Flooring at one layer is what ATOM's own compute_k_epsilon does for its wall
                // distance (max(height - y_mount, L_atm)), so this is its practice, not a new idea.
                const double dz_layer = std::max(1.0, m.layer_thickness_m(std::max(i - 1, 0)));
                const double z_safe   = std::max(z_i, dz_layer);
                const double z_frac = std::min(z_i, m.abl_height) / m.abl_height;

                // Physical TKE [m2/s2]: parabolic ABL profile with background floor
                const double k_abl  = vs * vs / std::sqrt(C_nue_l) * std::pow(1.0 - z_frac, 2);
                const double k_phys = std::max(k_abl, k_bg);

                double dis_nd = 0.0;
                if(turb_model == k_epsilon){
                    // eps = C_mu^0.75 k^1.5/(kappa z)
                    const double eps_phys = std::pow(C_nue_l, 0.75)
                        * std::pow(k_phys, 1.5) / (Karman * z_safe);
                    dis_nd = std::max(eps_phys, eps_bg) * nd_eps;
                } else {
                    // omega = C_mu^-0.25 k^0.5/(kappa z), floored at om_bg so nue=k/omega stays finite
                    const double om_phys = std::pow(C_nue_l, -0.25)
                        * std::sqrt(k_phys) / (Karman * z_safe);
                    dis_nd = std::max(om_phys, om_bg) * nd_omega;
                }

                const double dis_min_l = dis_min;
                m.tke.x[i][j][k]  = std::max(0.0, k_phys * inv_u02);
                m.dis.x[i][j][k]  = std::max(dis_min_l, dis_nd);
                m.tken.x[i][j][k] = m.tke.x[i][j][k];
                m.disn.x[i][j][k] = m.dis.x[i][j][k];
            }
        }
    }
}

// -----------------------------------------------------------------------
// Global nue clamp, plus confinement of the eddy viscosity to the boundary layer:
// above abl_height nue* is tapered to zero over ~one grid layer.
template<class Planet>
inline void Turbulence<Planet>::clamp_nue(){
    const double nue_max = nue_max_phys() / (m.u_0 * L_atm_m());

    const bool taper_on = (abl_taper() != 0);
    int i_abl = 0;
    while(i_abl < m.im - 1 && height_m(i_abl) <= m.abl_height) i_abl++;
    const double taper_dz = std::max(1.0,
          height_m(std::min(i_abl, m.im - 1))
        - height_m(std::max(i_abl - 1, 0)));

    #pragma omp parallel for collapse(3) schedule(static)
    for(int i = 0; i < m.im; i++){
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                double f = 1.0;
                if(taper_on){
                    const double z_agl = height_m(i) - height_m(m.surface_index(j, k));
                    f = clamp3((m.abl_height + taper_dz - z_agl) / taper_dz, 0.0, 1.0);
                }
                m.nue.x[i][j][k] = clamp3(m.nue.x[i][j][k], 0.0, nue_max) * f;
            }
        }
    }
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::zero_land_cells(){
    #pragma omp parallel for collapse(3) schedule(static)
    for(int i = 0; i < m.im; i++){
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                if(is_land(i, j, k)){
                    m.tke.x[i][j][k]        = 0.0;
                    m.tken.x[i][j][k]       = 0.0;
                    m.dis.x[i][j][k]        = 0.0;
                    m.disn.x[i][j][k]       = 0.0;
                    m.nue.x[i][j][k]        = 0.0;
                    m.prod.x[i][j][k]       = 0.0;
                    m.tke_source.x[i][j][k] = 0.0;
                    m.dis_source.x[i][j][k] = 0.0;
                }
            }
        }
    }
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::compute_sources(){
    const double C_nue_l   = C_nue;
    const double dis_min_l = dis_min;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 1; j < m.jm - 1; j++){
        for(int k = 1; k < m.km - 1; k++){
            const int    i_mount = m.surface_index(j, k);
            const double y_mount = height_m(i_mount);               // surface height [m]

            const double vs      = m.vel_star.y[j][k];
            const double k_bg_nd = 0.01 * vs * vs / (std::sqrt(C_nue_l) * m.u_0 * m.u_0);

            for(int i = 1; i < m.im - 1; i++){
                if(is_land(i, j, k)) continue;

                // ---- geometry ----
                const double rm        = m.rad.z[i];
                // coord_stretching, like everywhere else. This module applied 1/(rm+1)
                // UNCONDITIONALLY while RungeKutta_Jup_Turb and PressureSolverJup both gate it
                // on the flag, which defaults to false. Every radial velocity gradient feeding
                // the production tensor was therefore divided by (rm+1) — 2.5 in the original
                // geometry, 501 with ATJUP_METRIC_RADIUS — and the production is dominated by
                // radial shear, so it entered squared.
                const double exp_rm    = m.coord_stretching ? 1.0 / (rm + 1.0) : 1.0;
                double sinthe          = std::sin(m.the.z[j]);
                if(sinthe == 0.0) sinthe = 1.0e-5;
                const double costhe    = std::cos(m.the.z[j]);
                const double cotanthe  = costhe / sinthe;
                const double inv_rm    = 1.0 / rm;
                const double rmsinthe  = rm * sinthe;
                const double inv_2dr   = 1.0 / (2.0 * m.dr);
                const double inv_2dthe = 1.0 / (2.0 * m.dthe);
                const double inv_2dphi = 1.0 / (2.0 * m.dphi);
                const double inv_rm2dthe       = inv_2dthe / rm;
                const double inv_rmsinthe2dphi = inv_2dphi / rmsinthe;

                // Enforce floors before any reads.
                m.tke.x[i][j][k] = std::max(0.0,       m.tke.x[i][j][k]);
                m.dis.x[i][j][k] = std::max(dis_min_l, m.dis.x[i][j][k]);

                // ---- velocity gradients ----
                const double dudr = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) * inv_2dr * exp_rm;
                const double dvdr = (m.v.x[i+1][j][k] - m.v.x[i-1][j][k]) * inv_2dr * exp_rm;
                const double dwdr = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) * inv_2dr * exp_rm;

                const double dudthe = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) * inv_rm2dthe;
                const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) * inv_rm2dthe;
                const double dwdthe = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) * inv_rm2dthe;

                const double dudphi = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) * inv_rmsinthe2dphi;
                const double dvdphi = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) * inv_rmsinthe2dphi;
                const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) * inv_rmsinthe2dphi;

                // ---- tke/dis gradients with Neumann BC at land faces ----
                const bool neumann_bot = (i > 0        && is_land(i-1, j, k));
                const bool neumann_jm1 = (j > 0        && is_land(i, j-1, k));
                const bool neumann_jp1 = (j < m.jm - 1 && is_land(i, j+1, k));
                const bool neumann_km1 = (k > 0        && is_land(i, j, k-1));
                const bool neumann_kp1 = (k < m.km - 1 && is_land(i, j, k+1));

                const double tke_im1 = neumann_bot ? m.tke.x[i][j][k] : m.tke.x[i-1][j][k];
                const double dis_im1 = neumann_bot ? m.dis.x[i][j][k] : m.dis.x[i-1][j][k];
                const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
                const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
                const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
                const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
                const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
                const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
                const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
                const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

                const double dtkedr   = (m.tke.x[i+1][j][k] - tke_im1) * inv_2dr * exp_rm;
                const double ddisdr   = (m.dis.x[i+1][j][k] - dis_im1) * inv_2dr * exp_rm;
                const double dtkedthe = (tke_jp1 - tke_jm1) * inv_rm2dthe;
                const double ddisdthe = (dis_jp1 - dis_jm1) * inv_rm2dthe;
                const double dtkedphi = (tke_kp1 - tke_km1) * inv_rmsinthe2dphi;
                const double ddisdphi = (dis_kp1 - dis_km1) * inv_rmsinthe2dphi;

                // ---- production tensor P_k ----
                // nue_t is recomputed here from the current tke/dis rather than read from
                // m.nue, which lags one call behind and is still 0 on the first init() call.
                double cnue;
                {
                    const double dis_here = std::max(m.dis.x[i][j][k], dis_min_l);
                    const double tke_here = std::max(m.tke.x[i][j][k], 0.0);
                    const double nue_max  = nue_max_phys() / (m.u_0 * L_atm_m());
                    cnue = (turb_model == k_epsilon)
                         ? C_nue_l * tke_here * tke_here / dis_here
                         : tke_here / dis_here;
                    cnue = std::min(cnue, nue_max);
                }
                // ---- velocity gradient tensor in the orthonormal spherical basis ----
                // The nine components carry curvature terms that the bare derivatives do not.
                // dudthe and dudphi above already include their 1/r and 1/(r sin) metric
                // factors, so what is missing is only the extra pieces:
                //   G_tt = (1/r) dv/dthe + u/r                G_rt = (1/r) du/dthe - v/r
                //   G_pp = (1/r sin) dw/dphi + u/r + v cot/r  G_rp = (1/r sin) du/dphi - w/r
                //   G_tp = (1/r sin) dv/dphi - w cot/r
                // and the divergence gains 2u/r + v cot/r. Their size is set entirely by the
                // metric radius: with ATJUP_METRIC_RADIUS they are O(1/500) of the radial
                // derivatives, without it they are comparable to the horizontal ones — which is
                // why omitting them was defensible in one geometry and not in the other.
                // ATJUP_TURB_CURV=0 drops them again for A/B.
                static const bool curv = [](){
                    const char* e = ATPhys::env_for(Planet::planet_tag(), "TURB_CURV"); return !e || atoi(e) != 0; }();
                const double u_r = curv ? m.u.x[i][j][k] * inv_rm : 0.0;
                const double v_r = curv ? m.v.x[i][j][k] * inv_rm : 0.0;
                const double w_r = curv ? m.w.x[i][j][k] * inv_rm : 0.0;
                const double v_cot_r = v_r * cotanthe;
                const double w_cot_r = w_r * cotanthe;

                const double g_rr = dudr;
                const double g_rt = dudthe - v_r;
                const double g_rp = dudphi - w_r;
                const double g_tr = dvdr;
                const double g_tt = dvdthe + u_r;
                const double g_tp = dvdphi - w_cot_r;
                const double g_pr = dwdr;
                const double g_pt = dwdthe;
                const double g_pp = dwdphi + u_r + v_cot_r;

                const double der = 0.66667 * (g_rr + g_tt + g_pp);

                m.prod.x[i][j][k] = std::max(0.0,
                      (cnue * (2.0 * g_rr - der) - 0.66667 * m.tke.x[i][j][k]) * g_rr
                    + (cnue * (g_rt + g_tr))                                   * g_rt
                    + (cnue * (g_rp + g_pr))                                   * g_rp
                    + (cnue * (2.0 * g_tt - der) - 0.66667 * m.tke.x[i][j][k]) * g_tt
                    + (cnue * (g_tr + g_rt))                                   * g_tr
                    + (cnue * (g_tp + g_pt))                                   * g_tp
                    + (cnue * (2.0 * g_pp - der) - 0.66667 * m.tke.x[i][j][k]) * g_pp
                    + (cnue * (g_pr + g_rp))                                   * g_pr
                    + (cnue * (g_pt + g_tp))                                   * g_pt);

                // ---- vorticity magnitude Omega = |curl u| ----
                // Antisymmetric part of the SAME tensor, so it carries the curvature too.
                const double W12 = g_rt - g_tr;
                const double W13 = g_rp - g_pr;
                const double W23 = g_tp - g_pt;
                const double Omega = std::sqrt(W12*W12 + W13*W13 + W23*W23);

                if(turb_model == k_epsilon){
                    compute_k_epsilon(i, j, k, y_mount,
                        dudr, dvdr, dwdr, dudthe, dvdthe, dwdthe, dudphi, dvdphi, dwdphi,
                        Omega);
                }
                else if(turb_model == k_omega){
                    compute_k_omega(i, j, k,
                        dtkedr, ddisdr, dtkedthe, ddisdthe, dtkedphi, ddisdphi,
                        Omega, Omega);
                }
                else if(turb_model == k_omega_SST){
                    compute_k_omega_SST(i, j, k, y_mount, rm, sinthe,
                        dtkedr, ddisdr, dtkedthe, ddisdthe, dtkedphi, ddisdphi,
                        Omega, Omega);
                }

                // Caps: nue_max is a 1000 m2/s physical cap made dimensionless;
                // the source caps bound the P-Y balance and stop the 1/sin^2(theta)
                // amplification in the cross-diffusion term blowing up at the poles.
                const double nue_max = nue_max_phys() / (m.u_0 * L_atm_m());
                m.tke.x[i][j][k] = std::max(k_bg_nd,   m.tke.x[i][j][k]);
                m.dis.x[i][j][k] = std::max(dis_min_l, m.dis.x[i][j][k]);
                m.nue.x[i][j][k] = clamp3(m.nue.x[i][j][k], 0.0, nue_max);

                const double dis_src_max = 20.0 * m.dis.x[i][j][k] * m.dis.x[i][j][k];
                m.dis_source.x[i][j][k] = clamp3(m.dis_source.x[i][j][k], -dis_src_max, dis_src_max);
                const double tke_src_max = 20.0 * m.tke.x[i][j][k] * m.dis.x[i][j][k];
                m.tke_source.x[i][j][k] = clamp3(m.tke_source.x[i][j][k], -tke_src_max, tke_src_max);
            }
        }
    }
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::compute_k_epsilon(int i, int j, int k, double y_mount,
    double dudr, double dvdr, double dwdr,
    double dudthe, double dvdthe, double dwdthe,
    double dudphi, double dvdphi, double dwdphi,
    double /*Omega*/)
{
    const double C_eps_1 = 1.35;
    const double C_eps_2 = 1.80;
    const double C_nue_l = C_nue, dis_min_l = dis_min;
    const double L_m = L_atm_m();

    // Physical wall distance [m], floored at one grid length so y_star >= 1 keeps
    // the 1/y_star^2 in the Chien D term bounded over mountain tops.
    const double y_phys = std::max(height_m(i) - y_mount, L_m);
    const double y_star = y_phys / L_m;

    m.nue.x[i][j][k] = C_nue_l * m.tke.x[i][j][k] * m.tke.x[i][j][k]
                      / std::max(m.dis.x[i][j][k], dis_min_l);

    const double d_plus = y_phys * m.vel_star.y[j][k] / nue_air;
    const double f_nue  = 1.0 - std::exp(-0.0115 * d_plus);

    const double Re_T = m.tke.x[i][j][k] * m.tke.x[i][j][k] * m.u_0 * L_m
                        / std::max(m.dis.x[i][j][k] * nue_air, dis_min_l * nue_air);
    const double f_2  = 1.0 - 0.4 / 1.8 * std::exp(-Re_T * Re_T / 36.0);

    m.nue.x[i][j][k] *= f_nue;

    const double y_star2 = y_star * y_star;
    const double L_k = -2.0 * m.tke.x[i][j][k] / y_star2;
    const double L_w = -2.0 * m.dis.x[i][j][k] / y_star2 * std::exp(-d_plus / 2.0);

    const double P_k = m.prod.x[i][j][k];
    const double Y_k = m.dis.x[i][j][k];
    const double tke_safe = std::max(m.tke.x[i][j][k], dis_min_l);
    const double P_w = C_eps_1 * m.dis.x[i][j][k] / tke_safe * P_k;
    const double Y_w = C_eps_2 * f_2 * m.dis.x[i][j][k] * m.dis.x[i][j][k] / tke_safe;

    m.tke_source.x[i][j][k] = P_k - Y_k + L_k / re_turb;
    m.dis_source.x[i][j][k] = P_w - Y_w + L_w / re_turb;
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::compute_k_omega(int i, int j, int k,
    double dtkedr, double ddisdr,
    double dtkedthe, double ddisdthe,
    double dtkedphi, double ddisdphi,
    double Omega_mag, double /*W_unused*/)
{
    const double bet_star = 0.09;
    const double gam      = 0.52;
    const double C_lim    = 0.875;
    const double bet_0    = 0.0708;

    const double rm       = m.rad.z[i];
    double sinthe         = std::sin(m.the.z[j]);
    if(sinthe == 0.0) sinthe = 1.0e-5;
    const double rmsinthe = rm * sinthe;
    const double exp_rm   = 1.0 / (rm + 1.0);

    const bool neumann_bot = (i > 0        && is_land(i-1, j, k));
    const bool neumann_jm1 = (j > 0        && is_land(i, j-1, k));
    const bool neumann_jp1 = (j < m.jm - 1 && is_land(i, j+1, k));
    const bool neumann_km1 = (k > 0        && is_land(i, j, k-1));
    const bool neumann_kp1 = (k < m.km - 1 && is_land(i, j, k+1));

    const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
    const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
    const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
    const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
    const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
    const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
    const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
    const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

    const double dtkedr_neu = neumann_bot
        ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm : dtkedr;
    const double ddisdr_neu = neumann_bot
        ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm : ddisdr;
    const double dtkedthe_neu = (tke_jp1 - tke_jm1) / (rm * 2.0 * m.dthe);
    const double ddisdthe_neu = (dis_jp1 - dis_jm1) / (rm * 2.0 * m.dthe);
    const double dtkedphi_neu = (tke_kp1 - tke_km1) / (rmsinthe * 2.0 * m.dphi);
    const double ddisdphi_neu = (dis_kp1 - dis_km1) / (rmsinthe * 2.0 * m.dphi);

    const double dudr   = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) / (2.0 * m.dr) / (rm + 1.0);
    const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) / (rm * 2.0 * m.dthe);
    const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);

    const double dudthe_u = (m.u.x[i][j+1][k] - m.u.x[i][j-1][k]) / (rm       * 2.0 * m.dthe);
    const double dudphi_u = (m.u.x[i][j][k+1] - m.u.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);
    const double dvdr_u   = (m.v.x[i+1][j][k] - m.v.x[i-1][j][k]) / (2.0 * m.dr) / (rm + 1.0);
    const double dwdr_u   = (m.w.x[i+1][j][k] - m.w.x[i-1][j][k]) / (2.0 * m.dr) / (rm + 1.0);
    const double dvdphi_u = (m.v.x[i][j][k+1] - m.v.x[i][j][k-1]) / (rmsinthe * 2.0 * m.dphi);
    const double dwdthe_u = (m.w.x[i][j+1][k] - m.w.x[i][j-1][k]) / (rm       * 2.0 * m.dthe);

    const double S11 = dudr, S22 = dvdthe, S33 = dwdphi;
    const double S12 = 0.5 * (dudthe_u + dvdr_u);
    const double S13 = 0.5 * (dudphi_u + dwdr_u);
    const double S23 = 0.5 * (dvdphi_u + dwdthe_u);
    const double S_mag = std::sqrt(2.0 * (S11*S11 + S22*S22 + S33*S33
                                        + 2.0*(S12*S12 + S13*S13 + S23*S23)));

    const double w_hat = std::max(m.dis.x[i][j][k], C_lim * S_mag / std::sqrt(bet_star));

    double sig_d = dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu;
    sig_d = (sig_d <= 0.0) ? 0.0 : 0.125;

    // Wilcox (2006) vorticity-strain correction chi_omega
    const double chi_w = std::fabs(Omega_mag * Omega_mag * S_mag
        / std::pow(bet_star * std::max(m.dis.x[i][j][k], 1.0e-20), 3));
    const double f_bet  = (1.0 + 85.0 * chi_w) / (1.0 + 100.0 * chi_w);
    const double bet_wc = bet_0 * f_bet;

    m.nue.x[i][j][k] = m.tke.x[i][j][k] / std::max(w_hat, 1.0e-20);

    const double P_k = std::min(m.prod.x[i][j][k],
        20.0 * bet_star * m.tke.x[i][j][k] * m.dis.x[i][j][k]);
    const double Y_k = bet_star * m.tke.x[i][j][k] * m.dis.x[i][j][k];
    const double P_w = gam * m.dis.x[i][j][k] / std::max(m.tke.x[i][j][k], 1.0e-20) * P_k;
    const double Y_w = bet_wc * m.dis.x[i][j][k] * m.dis.x[i][j][k];
    const double D_w = sig_d / std::max(m.dis.x[i][j][k], 1.0e-20)
        * (dtkedr_neu * ddisdr_neu + dtkedthe_neu * ddisdthe_neu + dtkedphi_neu * ddisdphi_neu);

    // Only the production-destruction balance is stored; turbulent diffusion of tke
    // belongs to the RHS and must not be double-counted here.
    m.tke_source.x[i][j][k] = P_k - Y_k;
    m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
}

// -----------------------------------------------------------------------
template<class Planet>
inline void Turbulence<Planet>::compute_k_omega_SST(int i, int j, int k,
    double y_mount, double rm, double sinthe,
    double dtkedr, double ddisdr,
    double dtkedthe, double ddisdthe,
    double dtkedphi, double ddisdphi,
    double Omega, double /*W*/)
{
    // Menter 1994 ABL constants
    const double a1       = 0.31;
    const double bet_star = 0.09;      // beta* destruction coefficient, NOT C_mu = 0.028
    const double bet1     = 0.0333;
    const double bet2     = 0.0368;
    const double gam1     = 0.413;
    const double gam2     = 0.2;
    const double sig_w2   = 1.168;

    const double L_m = L_atm_m();
    const double nue_air_nd = nue_air / (m.u_0 * L_m);

    // Same one-grid-layer floor as init_fields: the blending functions divide by y_star, so a
    // 1e-6 m wall distance in the surface layer makes arg1/arg2 diverge.
    const double dz_layer = std::max(1.0, m.layer_thickness_m(std::max(i - 1, 0)));
    const double y_phys = std::max(height_m(i) - y_mount, dz_layer);
    const double y_star = y_phys / L_m;

    const double rmsinthe = rm * sinthe;
    const double exp_rm   = 1.0 / (rm + 1.0);

    const bool neumann_bot = (i > 0        && is_land(i-1, j, k));
    const bool neumann_jm1 = (j > 0        && is_land(i, j-1, k));
    const bool neumann_jp1 = (j < m.jm - 1 && is_land(i, j+1, k));
    const bool neumann_km1 = (k > 0        && is_land(i, j, k-1));
    const bool neumann_kp1 = (k < m.km - 1 && is_land(i, j, k+1));

    const double tke_jm1 = neumann_jm1 ? m.tke.x[i][j][k] : m.tke.x[i][j-1][k];
    const double tke_jp1 = neumann_jp1 ? m.tke.x[i][j][k] : m.tke.x[i][j+1][k];
    const double dis_jm1 = neumann_jm1 ? m.dis.x[i][j][k] : m.dis.x[i][j-1][k];
    const double dis_jp1 = neumann_jp1 ? m.dis.x[i][j][k] : m.dis.x[i][j+1][k];
    const double tke_km1 = neumann_km1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k-1];
    const double tke_kp1 = neumann_kp1 ? m.tke.x[i][j][k] : m.tke.x[i][j][k+1];
    const double dis_km1 = neumann_km1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k-1];
    const double dis_kp1 = neumann_kp1 ? m.dis.x[i][j][k] : m.dis.x[i][j][k+1];

    const double dtkedr_neu = neumann_bot
        ? (m.tke.x[i+1][j][k] - m.tke.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm : dtkedr;
    const double ddisdr_neu = neumann_bot
        ? (m.dis.x[i+1][j][k] - m.dis.x[i][j][k]) * (1.0 / (2.0 * m.dr)) * exp_rm : ddisdr;
    const double dtkedthe_neu = (tke_jp1 - tke_jm1) / (rm       * 2.0 * m.dthe);
    const double ddisdthe_neu = (dis_jp1 - dis_jm1) / (rm       * 2.0 * m.dthe);
    const double dtkedphi_neu = (tke_kp1 - tke_km1) / (rmsinthe * 2.0 * m.dphi);
    const double ddisdphi_neu = (dis_kp1 - dis_km1) / (rmsinthe * 2.0 * m.dphi);

    const double grad_dot = dtkedr_neu * ddisdr_neu
                          + dtkedthe_neu * ddisdthe_neu
                          + dtkedphi_neu * ddisdphi_neu;

    const double CD_kw = std::max(
        2.0 * sig_w2 / std::max(m.dis.x[i][j][k], 1.0e-20) * grad_dot, 1.0e-20);

    const double dis_safe = std::max(m.dis.x[i][j][k], 1.0e-20);
    const double arg1 = std::min(
        std::max(std::sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
                 500.0 * nue_air_nd / (y_star * y_star * dis_safe)),
        4.0 * sig_w2 * m.tke.x[i][j][k] / (CD_kw * y_star * y_star));
    const double arg2 = std::max(
        2.0 * std::sqrt(m.tke.x[i][j][k]) / (bet_star * dis_safe * y_star),
        500.0 * nue_air_nd / (y_star * y_star * dis_safe));

    const double F1 = std::tanh(std::pow(arg1, 4));
    const double F2 = std::tanh(std::pow(arg2, 2));

    m.nue.x[i][j][k] = a1 * m.tke.x[i][j][k] / std::max(a1 * dis_safe, Omega * F2);

    const double tke_safe = std::max(m.tke.x[i][j][k], 1.0e-20);
    const double P_k = std::min(m.prod.x[i][j][k], 20.0 * bet_star * tke_safe * dis_safe);
    const double Y_k = bet_star * tke_safe * dis_safe;
    const double P_w = blend(gam1, gam2, F1) * P_k * dis_safe / tke_safe;
    const double Y_w = blend(bet1, bet2, F1) * dis_safe * dis_safe;
    const double D_w = 2.0 * (1.0 - F1) * sig_w2 / dis_safe * grad_dot;

    m.tke_source.x[i][j][k] = P_k - Y_k;
    m.dis_source.x[i][j][k] = P_w - Y_w + D_w;
}

// -----------------------------------------------------------------------
// HISTORY: this file used to carry an advance() that integrated k* and dis* from the source
// terms with a Patankar splitting, because ATJUP's Runge-Kutta had no turbulence equations and
// the sources computed above would otherwise have been discarded every iteration. That local
// update could not transport k*/dis*: they evolved purely pointwise, apart from the SST
// cross-diffusion term already inside dis_source. Both equations now live in RHS_Jup_Turb.cpp
// with their advection and turbulent-diffusion terms and are integrated by
// RungeKutta_Jup_Turb.cpp, which is ATOM's arrangement, so advance() has been removed.

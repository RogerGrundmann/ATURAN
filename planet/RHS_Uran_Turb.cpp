/*
 * Atmosphere General Circulation Modell(ATURAN) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with additional transport equations for CH4, H2O, H2S, NH3, NH4SH concentrations
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include "cUranusModel.h"

using namespace std;


void cUranusModel::RHSUran(int i, int j, int k, const CellGeometry& geo){

    // All geometric quantities come from the precomputed struct —
    // NO sin(), cos(), division, or reciprocal computation here.
    const double rm                   = geo.rm;
    const double sinthe               = geo.sinthe;
    const double costhe               = geo.costhe;
    const double cotanthe             = geo.cotanthe;
    const double inv_rm               = geo.inv_rm;
    const double inv_rm2              = geo.inv_rm2;
    const double inv_rmsinthe         = geo.inv_rmsinthe;
    const double inv_rm2sinthe        = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2       = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dr   = geo.inv_2dr;
    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dr2   = geo.inv_dr2;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;

    // Cache local cell values
    const double u_ijk = u.x[i][j][k];
    const double v_ijk = v.x[i][j][k];
    const double w_ijk = w.x[i][j][k];

    // ---- First-order derivative storage ----
    double dudr, dvdr, dwdr, dtdr, dpdr;
    double dch4dr, dch4cdr, dch4idr;
    double dh2odr, dh2ocdr, dh2oidr;
    double dh2sdr, dh2scdr, dh2sidr;
    double dnh3dr, dnh3cdr, dnh3idr, dnh4shdr;

    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe;
    double dch4dthe, dch4cdthe, dch4idthe;
    double dh2odthe, dh2ocdthe, dh2oidthe;
    double dh2sdthe, dh2scdthe, dh2sidthe;
    double dnh3dthe, dnh3cdthe, dnh3idthe, dnh4shdthe;

    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi;
    double dch4dphi, dch4cdphi, dch4idphi;
    double dh2odphi, dh2ocdphi, dh2oidphi;
    double dh2sdphi, dh2scdphi, dh2sidphi;
    double dnh3dphi, dnh3cdphi, dnh3idphi, dnh4shdphi;

    // ---- Second-order derivative storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2;
    double d2ch4dr2, d2ch4cdr2, d2ch4idr2;
    double d2h2odr2, d2h2ocdr2, d2h2oidr2;
    double d2h2sdr2, d2h2scdr2, d2h2sidr2;
    double d2nh3dr2, d2nh3cdr2, d2nh3idr2, d2nh4shdr2;

    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2;
    double d2ch4dthe2, d2ch4cdthe2, d2ch4idthe2;
    double d2h2odthe2, d2h2ocdthe2, d2h2oidthe2;
    double d2h2sdthe2, d2h2scdthe2, d2h2sidthe2;
    double d2nh3dthe2, d2nh3cdthe2, d2nh3idthe2, d2nh4shdthe2;

    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2;
    double d2ch4dphi2, d2ch4cdphi2, d2ch4idphi2;
    double d2h2odphi2, d2h2ocdphi2, d2h2oidphi2;
    double d2h2sdphi2, d2h2scdphi2, d2h2sidphi2;
    double d2nh3dphi2, d2nh3cdphi2, d2nh3idphi2, d2nh4shdphi2;


    // ===== R-direction derivatives (central differences) =====
    #define COMPUTE_DR(FIELD, d1, d2) \
        d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dr; \
        d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dr2;

    COMPUTE_DR(u,         dudr,    d2udr2)
    COMPUTE_DR(v,         dvdr,    d2vdr2)
    COMPUTE_DR(w,         dwdr,    d2wdr2)
    COMPUTE_DR(t,         dtdr,    d2tdr2)
    COMPUTE_DR(ch4,       dch4dr,  d2ch4dr2)
    COMPUTE_DR(ch4_cloud, dch4cdr, d2ch4cdr2)
    COMPUTE_DR(ch4_ice,   dch4idr, d2ch4idr2)
    COMPUTE_DR(h2o,       dh2odr,  d2h2odr2)
    COMPUTE_DR(h2o_cloud, dh2ocdr, d2h2ocdr2)
    COMPUTE_DR(h2o_ice,   dh2oidr, d2h2oidr2)
    COMPUTE_DR(h2s,       dh2sdr,  d2h2sdr2)
    COMPUTE_DR(h2s_cloud, dh2scdr, d2h2scdr2)
    COMPUTE_DR(h2s_ice,   dh2sidr, d2h2sidr2)
    COMPUTE_DR(nh3,       dnh3dr,  d2nh3dr2)
    COMPUTE_DR(nh3_cloud, dnh3cdr, d2nh3cdr2)
    COMPUTE_DR(nh3_ice,   dnh3idr, d2nh3idr2)
    COMPUTE_DR(nh4sh,     dnh4shdr,d2nh4shdr2)
    double dtkedr, ddisdr, d2tkedr2, d2disdr2;
    COMPUTE_DR(tke,       dtkedr,  d2tkedr2)
    COMPUTE_DR(dis,       ddisdr,  d2disdr2)
    dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr;
    #undef COMPUTE_DR


    // ===== Theta-direction derivatives (central differences) =====
    #define COMPUTE_DTHE(FIELD, d1, d2) \
        d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dthe; \
        d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dthe2;

    COMPUTE_DTHE(u,         dudthe,    d2udthe2)
    COMPUTE_DTHE(v,         dvdthe,    d2vdthe2)
    COMPUTE_DTHE(w,         dwdthe,    d2wdthe2)
    COMPUTE_DTHE(t,         dtdthe,    d2tdthe2)
    COMPUTE_DTHE(ch4,       dch4dthe,  d2ch4dthe2)
    COMPUTE_DTHE(ch4_cloud, dch4cdthe, d2ch4cdthe2)
    COMPUTE_DTHE(ch4_ice,   dch4idthe, d2ch4idthe2)
    COMPUTE_DTHE(h2o,       dh2odthe,  d2h2odthe2)
    COMPUTE_DTHE(h2o_cloud, dh2ocdthe, d2h2ocdthe2)
    COMPUTE_DTHE(h2o_ice,   dh2oidthe, d2h2oidthe2)
    COMPUTE_DTHE(h2s,       dh2sdthe,  d2h2sdthe2)
    COMPUTE_DTHE(h2s_cloud, dh2scdthe, d2h2scdthe2)
    COMPUTE_DTHE(h2s_ice,   dh2sidthe, d2h2sidthe2)
    COMPUTE_DTHE(nh3,       dnh3dthe,  d2nh3dthe2)
    COMPUTE_DTHE(nh3_cloud, dnh3cdthe, d2nh3cdthe2)
    COMPUTE_DTHE(nh3_ice,   dnh3idthe, d2nh3idthe2)
    COMPUTE_DTHE(nh4sh,     dnh4shdthe,d2nh4shdthe2)
    double dtkedthe, ddisdthe, d2tkedthe2, d2disdthe2;
    COMPUTE_DTHE(tke,       dtkedthe, d2tkedthe2)
    COMPUTE_DTHE(dis,       ddisdthe, d2disdthe2)
    dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
    #undef COMPUTE_DTHE


    // ===== Phi-direction derivatives (central differences) =====
    #define COMPUTE_DPHI(FIELD, d1, d2) \
        d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dphi; \
        d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dphi2;

    COMPUTE_DPHI(u,         dudphi,    d2udphi2)
    COMPUTE_DPHI(v,         dvdphi,    d2vdphi2)
    COMPUTE_DPHI(w,         dwdphi,    d2wdphi2)
    COMPUTE_DPHI(t,         dtdphi,    d2tdphi2)
    COMPUTE_DPHI(ch4,       dch4dphi,  d2ch4dphi2)
    COMPUTE_DPHI(ch4_cloud, dch4cdphi, d2ch4cdphi2)
    COMPUTE_DPHI(ch4_ice,   dch4idphi, d2ch4idphi2)
    COMPUTE_DPHI(h2o,       dh2odphi,  d2h2odphi2)
    COMPUTE_DPHI(h2o_cloud, dh2ocdphi, d2h2ocdphi2)
    COMPUTE_DPHI(h2o_ice,   dh2oidphi, d2h2oidphi2)
    COMPUTE_DPHI(h2s,       dh2sdphi,  d2h2sdphi2)
    COMPUTE_DPHI(h2s_cloud, dh2scdphi, d2h2scdphi2)
    COMPUTE_DPHI(h2s_ice,   dh2sidphi, d2h2sidphi2)
    COMPUTE_DPHI(nh3,       dnh3dphi,  d2nh3dphi2)
    COMPUTE_DPHI(nh3_cloud, dnh3cdphi, d2nh3cdphi2)
    COMPUTE_DPHI(nh3_ice,   dnh3idphi, d2nh3idphi2)
    COMPUTE_DPHI(nh4sh,     dnh4shdphi,d2nh4shdphi2)
    double dtkedphi, ddisdphi, d2tkedphi2, d2disdphi2;
    COMPUTE_DPHI(tke,       dtkedphi, d2tkedphi2)
    COMPUTE_DPHI(dis,       ddisdphi, d2disdphi2)
    dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
    #undef COMPUTE_DPHI


    // ===== Coriolis and centrifugal forces =====
    // ===== TWO OF THE THREE CORIOLIS COMPONENTS DEFLECTED TO THE LEFT =====
    //
    // Ported from ATJUP ec698e1, which ATURAN never received. With theta the COLATITUDE and v the
    // theta-component (so +v is southward), -2*Omega x u is
    //
    //     a_r = +2*Omega*sin*w    a_theta = +2*Omega*cos*w    a_phi = -2*Omega*(cos*v + sin*u)
    //
    // and all three enter rhs_* through a minus, so each variable below must hold -a.
    //
    // Coriolis_rad already did. Coriolis_the held +2*Omega*cos*w and so delivered -a_theta: an
    // eastward wind in the northern hemisphere was turned NORTHWARD. The physical check is
    // unambiguous — eastward flow at northern mid-latitudes must go RIGHT, which is south, which
    // is +e_theta. Coriolis_phi held +2*Omega*(-cos*v + sin*u): the sin*u half was right and the
    // cos*v half carried the wrong sign, so southward flow in the north was turned EAST instead of
    // west, again to the left. That mixture is why it survived — Coriolis_phi was not plus or
    // minus any consistent expression, so neither sign looked obviously right.
    //
    // THIS IS THE COMPONENT THAT CAN MOVE THE CIRCULATION. Unlike the centrifugal force corrected
    // in 4201957, Coriolis is velocity-dependent and NOT curl-free, so the pressure projection
    // cannot absorb it — and it is Coriolis that closes an overturning cell.
    //
    // ATURAN_CORIOLIS_LEGACY=1 restores the two pre-port signs.
    //
    // NOT FIXED HERE, and it is a separate defect in the same term: Coriolis_phi enters rhs_w
    // WITHOUT the scale_Cor = L_atm/u_0 factor that rhs_u and rhs_v apply to the other two
    // components, so the zonal component is 3600x weaker than its siblings. That belongs to the
    // unit-system port (ATJUP af0446d, which carries one nd_cor on all three) and folding a 3600x
    // magnitude change into a sign fix would make neither measurable.
    static const bool cor_legacy = [](){
        const char* e = getenv("ATURAN_CORIOLIS_LEGACY"); return e && atoi(e) != 0; }();
    double Coriolis_rad  = -2.0 * omega * sinthe * w_ijk;
    double Coriolis_the  = (cor_legacy ? +2.0 : -2.0) * omega * costhe * w_ijk;
    double Coriolis_phi  = +2.0 * omega * ((cor_legacy ? -1.0 : +1.0) * costhe * v_ijk
                                           + sinthe * u_ijk);

    // ===== THE CENTRIFUGAL FORCE POINTS AWAY FROM THE ROTATION AXIS =====
    //
    // Ported from ATJUP 8649675, which ATURAN never received. Centrifugal acceleration is
    // Omega^2 * s * s_hat with s = r*sin(theta) the distance from the ROTATION AXIS and
    // s_hat = sin(theta)*e_r + cos(theta)*e_theta pointing away from it:
    //
    //     a_r     = +Omega^2 * r * sin^2(theta)
    //     a_theta = +Omega^2 * r * sin(theta) * cos(theta)
    //
    // All three parts were wrong here, in the same three ways ATJUP records. The radial part
    // carried NO sin^2, so it had full strength at the poles where it must vanish. The meridional
    // part had |sin| where sin*cos belongs, which is neither the right magnitude nor
    // equator-directed — and the absolute value destroyed the hemispheric antisymmetry, so both
    // hemispheres were pushed the same way. And both entered rhs_* through a MINUS, pointing the
    // force TOWARD the axis instead of away from it; they now enter through a plus, as ATJUP's do.
    //
    // The equator-directed meridional part only works because cos(theta) already changes sign at
    // the equator here — ATJUP needed 53e75b2 first for that, and ATURAN's costhe_tbl is built as
    // cos(the.z[j]) with no fabs, so that prerequisite was already satisfied.
    //
    // STILL WRONG, and deliberately left for step 3: `sinthe` is the CLAMPED metric value, floored
    // at sinthe_min() = 0.4, so sin^2 reads 0.16 at the poles rather than 0. That is a metric guard
    // leaking into a body force. ATJUP fixed it in 7782207 by deriving sinthe_true from costhe.
    // Fixing it here would conflate two ports in one measurement.
    //
    // ATURAN_CENT_LEGACY=1 restores the pre-port form AND the pre-port entry sign together, since
    // neither is meaningful without the other.
    static const bool cent_legacy = [](){
        const char* e = getenv("ATURAN_CENT_LEGACY"); return e && atoi(e) != 0; }();
    double centrifugal_rad, centrifugal_the;
    if(cent_legacy){
        centrifugal_rad = omega * omega * rm;
        centrifugal_the = omega * omega * rm * fabs(sinthe);
    } else {
        centrifugal_rad = omega * omega * rm * sinthe * sinthe;
        centrifugal_the = omega * omega * rm * sinthe * costhe;
    }

    double coeff_energy_p = u_0 * u_0 / (cp_mix * t_ref);

    // Non-dimensional scaling factors (velocities = v/u_0, lengths = r/L_atm)
    const double scale_p   = 1.0e5 / (r_mix * u_0 * u_0);
    const double scale_Cor = L_atm / u_0;
    const double scale_cen = L_atm * L_atm / (u_0 * u_0);


    // ===== Transport terms (advection) =====
    double v_invrm = v_ijk * inv_rm;
    double w_invrs = w_ijk * inv_rmsinthe;

    double pressure_t = coeff_energy_p * scale_p
        * (u_ijk * dpdr + v_invrm * dpdthe + w_invrs * dpdphi);

    double transport_t = u_ijk * dtdr + v_invrm * dtdthe + w_invrs * dtdphi;

    double transport_u = u_ijk * dudr + v_invrm * dudthe + w_invrs * dudphi
        - (v_ijk * v_ijk + w_ijk * w_ijk) * inv_rm;
    double transport_v = u_ijk * dvdr + v_invrm * dvdthe + w_invrs * dvdphi
        + (u_ijk * v_ijk - w_ijk * w_ijk * cotanthe) * inv_rm;
    double transport_w = u_ijk * dwdr + v_invrm * dwdthe + w_invrs * dwdphi
        + (w_ijk * u_ijk + v_ijk * w_ijk * cotanthe) * inv_rm;

    double transport_ch4       = u_ijk * dch4dr  + v_invrm * dch4dthe  + w_invrs * dch4dphi;
    double transport_ch4_cloud = u_ijk * dch4cdr + v_invrm * dch4cdthe + w_invrs * dch4cdphi;
    double transport_ch4_ice   = u_ijk * dch4idr + v_invrm * dch4idthe + w_invrs * dch4idphi;

    double transport_h2o       = u_ijk * dh2odr  + v_invrm * dh2odthe  + w_invrs * dh2odphi;
    double transport_h2o_cloud = u_ijk * dh2ocdr + v_invrm * dh2ocdthe + w_invrs * dh2ocdphi;
    double transport_h2o_ice   = u_ijk * dh2oidr + v_invrm * dh2oidthe + w_invrs * dh2oidphi;

    double transport_h2s       = u_ijk * dh2sdr  + v_invrm * dh2sdthe  + w_invrs * dh2sdphi;
    double transport_h2s_cloud = u_ijk * dh2scdr + v_invrm * dh2scdthe + w_invrs * dh2scdphi;
    double transport_h2s_ice   = u_ijk * dh2sidr + v_invrm * dh2sidthe + w_invrs * dh2sidphi;

    double transport_nh3       = u_ijk * dnh3dr  + v_invrm * dnh3dthe  + w_invrs * dnh3dphi;
    double transport_nh3_cloud = u_ijk * dnh3cdr + v_invrm * dnh3cdthe + w_invrs * dnh3cdphi;
    double transport_nh3_ice   = u_ijk * dnh3idr + v_invrm * dnh3idthe + w_invrs * dnh3idphi;

    double transport_nh4sh     = u_ijk * dnh4shdr + v_invrm * dnh4shdthe + w_invrs * dnh4shdphi;


    // ===== Diffusion terms =====
    double two_inv_rm = 2.0 * inv_rm;
    double v_metric   = (1.0 + costhe / geo.sinthe2) * inv_rm2;

    double diffusion_t = d2tdr2 + dtdr * two_inv_rm + d2tdthe2 * inv_rm2
        + dtdthe * costhe_inv_rm2sinthe + d2tdphi2 * inv_rm2sinthe2;

    double diffusion_u = d2udr2 + 2.0 * u_ijk * inv_rm2 + d2udthe2 * inv_rm2
        + 4.0 * dudr * inv_rm + dudthe * costhe_inv_rm2sinthe
        + d2udphi2 * inv_rm2sinthe2;
    double diffusion_v = d2vdr2 + dvdr * two_inv_rm + d2vdthe2 * inv_rm2
        + dvdthe * costhe_inv_rm2sinthe - v_metric * v_ijk
        + d2vdphi2 * inv_rm2sinthe2
        + 2.0 * dudthe * inv_rm2
        - dwdphi * 2.0 * costhe * inv_rm2sinthe2;
    double diffusion_w = d2wdr2 + dwdr * two_inv_rm + d2wdthe2 * inv_rm2
        + dwdthe * costhe_inv_rm2sinthe - v_metric * w_ijk
        + d2wdphi2 * inv_rm2sinthe2
        + 2.0 * dudphi * inv_rm2sinthe
        + dvdphi * 2.0 * costhe * inv_rm2sinthe2;

    double diffusion_ch4 = d2ch4dr2 + dch4dr * two_inv_rm + d2ch4dthe2 * inv_rm2
        + dch4dthe * costhe_inv_rm2sinthe + d2ch4dphi2 * inv_rm2sinthe2;
    double diffusion_ch4_cloud = d2ch4cdr2 + dch4cdr * two_inv_rm + d2ch4cdthe2 * inv_rm2
        + dch4cdthe * costhe_inv_rm2sinthe + d2ch4cdphi2 * inv_rm2sinthe2;
    double diffusion_ch4_ice = d2ch4idr2 + dch4idr * two_inv_rm + d2ch4idthe2 * inv_rm2
        + dch4idthe * costhe_inv_rm2sinthe + d2ch4idphi2 * inv_rm2sinthe2;

    double diffusion_h2o = d2h2odr2 + dh2odr * two_inv_rm + d2h2odthe2 * inv_rm2
        + dh2odthe * costhe_inv_rm2sinthe + d2h2odphi2 * inv_rm2sinthe2;
    double diffusion_h2o_cloud = d2h2ocdr2 + dh2ocdr * two_inv_rm + d2h2ocdthe2 * inv_rm2
        + dh2ocdthe * costhe_inv_rm2sinthe + d2h2ocdphi2 * inv_rm2sinthe2;
    double diffusion_h2o_ice = d2h2oidr2 + dh2oidr * two_inv_rm + d2h2oidthe2 * inv_rm2
        + dh2oidthe * costhe_inv_rm2sinthe + d2h2oidphi2 * inv_rm2sinthe2;

    double diffusion_h2s = d2h2sdr2 + dh2sdr * two_inv_rm + d2h2sdthe2 * inv_rm2
        + dh2sdthe * costhe_inv_rm2sinthe + d2h2sdphi2 * inv_rm2sinthe2;
    double diffusion_h2s_cloud = d2h2scdr2 + dh2scdr * two_inv_rm + d2h2scdthe2 * inv_rm2
        + dh2scdthe * costhe_inv_rm2sinthe + d2h2scdphi2 * inv_rm2sinthe2;
    double diffusion_h2s_ice = d2h2sidr2 + dh2sidr * two_inv_rm + d2h2sidthe2 * inv_rm2
        + dh2sidthe * costhe_inv_rm2sinthe + d2h2sidphi2 * inv_rm2sinthe2;

    double diffusion_nh3 = d2nh3dr2 + dnh3dr * two_inv_rm + d2nh3dthe2 * inv_rm2
        + dnh3dthe * costhe_inv_rm2sinthe + d2nh3dphi2 * inv_rm2sinthe2;
    double diffusion_nh3_cloud = d2nh3cdr2 + dnh3cdr * two_inv_rm + d2nh3cdthe2 * inv_rm2
        + dnh3cdthe * costhe_inv_rm2sinthe + d2nh3cdphi2 * inv_rm2sinthe2;
    double diffusion_nh3_ice = d2nh3idr2 + dnh3idr * two_inv_rm + d2nh3idthe2 * inv_rm2
        + dnh3idthe * costhe_inv_rm2sinthe + d2nh3idphi2 * inv_rm2sinthe2;

    double diffusion_nh4sh = d2nh4shdr2 + dnh4shdr * two_inv_rm + d2nh4shdthe2 * inv_rm2
        + dnh4shdthe * costhe_inv_rm2sinthe + d2nh4shdphi2 * inv_rm2sinthe2;


    // ===== RHS assembly =====
    // Non-dimensional scaling factors (all velocities are v/u_0, lengths r/L_atm):
    //   pressure:    × 1e5 / (r_mix · u_0²)        converts bar → dimensionless
    //   buoyancy:    × L_atm / u_0²                 converts g·(1−ρ/ρ_ref) → 1/Fr²
    //   Coriolis:    × L_atm / u_0                  converts [1/s] → 1/Ro
    //   centrifugal: × L_atm² / u_0²                converts [1/s²] → dimensionless
    double dpdr_term   = dpdr   * scale_p;
    double dpdthe_term = dpdthe * scale_p * inv_rm;
    double dpdphi_term = dpdphi * scale_p * inv_rmsinthe;

    // Einstein viscosity: mu_eff = mu_0*(1+2.5*phi)  =>  re_eff = re/(1+2.5*phi)
    const double phi =
          (h2o_cloud.x[i][j][k] + h2o_ice.x[i][j][k]) * r_mix / rho_cond_h2o
        + (ch4_cloud.x[i][j][k] + ch4_ice.x[i][j][k]) * r_mix / rho_cond_ch4
        + (h2s_cloud.x[i][j][k] + h2s_ice.x[i][j][k]) * r_mix / rho_cond_h2s
        + (nh3_cloud.x[i][j][k] + nh3_ice.x[i][j][k]) * r_mix / rho_cond_nh3
        + nh4sh.x[i][j][k]                            * r_mix / rho_cond_nh4sh;
    const double re_eff = re / (1.0 + 2.5 * phi);

    // ===== Turbulent (eddy) diffusion from the closure — stage three =====
    //
    // Until now the closure filled nue* and nothing read it: k* and dis* were transported, a
    // viscosity was computed from them, and no equation felt it. This is the step that closes that
    // loop. ATURAN_TURB_COUPLING is a MULTIPLIER, not a flag — default 0.0 leaves every equation
    // exactly as it was, and 1.0 applies the closure's own nue* in full, so the port can be walked
    // in rather than switched on.
    //
    // nue_t is added ALONGSIDE the existing molecular term rather than replacing it; re_eff
    // already carries Uranus's correction for condensate loading, which is a separate physical
    // effect and stays. Pr_t = 0.9 is the turbulent Prandtl number relating the scalar eddy
    // diffusivity to the momentum one.
    //
    // THE CAVEAT ATSAT FOUND APPLIES HERE WITH LESS FORCE BUT THE SAME SIGN: its nue* came out
    // ~77x SMALLER than the molecular background 1/re, so coupling changed almost nothing.
    // ATURAN's seeded nue* is 0.000165 against 1/re = 0.001 — 6.1x smaller. A coupling that does
    // almost nothing is worth knowing about before it is trusted.
    static const double turb_coupling = [](){
        const char* e = getenv("ATURAN_TURB_COUPLING"); return e ? atof(e) : 0.0; }();
    constexpr double Pr_t = 0.9;
    const double nue_t   = (turb_coupling != 0.0 && std::isfinite(nue.x[i][j][k]))
                         ? turb_coupling * std::max(0.0, nue.x[i][j][k]) : 0.0;
    const double nue_t_s = nue_t / Pr_t;      // scalar (heat / species) eddy diffusivity

    // ATURAN_THERMAL_MASSFLUX scales the diffusive-enthalpy sink below. DEFAULT 1.0, so the
    // model is unchanged unless it is set; 0.0 removes the term, which is what ATJUP's rhs_t does
    // permanently (it carries radiation_t + precip_t there and no thermalmassflux term at all).
    //
    // IT IS A KNOB BECAUSE THE TERM DOMINATES rhs_t ON THE ICE GIANTS BY ORDERS OF MAGNITUDE.
    // Measured at i=30, j=90, k=180, single-threaded, against the other three contributions:
    //
    //      model    pressure   transport   diffusion    chemical    thermalmassflux
    //      ATSAT     +1.3e-03   +1.6e-01    -4.3e-03    +1.6e-01      -0.163
    //      ATNEPT    +9.4e-07   +7.9e-03    -1.4e-01    +3.8e+01     -37.874
    //      ATURAN    -1.9e-07   +7.8e-04    +3.5e-02    -1.2e+02    +119.450
    //
    // On ATSAT it BALANCES the transport term, which is what a physical enthalpy flux should do.
    // On the ice giants it is three to six orders above every other term and simply sets rhs_t.
    //
    // The dominant part is the DIFFUSIVE-ENTHALPY product, not the reaction: measured at that cell
    // the reaction part is exactly zero (all w_* = 0) and the whole of thermalmassflux comes from
    //     (j_nh3*cp_nh3 + j_h2s*cp_h2s + j_nh4sh*cp_nh4sh) * (dtdr + |dtdthe|/rm + dtdphi/rmsinthe)
    // a flux multiplied by a TEMPERATURE GRADIENT — so a steeper gradient drives more cooling,
    // which steepens the gradient. On ATURAN that feedback nucleates a cold front near 1.5 bar and
    // advances it one grid cell per iteration until cells reach t_min_K() = 7.5 K, leaving a hole
    // at ~2.7 bar with a decoupled warm blob above it whose top is what the photosphere diagnostic
    // then reads. See the measurement in the commit that added this knob.
    static const double tmf_scale = [](){
        const char* e = getenv("ATURAN_THERMAL_MASSFLUX"); return e ? atof(e) : 1.0; }();

    // ===== Radiative heating source (ATURAN_RAD_COUPLING, default 0 = off, bit-identical) =====
    //
    // ATJUP's term, ported unchanged in form (RHS_Jup_Turb.cpp). Converts the diagnostic radiative
    // flux divergence Q_rad [W/m3], which Radiation.h fills, into a nondimensional temperature
    // tendency: physically dT/dt = Q_rad/(rho*cp), nondimensionalised by this model's energy-
    // equation scaling (radial length L_rad, velocity u_0, temperature t_ref):
    //
    //      radiation_t = rad_coupling * Q_rad * L_rad / (rho * cp_mix * u_0 * t_ref)
    //
    // rho is the LOCAL density, not r_mix: Q_rad comes from real radiative fluxes, so dividing by
    // the local density is the correct conversion, and it is what makes the thin, cold upper
    // atmosphere respond strongly instead of being held by the deep column's density. rho_mix
    // already carries p_stat*1e5/(R_mix*T) (ChemistryUran.h), so it is read rather than recomputed;
    // the ideal-gas fallback keeps the term alive if rho_mix has not been filled yet.
    //
    // L_rad is L_atm converted to METRES. L_atm is in km on all four models and this is the same
    // conversion computeHydrostaticPressure() needs — the place the split's port dropped it.
    //
    // WHY THIS TERM MATTERS HERE, measured on this model rather than assumed. Over 224 iterations
    // the equatorial column redistributes without heating: the deep loses 58.5 K (405.11 -> 346.63)
    // and the top gains 61.9 K (77.15 -> 139.07) while the COLUMN MEAN barely moves
    // (245.99 -> 243.87, -0.9%). So there is no heating excess to cancel; what is missing is
    // anything that ties the top of the column to the planet's energy budget. With Q_rad absent
    // from rhs_t the top has no thermal anchor and diffusion drags it to the column mean, which is
    // why the photosphere reads 136 K against a T_eff(in) of 59 K and emits 30x the budget.
    //
    // rad_coupling = 1.0 IS THE PHYSICALLY CORRECT VALUE — the expression is the exact
    // nondimensional form of dT/dt = Q/(rho*cp) under this scaling. Values >> 1 are NOT a scaling
    // correction but a deliberate ACCELERATION factor, and ATJUP's note records why one is
    // tempting: the radiative relaxation time vastly exceeds the step, so equilibration at
    // coupling 1.0 needs far more iterations than a 224-step run provides. On Uranus, colder and
    // with a longer radiative time constant than Jupiter, that separation is wider still. Expect
    // coupling 1.0 to move this model very little in 224 iterations; that is a statement about the
    // run length, not about the term.
    static const double rad_coupling = [](){
        const char* e = getenv("ATURAN_RAD_COUPLING"); return e ? atof(e) : 0.0; }();
    double radiation_t = 0.0;
    if(rad_coupling != 0.0){
        const double T_phys = t.x[i][j][k] * t_ref;            // [K]
        const double P_phys = p_stat.x[i][j][k] * 1.0e5;       // p_stat is in bar -> [Pa]
        const double rho_f  = rho_mix.x[i][j][k];
        const double rho    = (rho_f > 0.0 && std::isfinite(rho_f))
                            ? rho_f
                            : ((T_phys > 1.0) ? P_phys / (R_mix * T_phys) : 0.0);  // [kg/m3]
        const double L_rad  = L_atm * 1.0e3;                   // atmosphere thickness [m]
        if(rho > 0.0 && cp_mix > 0.0){
            radiation_t = rad_coupling * Q_rad.x[i][j][k] * L_rad
                        / (rho * cp_mix * u_0 * t_ref);
            // Explicit-scheme stability limiter, ATJUP's value and reasoning: in a very low-density
            // cell the 1/rho factor can make the tendency blow up and destabilise the pole. Guard
            // non-finite FIRST — the cap's >/< tests are both false for NaN and would let it
            // through — then cap, which preserves the sign and only bites on a runaway.
            constexpr double rad_t_max = 0.5;
            if(!std::isfinite(radiation_t)) radiation_t = 0.0;
            else if(radiation_t >  rad_t_max) radiation_t =  rad_t_max;
            else if(radiation_t < -rad_t_max) radiation_t = -rad_t_max;
        }
    }

    rhs_t.x[i][j][k] =
        + pressure_t
        - transport_t
        + diffusion_t / (re * pr) + diffusion_t * nue_t_s
        - tmf_scale * chemical_reaction * thermalmassflux.x[i][j][k]
        + radiation_t;

    // Item 1's instrument: the same five terms, recorded separately for one column at stage 0.
    // Signs are as they enter the sum above, so the five add to tbud_tot exactly. See the note on
    // the members in cUranusModel.h. Inert unless ATURAN_TBUDGET is set.
    if(tbud_stage == 0 && j == tbudget_column_j() && k == tbudget_column_k()
       && !tbud_tot.empty() && i >= 0 && i < (int)tbud_tot.size()){
        tbud_pres[i]  = pressure_t;
        tbud_trans[i] = -transport_t;
        tbud_diff[i]  = diffusion_t / (re * pr) + diffusion_t * nue_t_s;
        tbud_tmf[i]   = -tmf_scale * chemical_reaction * thermalmassflux.x[i][j][k];
        tbud_rad[i]   = radiation_t;
        tbud_tot[i]   = rhs_t.x[i][j][k];
    }

    // Sponge layer: quadratic Rayleigh damping over the top quarter of the domain.
    // frac = 0 at i_sponge_start, 1 at i=im-1 → damping rate = alpha_sponge * frac².
    const double frac_sp = std::max(0.0, (double)(i - (im - 1) * 3 / 4) / (double)((im - 1) / 4));
    const double sponge  = alpha_sponge * frac_sp * frac_sp;

    // ATURAN_BUOY_SCALE scales the buoyancy term in rhs_u below. DEFAULT 1.0 — the model is
    // unchanged unless it is set. Same name and semantics as ATJUP's knob.
    //
    // IT EXISTS BECAUSE THE TERM AS WRITTEN IS ~1e5 SMALLER THAN IT SHOULD BE. p_stat is in BAR;
    // the ideal-gas density p/(R*T) needs pascals, so the expression wants a factor 1e5. The very
    // same expression carries that 1e5 where it fills the DIAGNOSTIC BuoyancyForce array in
    // Forces(); only the momentum equation is missing it. Measured at the equator, iteration 1:
    //
    //      i    -dpdr        buoyancy(as written)   shortfall   buoyancy*1e5   vs -dpdr
    //     10   +7.77e-02        +2.23e-06            3.5e4x      +2.23e-01       2.87
    //     25   +4.51e-02        +1.06e-06            4.3e4x      +1.06e-01       2.35
    //     35   +2.55e-02        +4.28e-07            6.0e4x      +4.28e-02       1.68
    //
    // So buoyancy is not fighting the pressure gradient and losing — it is absent. Restored it
    // would be 1.7-2.9x the pressure gradient, a first-order force. That is why this model's radial
    // velocity is ~1e-4 where ATSAT's is ~1e-1, why there is no overturning circulation, and why
    // thermal diffusion flattens the column unopposed.
    //
    // WHY A KNOB AND NOT A FIX. The commented-out line directly below records that a correctly
    // scaled density-anomaly buoyancy "leads to ozillations in the upper region". The version left
    // in place is stable precisely BECAUSE it does nothing, so restoring the factor may reproduce
    // exactly what was fled from. Both siblings that behave replaced the raw rho*g rather than
    // rescaling it — ATSAT with an anomaly against buoy_ref_level, ATJUP with an exact hydrostatic
    // split — and ATJUP's nondimensionalisation carries 1.0e5*(L_atm*1.0e3)/(u_0*u_0), which is not
    // the same factor this expression is short by. Measure before believing any of it.
    static const double buoy_scale = [](){
        const char* e = getenv("ATURAN_BUOY_SCALE"); return e ? atof(e) : 1.0; }();

    // ---- Hydrostatic split (ATURAN_HYDRO_SPLIT, default 0 = off, bit-identical) ----
    //
    // With the split on the radial buoyancy is carried by p_hydro instead of appearing in rhs_u,
    // and what enters the momentum equation from it is the HORIZONTAL gradient of that pressure.
    // Those two changes are the whole of it: rhs_u loses its buoyancy term, rhs_v and rhs_w gain a
    // much smaller one, and the difference between those magnitudes is the point.
    //
    // ATURAN_BUOY_REF=1 without the split subtracts buoy_ref_level from the buoyancy in place —
    // ATSAT's treatment — so the anomaly and the split can be attributed separately.
    static const int hydro_split = [](){
        const char* e = getenv("ATURAN_HYDRO_SPLIT"); return e ? atoi(e) : 0; }();
    static const int buoy_ref_on = [](){
        const char* e = getenv("ATURAN_BUOY_REF"); return e ? atoi(e) : 0; }();

    // Buoyancy as it enters rhs_u. Unchanged by default; zero when the split carries it.
    double buoyancy_u = 0.0;
    if(hydro_split == 0 && t.x[i][j][k] > 0.0){
        const double b_raw = g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])
                           / (r_mix * R_mix * t.x[i][j][k] * t_ref);
        buoyancy_u = buoy_scale * buoyancy * (L_atm / (u_0 * u_0))
                   * (buoy_ref_on != 0 ? (b_raw - buoy_ref_level[i]) : b_raw);
    }
    double dphdthe_term = 0.0, dphdphi_term = 0.0;
    if(hydro_split != 0){
        dphdthe_term = (p_hydro.x[i][j+1][k] - p_hydro.x[i][j-1][k]) * inv_2dthe * inv_rm;
        dphdphi_term = (p_hydro.x[i][j][k+1] - p_hydro.x[i][j][k-1]) * inv_2dphi * inv_rmsinthe;
    }

    rhs_u.x[i][j][k] =
        - dpdr_term
//        + buoyancy * (L_atm / (u_0 * u_0)) * g * (1.0 - rho_mix.x[i][j][k] / r_mix) // leads to ozillations in the upper region
        + buoyancy_u                                          // see the hydrostatic-split note above
        - transport_u
        + diffusion_u / re_eff + diffusion_u * nue_t
        - Coriolis    * scale_Cor * Coriolis_rad
        // PLUS: the force points away from the axis. See the block where it is formed.
        + (cent_legacy ? -1.0 : +1.0) * centrifugal * scale_cen * centrifugal_rad
        - sponge * u_ijk;

    rhs_v.x[i][j][k] =
        - dpdthe_term
        - dphdthe_term
        - transport_v
        + diffusion_v / re_eff + diffusion_v * nue_t
        - Coriolis    * scale_Cor * Coriolis_the
        + (cent_legacy ? -1.0 : +1.0) * centrifugal * scale_cen * centrifugal_the;

    rhs_w.x[i][j][k] =
        - dpdphi_term
        - dphdphi_term
        - transport_w
        + diffusion_w / re_eff + diffusion_w * nue_t
        - Coriolis    * Coriolis_phi;

    rhs_ch4.x[i][j][k] =
        - transport_ch4
        + diffusion_ch4 / (sc_ch4 * re) + diffusion_ch4 * nue_t_s;

    rhs_ch4_cloud.x[i][j][k] =
        - transport_ch4_cloud
        + diffusion_ch4_cloud / (sc_ch4 * re) + diffusion_ch4_cloud * nue_t_s;

    rhs_ch4_ice.x[i][j][k] =
        - transport_ch4_ice
        + diffusion_ch4_ice / (sc_ch4 * re) + diffusion_ch4_ice * nue_t_s;

    rhs_h2o.x[i][j][k] =
        - transport_h2o
        + diffusion_h2o / (sc_h2o * re) + diffusion_h2o * nue_t_s;

    rhs_h2o_cloud.x[i][j][k] =
        - transport_h2o_cloud
        + diffusion_h2o_cloud / (sc_h2o * re) + diffusion_h2o_cloud * nue_t_s;

    rhs_h2o_ice.x[i][j][k] =
        - transport_h2o_ice
        + diffusion_h2o_ice / (sc_h2o * re) + diffusion_h2o_ice * nue_t_s;

    rhs_h2s.x[i][j][k] =
        - transport_h2s
        + diffusion_h2s / (sc_h2s * re) + diffusion_h2s * nue_t_s
        + chemical_reaction * massflux_h2s.x[i][j][k];

    rhs_h2s_cloud.x[i][j][k] =
        - transport_h2s_cloud
        + diffusion_h2s_cloud / (sc_h2s * re) + diffusion_h2s_cloud * nue_t_s;

    rhs_h2s_ice.x[i][j][k] =
        - transport_h2s_ice
        + diffusion_h2s_ice / (sc_h2s * re) + diffusion_h2s_ice * nue_t_s;

    rhs_nh3.x[i][j][k] =
        - transport_nh3
        + diffusion_nh3 / (sc_nh3 * re) + diffusion_nh3 * nue_t_s
        + chemical_reaction * massflux_nh3.x[i][j][k];

    rhs_nh3_cloud.x[i][j][k] =
        - transport_nh3_cloud
        + diffusion_nh3_cloud / (sc_nh3 * re) + diffusion_nh3_cloud * nue_t_s;

    rhs_nh3_ice.x[i][j][k] =
        - transport_nh3_ice
        + diffusion_nh3_ice / (sc_nh3 * re) + diffusion_nh3_ice * nue_t_s;

    // Stokes terminal velocity for NH4SH crystals falling in the -r direction.
    // v_stokes [m/s] = (2/9) * r_p² * (rho_crystal - rho_mix) * g / mue_mix
    // Divided by u_0 to get the non-dimensional sedimentation velocity;
    // positive sign because downward settling ≡ negative radial velocity,
    // giving +v_sed * dq/dr in the concentration equation.
    const double v_stokes_nh4sh =
        (2.0 / 9.0) * r_p_nh4sh * r_p_nh4sh
        * (rho_cond_nh4sh - rho_mix.x[i][j][k]) * g / mue_mix;

    rhs_nh4sh.x[i][j][k] =
        - transport_nh4sh
        + fluxlim_nh4sh.x[i][j][k]
        + diffusion_nh4sh / (sc_nh4sh * re) + diffusion_nh4sh * nue_t_s
        + chemical_reaction * massflux_nh4sh.x[i][j][k]
        + (v_stokes_nh4sh / u_0) * dnh4shdr;

    /*
     * k* and dis* become PROGNOSTIC here — stage two of the turbulence port.
     *
     * Transport + diffusion + the closure's own source terms, the standard k-omega transport
     * equations. sig_k and sig_w are the transport Prandtl numbers of the closure, not derived
     * from it; TurbulenceUran carries its own sig_w2 for the SST cross-diffusion.
     *
     * With the closure OFF this block is skipped and both tendencies are zeroed, so the integrator
     * leaves k* and dis* exactly where they were — which is what keeps ATURAN_TURB unset
     * byte-identical.
     */
    if(turb_active){
        constexpr double sig_k = 0.85, sig_w = 0.5;
        const double nue_here = nue.x[i][j][k];

        const double transport_tke = u.x[i][j][k] * dtkedr
                                   + v.x[i][j][k] * dtkedthe * geo.inv_rm
                                   + w.x[i][j][k] * dtkedphi * geo.inv_rmsinthe;
        const double transport_dis = u.x[i][j][k] * ddisdr
                                   + v.x[i][j][k] * ddisdthe * geo.inv_rm
                                   + w.x[i][j][k] * ddisdphi * geo.inv_rmsinthe;

        const double diffusion_tke = d2tkedr2 + dtkedr * 2.0 * geo.inv_rm
                                   + d2tkedthe2 * geo.inv_rm2
                                   + dtkedthe * geo.costhe_inv_rm2sinthe
                                   + d2tkedphi2 * geo.inv_rm2sinthe2;
        const double diffusion_dis = d2disdr2 + ddisdr * 2.0 * geo.inv_rm
                                   + d2disdthe2 * geo.inv_rm2
                                   + ddisdthe * geo.costhe_inv_rm2sinthe
                                   + d2disdphi2 * geo.inv_rm2sinthe2;

        rhs_tke.x[i][j][k] = - transport_tke
                           + diffusion_tke * (1.0/re + nue_here/sig_k)
                           + tke_source.x[i][j][k];
        rhs_dis.x[i][j][k] = - transport_dis
                           + diffusion_dis * (1.0/re + nue_here/sig_w)
                           + dis_source.x[i][j][k];
    } else {
        rhs_tke.x[i][j][k] = 0.0;
        rhs_dis.x[i][j][k] = 0.0;
    }

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr_term;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe_term;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi_term;
}
/*
*
*/

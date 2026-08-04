/*
 * The four boundary-condition FIELD LISTS, for the shared BoundaryConditions<Planet>.
 *
 * Transcribed VERBATIM from the arrays inside BC_Uran.h's bcRadius/bcTheta/bcPhi, which already
 * used this shape — a list walked once — so nothing here is a new decision. They are moved out to
 * the model because the shared header asks the planet which fields it has, and no two planets
 * answer alike: ATURAN's radius list carries h2s_cloud and h2s_ice, which ATSAT's does not, and
 * omits rho_mix, which ATNEPT's carries.
 *
 * The theta pass is two lists on purpose. Scalars and tracers are extrapolated to the poles; the
 * mass and diffusive fluxes are ZEROED there instead, because they are singular on the axis.
 *
 * THE PHI LIST IS SHORTER THAN THE RADIUS LIST, and deliberately so rather than by oversight:
 * BC_Uran.h omits j_nh4sh, jT_nh4sh and difflux_nh4sh from its longitudinal pass while including
 * them radially. That is transcribed as it stands. Whether it is intended is a question for the
 * model, not for a commit that is only moving the lists.
 */
#include "cUranusModel.h"

std::vector<Array*> cUranusModel::bc_fields_radius(){
    return {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s, &h2s_cloud, &h2s_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,   &j_nh3,   &j_nh4sh,
        &jT_h2s,  &jT_nh3,  &jT_nh4sh,
        &w_h2s,   &w_nh3,   &w_nh4sh,
        &massflux_h2s,  &massflux_nh3, &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,   &difflux_nh3,  &difflux_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
}

std::vector<Array*> cUranusModel::bc_fields_theta_extrap(){
    return {
        &t, &u,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s, &h2s_cloud, &h2s_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,   &j_nh3,   &j_nh4sh,
        &jT_h2s,  &jT_nh3,  &jT_nh4sh,
        &w_h2s,   &w_nh3,   &w_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
}

std::vector<Array*> cUranusModel::bc_fields_theta_zero(){
    return {
        // v and w FIRST, and they are the reason this list is not just the flux fields: BC_Uran
        // zeroes the two tangential velocities at the poles with four hand-written lines that sit
        // OUTSIDE its zero_at_poles array, so transcribing only that array would leave them
        // extrapolated instead. ATNEPT hit exactly this and its acceptance test caught it.
        &v, &w,
        &massflux_h2s, &massflux_nh3, &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,  &difflux_nh3, &difflux_nh4sh
    };
}

std::vector<Array*> cUranusModel::bc_fields_phi(){
    return {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s, &h2s_cloud, &h2s_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,   &j_nh3,
        &jT_h2s,  &jT_nh3,
        &w_h2s,   &w_nh3,  &w_nh4sh,
        &massflux_h2s,  &massflux_nh3, &massflux_nh4sh,
        &fluxlim_nh4sh,
        &difflux_h2s,   &difflux_nh3,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce,
        &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
}

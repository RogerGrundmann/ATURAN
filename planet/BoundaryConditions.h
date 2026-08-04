/*
 * SHARED PHYSICS — the three domain boundary conditions. MUST BE BYTE-IDENTICAL IN EVERY MODEL
 * THAT USES IT; `make check-shared` verifies that against planet/SHARED.md5.
 *
 * bcRadius (the deep boundary i=0 and the model top i=im-1), bcTheta (the two poles) and bcPhi
 * (the seam at k=0/km-1) were the same three routines in ATJUP and ATSAT: the same loops applying
 * one extrapolation formula to a list of transported fields, plus the same five hardening knobs.
 * Sharing them was impossible while ATSAT's copy was written longhand — one named block per field
 * per boundary, ~110 copies of the same two lines — and became possible when it was restructured
 * into field lists on 2026-08-02.
 *
 * ===== WHAT THE MODEL SUPPLIES, AND WHY EACH ONE IS A MODEL FACT =====
 *
 *   bc_fields_radius()          The lists themselves. 37 of 38 names are common, but the tails
 *   bc_fields_theta_extrap()    are genuinely per-planet: ATJUP carries six prognostic turbulence
 *   bc_fields_theta_zero()      fields in its main lists where ATSAT gives them their own
 *   bc_fields_phi()             treatment, ATSAT carries fluxlim_nh4sh where ATJUP does not, and
 *                               ATJUP's phi list omits j_nh4sh and jT_nh4sh. Those are the models
 *                               disagreeing about which fields exist and are transported, which
 *                               no shared file can decide for them.
 *
 *                               ONE OF THESE IS A REAL PHYSICS DIFFERENCE, not bookkeeping: at
 *                               the poles ATSAT EXTRAPOLATES massflux_* and difflux_* while ATJUP
 *                               ZEROES them. Both behaviours survive here untouched, because
 *                               deciding between them is a modelling question and this is a port.
 *
 *   bc_turb_fields()            ATSAT's turbulence treatment, which ATJUP does not have: tke, dis
 *   bc_turb_floors()            and nue extrapolated with the 2-point form and clamped to k*>=0,
 *   bc_turb_active()            dis*>=dis_min, gated on the closure being on. The cubic amplifies
 *                               an alternating error 7x per call and dis* sits in denominators
 *                               throughout the closure, which is why it must not be the cubic.
 *                               ATJUP returns an EMPTY list and keeps its turbulence fields in
 *                               the main lists, exactly as before. Promoting ATSAT's treatment to
 *                               ATJUP is the better physics and is deliberately NOT done here: it
 *                               would change ATJUP's results, and a sharing commit that changes
 *                               results cannot be checked.
 *
 *   Planet::bc_margin()         Which rows the loops cover. ATSAT works the interior
 *                               (j=1..jm-2, k=1..km-2) and leaves the corner rows to the other two
 *                               routines; ATJUP covers the full range. Returns 1 and 0.
 *
 *   Planet::bc_default_form()   THE reason this file can exist at all — see below.
 *
 *   Planet::bc_default_*()      Each knob's default. ATJUP ships rigid_lid and top_taper ON, on
 *                               Jupiter evidence; ATSAT ships everything OFF, because on Saturn
 *                               none of it is measured yet.
 *
 * ===== THE FORM ENCODING, WHICH IS WHERE THE TWO MODELS NEARLY COLLIDED =====
 *
 * Both models already had BC_POLE_COPY and BC_RADIUS_COPY, and THE SAME VALUE MEANT DIFFERENT
 * THINGS: ATSAT's 0 is the 3-point cubic, ATJUP's 0 is the 2-point Neumann. Unifying on either
 * would have silently changed the other model's default. So 0 does not name a formula here:
 *
 *     0  DEFAULT  — whatever Planet::bc_default_form() says (ATSAT cubic, ATJUP Neumann)
 *     1  COPY     — f[s] = f[a], first-order, cannot overshoot
 *     2  NEUMANN  — f[s] = (4/3)f[a] - (1/3)f[b]
 *     3  CUBIC    — f[s] = f[c] - 3f[b] + 3f[a]
 *
 * Every value either model's users have ever set keeps its exact previous meaning: ATSAT's 1 and
 * 2 were already copy and Neumann, ATJUP's 1 was already copy, and each model's 0 resolves to
 * the form it already used. Nothing had to be renumbered.
 */

#pragma once

#include "ATPhys.h"

#include <vector>
#include <algorithm>
#include <cstdlib>

class Array;


namespace BCForm {
    enum { DEFAULT = 0, COPY = 1, NEUMANN = 2, CUBIC = 3 };

    // `a`, `b`, `c` are the first three interior neighbours inward from the boundary slot.
    inline double apply(int form, int model_default,
                        double a, double b, double c, double c43, double c13){
        const int f = (form == DEFAULT) ? model_default : form;
        if(f == COPY)    return a;
        if(f == NEUMANN) return c43 * a - c13 * b;
        return c - 3.0 * b + 3.0 * a;      // CUBIC
    }
}


template<class Planet>
class BoundaryConditions {
public:
    explicit BoundaryConditions(Planet& model) : m(model) {}

    void bcRadius();
    void bcTheta();
    void bcPhi();

private:
    Planet& m;

    // Each knob is read once and cached. The static lives in a template member function, so
    // there is one per planet instantiation, which is what we want.
    static int knob(const char* name, int dflt){
        return ATPhys::env_int(Planet::planet_tag(), name, dflt);
    }
    static int rigid_lid(){
        static const int v = knob("BC_RIGID_LID", Planet::bc_default_rigid_lid()); return v; }
    static int lid_u(){
        static const int v = knob("BC_LID_U", rigid_lid()); return v; }
    static int top_taper(){
        static const int v = knob("BC_TOP_TAPER", Planet::bc_default_top_taper()); return v; }
    static int pole_copy(){
        static const int v = knob("BC_POLE_COPY", Planet::bc_default_pole_copy()); return v; }
    static int radius_copy(){
        static const int v = knob("BC_RADIUS_COPY", Planet::bc_default_radius_copy()); return v; }
    static int seam_damp(){
        static const int v = knob("BC_SEAM_DAMP", 0); return v; }
    static int t_lid_pin(){
        static const int v = knob("BC_T_LID_PIN", 0); return v; }
    static double t_lid_K(){
        static const double v = ATPhys::env_double(Planet::planet_tag(), "BC_T_LID_K", 0.0);
        return v; }
};


// ---------------------------------------------------------------------------
// Deep boundary i=0 and model top i=im-1.
// ---------------------------------------------------------------------------
template<class Planet>
void BoundaryConditions<Planet>::bcRadius()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;
    const int mg   = Planet::bc_margin();
    const int fdef = Planet::bc_default_form();
    const int form = radius_copy();

    std::vector<Array*> fields = m.bc_fields_radius();
    const int nf = (int)fields.size();

    std::vector<Array*> turb   = m.bc_turb_fields();
    std::vector<double> floors = m.bc_turb_floors();
    const int nt = (int)turb.size();
    const bool turb_bc = m.bc_turb_active();

    // --- (6) Lid temperature pin. <TAG>_BC_T_LID_PIN, with <TAG>_BC_T_LID_K for the
    // prescribed-temperature mode. The snapshot is taken on the FIRST call, before the
    // extrapolation below overwrites t.x[im-1]: the first bcRadius() runs after all
    // initialisation, so it captures the initial condition. Serial, outside any parallel region,
    // and only when the knob is on.
    const bool do_t_pin = t_lid_pin() != 0;
    if(do_t_pin && (int)m.t_top_init.size() != jm){
        const double t_K = t_lid_K();
        m.t_top_init.assign(jm, std::vector<double>(km, 0.0));
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                m.t_top_init[j][k] = (t_K > 0.0) ? (t_K / m.t_ref) : m.t.x[im-1][j][k];
    }
    const bool pin_t_top = do_t_pin && ((int)m.t_top_init.size() == jm);

    #pragma omp parallel for
    for(int j = mg; j < jm-mg; j++){
        for(int k = mg; k < km-mg; k++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[0][j][k] = BCForm::apply(form, fdef,
                    F.x[1][j][k], F.x[2][j][k], F.x[3][j][k], c43, c13);
                F.x[im-1][j][k] = BCForm::apply(form, fdef,
                    F.x[im-2][j][k], F.x[im-3][j][k], F.x[im-4][j][k], c43, c13);
            }

            // The turbulence fields take the 2-point form and a floor regardless of the knob
            // above — see the note at the top of this file. Empty list on ATJUP.
            if(turb_bc){
                for(int f = 0; f < nt; f++){
                    Array& F = *turb[f];
                    F.x[0][j][k] = std::max(floors[f],
                        c43 * F.x[1][j][k] - c13 * F.x[2][j][k]);
                    F.x[im-1][j][k] = std::max(floors[f],
                        c43 * F.x[im-2][j][k] - c13 * F.x[im-3][j][k]);
                }
            }
        }
    }

    // The three corrections below run in their own FULL-RANGE loops rather than inside the loop
    // above, because that loop may be interior-only (bc_margin) and a lid, a pin or a taper with
    // four untreated edges is none of those things. Their order relative to each other and to the
    // extrapolation is preserved; each touches a different quantity and reads no neighbour, so
    // running them per-array rather than per-cell changes nothing.

    // --- (1) Rigid radial walls on the wall-normal velocity u. ---
    // The modelled shell is closed: no mass crosses the model top, and none crosses the deep
    // boundary either. The extrapolation above is correct only for the TANGENTIAL v,w; applied to
    // u it admits a spurious mass flux through a closed boundary. It also bears on well-posedness,
    // the pressure solve being all-Neumann: with nothing pinning the wall-normal velocity the
    // column-mean vertical velocity is an undetermined, freely drifting constant.
    //
    // un is NOT written here: restoreVar(1.0) copies u -> un after all three BC routines run, so
    // the next Runge-Kutta step already starts from the wall value.
    if(lid_u() != 0){
      #pragma omp parallel for
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                m.u.x[0][j][k]    = 0.0;
                m.u.x[im-1][j][k] = 0.0;
            }
        }
    }

    // (6) Override the lid temperature extrapolation with the pinned value.
    if(pin_t_top){
      #pragma omp parallel for
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                m.t.x[im-1][j][k] = m.t_top_init[j][k];
    }

    // --- (2) Taper the HORIZONTAL velocities to a quiet grid ceiling. ---
    // Both models ramp v,w to zero between the tropopause and the model top in the initial
    // condition, and the extrapolation above copies the interior value straight back onto the lid
    // and undoes it. Two forms:
    //
    //   =1  multiply in place by 2/3, 1/3, 0 over the top three layers.
    //   =2  re-derive from i=im-4, the first layer the taper does not touch.
    //
    // =1 COMPOUNDS, contrary to the comment it was ported with. RK4 forms v = vn + dt*rhs and
    // restoreVar copies v back into vn, so the factor multiplies a value it already multiplied.
    // Measured on ATSAT over 30 iterations, max|w|: i=38 0.5717 -> 0.000056, i=39 0.4362 -> 0,
    // i=40 0.1505 -> 0, i.e. (2/3)^30. It is a hard three-layer zero reached geometrically rather
    // than the graded ramp it reads as. =2 holds the intended profile because i=im-4 is refreshed
    // by RK4 every iteration. i=im-1 is outside the RK4 range in both forms.
    const int taper = top_taper();
    if(taper != 0){
        const int iml = im - 1;
      #pragma omp parallel for
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                m.v.x[iml][j][k] = 0.0;
                m.w.x[iml][j][k] = 0.0;
                if(taper == 2){
                    m.v.x[iml-1][j][k] = (1.0/3.0) * m.v.x[iml-3][j][k];
                    m.w.x[iml-1][j][k] = (1.0/3.0) * m.w.x[iml-3][j][k];
                    m.v.x[iml-2][j][k] = (2.0/3.0) * m.v.x[iml-3][j][k];
                    m.w.x[iml-2][j][k] = (2.0/3.0) * m.w.x[iml-3][j][k];
                }else{
                    m.v.x[iml-1][j][k] *= (1.0/3.0);
                    m.w.x[iml-1][j][k] *= (1.0/3.0);
                    m.v.x[iml-2][j][k] *= (2.0/3.0);
                    m.w.x[iml-2][j][k] *= (2.0/3.0);
                }
            }
        }
    }
}


// ---------------------------------------------------------------------------
// The two poles, j=0 and j=jm-1.
// ---------------------------------------------------------------------------
template<class Planet>
void BoundaryConditions<Planet>::bcTheta()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;
    const int mg   = Planet::bc_margin();
    const int fdef = Planet::bc_default_form();

    // At the spherical singularity sin(theta) -> 0 any extrapolation amplifies grid noise: the
    // 2-point form by 4/3 per call, the cubic by ~7. A plain copy has amplification exactly 1.0
    // and is the axisymmetric-pole assumption to first order, which is what the pole physically
    // is. On ATSAT the cubic's overshoot was measured large enough to make the model's global
    // max|u| a polar boundary value while the interior maximum was unchanged to six digits.
    const int form = pole_copy();

    std::vector<Array*> fields = m.bc_fields_theta_extrap();
    const int nf = (int)fields.size();
    std::vector<Array*> zeros = m.bc_fields_theta_zero();
    const int nz = (int)zeros.size();

    std::vector<Array*> turb   = m.bc_turb_fields();
    std::vector<double> floors = m.bc_turb_floors();
    const int nt = (int)turb.size();
    const bool turb_bc = m.bc_turb_active();

    #pragma omp parallel for
    for(int k = mg; k < km-mg; k++){
        for(int i = mg; i < im-mg; i++){
            for(int f = 0; f < nz; f++){
                Array& F = *zeros[f];
                F.x[i][0][k]    = 0.0;
                F.x[i][jm-1][k] = 0.0;
            }

            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[i][0][k] = BCForm::apply(form, fdef,
                    F.x[i][1][k], F.x[i][2][k], F.x[i][3][k], c43, c13);
                F.x[i][jm-1][k] = BCForm::apply(form, fdef,
                    F.x[i][jm-2][k], F.x[i][jm-3][k], F.x[i][jm-4][k], c43, c13);
            }

            if(turb_bc){
                for(int f = 0; f < nt; f++){
                    Array& F = *turb[f];
                    F.x[i][0][k] = std::max(floors[f],
                        c43 * F.x[i][1][k] - c13 * F.x[i][2][k]);
                    F.x[i][jm-1][k] = std::max(floors[f],
                        c43 * F.x[i][jm-2][k] - c13 * F.x[i][jm-3][k]);
                }
            }
        }
    }
}


// ---------------------------------------------------------------------------
// The phi seam, k=0 and k=km-1.
// ---------------------------------------------------------------------------
template<class Planet>
void BoundaryConditions<Planet>::bcPhi()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    // k=0 and k=km-1 are the SAME meridian, so this is not a boundary condition in the sense of
    // the two above: each face is extrapolated from its own side and the two are averaged and set
    // equal, because a jump between them is a discontinuity in the middle of the domain. The
    // 2-point form is used here in both models and is deliberately NOT under the form knob —
    // there is no wall to be first-order accurate about.
    std::vector<Array*> fields = m.bc_fields_phi();
    const int nf = (int)fields.size();

    std::vector<Array*> turb   = m.bc_turb_fields();
    std::vector<double> floors = m.bc_turb_floors();
    const int nt = (int)turb.size();
    const bool turb_bc = m.bc_turb_active();

    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                const double lo = c43 * F.x[i][j][1]    - c13 * F.x[i][j][2];
                const double hi = c43 * F.x[i][j][km-2] - c13 * F.x[i][j][km-3];
                F.x[i][j][0] = F.x[i][j][km-1] = 0.5 * (lo + hi);
            }

            // The floor is applied to the AVERAGE rather than to each face, so the two stay
            // exactly equal.
            if(turb_bc){
                for(int f = 0; f < nt; f++){
                    Array& F = *turb[f];
                    const double lo = c43 * F.x[i][j][1]    - c13 * F.x[i][j][2];
                    const double hi = c43 * F.x[i][j][km-2] - c13 * F.x[i][j][km-3];
                    F.x[i][j][0] = F.x[i][j][km-1] = std::max(floors[f], 0.5 * (lo + hi));
                }
            }
        }
    }

    // --- (4) Shapiro damping across the phi seam. Default 0 passes in both models. ---
    // The reconstruction above pins k=0/km-1 to the mean of their own neighbours, which drops the
    // discrete d2/dphi2 self-damping at k=1 and k=km-2 from -2 to -1.5 — a 25 % loss of numerical
    // zonal diffusion exactly at the seam. Combined with the 1/sin^2(theta) metric that leaves an
    // under-damped zonal mode which runs away in ATOM. Neither ATJUP nor ATSAT shows such a mode,
    // which is why this is off: a filter that perturbs u,v,w every iteration would only confound
    // later experiments. Set to 2 if zonal energy is ever seen accumulating at k = 0/1/km-2.
    //
    // Only u,v,w are damped; restoreVar copies them to un,vn,wn after all BCs. A solid neighbour
    // acts as no-flux, contributing the cell's own value.
    const int npass = seam_damp();
    if(npass > 0){
        constexpr double seam_coeff = 0.25;
        const int km2 = km - 2;
        const int km3 = km - 3;

        #pragma omp parallel for
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                auto is_fluid = [&](int kk){ return !m.is_solid(i, j, kk); };

                auto smooth_seam = [&](Array& F){
                    const double f_km2 = F.x[i][j][km2];
                    const double f_0   = F.x[i][j][0];
                    const double f_1   = F.x[i][j][1];

                    const bool fl_km2 = is_fluid(km2);
                    const bool fl_0   = is_fluid(0);
                    const bool fl_1   = is_fluid(1);

                    const double w_km2 = is_fluid(km3) ? F.x[i][j][km3] : f_km2;
                    const double e_km2 = fl_0 ? f_0   : f_km2;
                    const double w_0   = fl_km2 ? f_km2 : f_0;
                    const double e_0   = fl_1 ? f_1   : f_0;
                    const double w_1   = fl_0 ? f_0   : f_1;
                    const double e_1   = is_fluid(2) ? F.x[i][j][2] : f_1;

                    if(fl_km2)
                        F.x[i][j][km2] = f_km2 + seam_coeff*(w_km2 - 2.0*f_km2 + e_km2);
                    if(fl_0)
                        F.x[i][j][0] = F.x[i][j][km-1] =
                            f_0 + seam_coeff*(w_0 - 2.0*f_0 + e_0);
                    if(fl_1)
                        F.x[i][j][1] = f_1 + seam_coeff*(w_1 - 2.0*f_1 + e_1);
                };

                for(int p = 0; p < npass; p++){
                    smooth_seam(m.u);
                    smooth_seam(m.v);
                    smooth_seam(m.w);
                }
            }
        }
    }
}

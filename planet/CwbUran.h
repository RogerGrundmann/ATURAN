/*
 * Uranus Atmosphere Circulation Model (ATURAN)
 * ATURAN_CWB_DIAG=<n> -- the column budget of H2O, NH3, CH4 and H2S (vapour + cloud + ice each),
 * PRINT-ONLY, default off; printed every n iterations. ATSAT's ATSAT_CWB_DIAG (planet/CwbSat.h,
 * itself ATJUP's ATJUP_CWB_DIAG), ported 2026-10-12 and widened to this model's four condensables.
 *
 * WHY. After ATURAN_PDYN_INIT=0 the H2O, NH3 and H2S columns no longer fall between the saturation
 * adjustment's calls but rise slightly (+0.09 / +0.02 / +0.02 % in 224 iterations, satchk/giants/UP0), and
 * nothing in the model said in which term.
 *
 * This instrument marks each column after every stage of the time loop that can write the twelve
 * fields and charges the difference to that stage (the stages ATURAN_TB_DIAG marks). The
 * Runge-Kutta step is split further:
 *   rk_reset   it integrates from the n-level copies, which restoreVar refreshed at the END of
 *              the previous iteration -- every direct write made since then and not copied into
 *              them (ATURAN_RK_SCALAR_SYNC) is discarded here;
 *   rk terms   the RK4-weighted column integrals of each right-hand-side term of the three
 *              equations of a gas, accumulated per thread inside RHSUran: radial / meridional /
 *              zonal transport (advective form, -u.grad q), diffusion, the chemistry source
 *              (NH3 and H2S), the precipitation source (ATURAN_PRECIP_COUPLING);
 *   rk_rest    what the step changed beyond reset + terms: the floor of the final update, the
 *              limiters and round-off. An identity check when none of them binds.
 * q*div(u) is accumulated beside them and NOT applied: the transport is written in advective
 * form, and -u.grad q = -div(q u) + q div(u), so on a closed column the advective form creates or
 * destroys exactly the integral of q*div(u). Printed so the two can be compared.
 * With ATURAN_SPECIES_DIVU set the term IS applied (RHS_Uran_Turb.cpp) and appears as its own row,
 * "flux-form term -s q div(u)"; the not-applied row stays as the reference.
 *
 * Unit: g/m2, the mean over the sphere weighted by sin(colatitude), layer thickness in metres on
 * densities in kg/m3 -- the same sum the saturation adjustment's instrument prints. The top level
 * has no thickness and is not counted. The plane i = 0 is not integrated (bcRadius sets it) and
 * its part of the change is printed apart from the integrated cells.
 *
 * The state lives in a function-local static so that the three translation units share one
 * instance without a model member (C++11 has no inline variables).
 */

#pragma once

namespace UranCwb {

    enum { G_H2O = 0, G_NH3, G_CH4, G_H2S, NGAS };
    enum { T_RAD = 0, T_THE, T_PHI, T_DIFF, T_CHEM, T_PRECIP, T_DIVU, T_QDIV, NTERM };
    const int MAXTHR = 256;

    struct State {
        bool   on;                               // ATURAN_CWB_DIAG
        double stage_wgt;                        // RK4 weight of the stage being evaluated (1, 2, 2, 1)
        double acc[MAXTHR][NGAS][8];             // [thread][gas][term], weighted sums of the tendencies
    };

    inline State& st(){
        static State s = State();
        return s;
    }
}

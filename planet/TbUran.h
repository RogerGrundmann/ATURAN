/*
 * Uranus Atmosphere Circulation Model (ATURAN)
 * ATURAN_TB_DIAG=<n> -- the temperature budget by level, PRINT-ONLY, default off; printed every n
 * iterations.
 *
 * Ported from ATSAT (2026-10-11), where it showed that the thermal mass flux term of the
 * chemistry was the whole of the temperature drift and that the time step discarded the
 * saturation adjustment. The precipitation row stays empty: rhs_t has no such term here.
 *
 * This instrument marks the level-mean temperature (weighted by sin(colatitude), in K) after
 * every stage of the time loop that can write t and charges the difference to that stage, as
 * ATJUP_CWB_DIAG does for the water column. The Runge-Kutta step is split further:
 *   rk_reset   the step integrates from tn, which restoreVar refreshed at the END of the previous
 *              iteration, so every direct write to t made since then is discarded here;
 *   rk terms   the RK4-weighted level means of each term of rhs_t, accumulated per thread in
 *              RHSUran: pressure work, the three transport terms, diffusion, the thermal mass
 *              flux of the chemistry, radiation, precipitation;
 *   rk_rest    what the step changed beyond reset + terms: the limiter, and the rows and planes
 *              the step does not integrate (they keep t while the mean is over all of them).
 *
 * The state lives in a function-local static so that the three translation units share one
 * instance without a model member (C++11 has no inline variables).
 */

#pragma once

namespace UranTb {

    enum { T_PRESS = 0, T_TR_RAD, T_TR_THE, T_TR_PHI, T_DIFF, T_TMF, T_RAD, T_PRECIP, NTERM };
    const int MAXTHR = 256;
    const int MAXLEV = 96;

    struct State {
        bool   on;                               // ATURAN_TB_DIAG
        double stage_wgt;                        // RK4 weight of the stage being evaluated (1, 2, 2, 1)
        double acc[MAXTHR][NTERM][MAXLEV];       // [thread][term][level], weighted sums of the tendencies
    };

    inline State& st(){
        static State s = State();
        return s;
    }
}

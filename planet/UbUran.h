/*
 * Uranus Atmosphere Circulation Model (ATURAN)
 * ATURAN_UB_DIAG=<n> -- the budget of the radial wind u by level, PRINT-ONLY, default off; printed
 * every n iterations.
 *
 * Written 2026-10-13 for the continuity residuum under the rigid lid: the residuum the run ends
 * on is du/dr in the last cell, it grows because u grows at the equator, and the question was
 * whether that growth is made in the wall cells or by the pressure balance wearing off.
 *
 * Built like ATURAN_TB_DIAG (TbNept.h). The mean of u over a level (weighted by sin(colatitude),
 * in m/s) is marked after every stage of the time loop that can write u and the difference is
 * charged to that stage; the Runge-Kutta step is split into the RK4-weighted means of each term
 * of rhs_u, accumulated per thread in RHSUran. Two regions: the whole level, and the equatorial
 * band the residuum sits in.
 *
 * The state lives in a function-local static so that the three translation units share one
 * instance without a model member (C++11 has no inline variables).
 */

#pragma once

namespace UranUb {

    // the terms of rhs_u, with the sign they enter it
    enum { T_DPDR = 0, T_METRIC, T_ADV, T_DIFF, T_COR, T_CEN, T_BUOY, T_SPONGE, NTERM };
    enum { R_ALL = 0, R_EQ, NREG };
    const int MAXTHR = 256;
    const int MAXLEV = 96;

    struct State {
        bool   on;                                   // ATURAN_UB_DIAG
        double stage_wgt;                            // RK4 weight of the stage being evaluated (1, 2, 2, 1)
        int    j_lo, j_hi;                           // the equatorial band, rows j_lo..j_hi
        double acc[MAXTHR][NREG][NTERM][MAXLEV];     // [thread][region][term][level], weighted sums of the tendencies
    };

    inline State& st(){
        static State s = State();
        return s;
    }
}

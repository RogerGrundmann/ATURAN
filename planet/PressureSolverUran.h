/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone pressure-solver class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * Header-only: all method bodies are inline.
*/

#pragma once

#include "cUranusModel.h"
#include "Utils.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class PressureSolverUran {
public:
    explicit PressureSolverUran(cUranusModel& model)
        : m(model)
    {}

    void run()
    {
        using namespace std;
        cout << endl << "      ATURAN: computePressure" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double c43 = m.c43, c13 = m.c13;

        // precompute sin(the) table — only depends on j, avoids redundant sin() calls
        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            // Same floor the integrator and the shared solver use; see cUranusModel.h.
            if (sinthe_table[j] < m.sinthe_min()) sinthe_table[j] = m.sinthe_min();
        }

        // Boundary conditions for aux — r-direction
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {
                m.aux_u.x[0][j][k]      = c43 * m.aux_u.x[1][j][k]      - c13 * m.aux_u.x[2][j][k];
                m.aux_u.x[m.im-1][j][k] = c43 * m.aux_u.x[m.im-2][j][k] - c13 * m.aux_u.x[m.im-3][j][k];
                m.aux_v.x[0][j][k]      = c43 * m.aux_v.x[1][j][k]      - c13 * m.aux_v.x[2][j][k];
                m.aux_v.x[m.im-1][j][k] = c43 * m.aux_v.x[m.im-2][j][k] - c13 * m.aux_v.x[m.im-3][j][k];
                m.aux_w.x[0][j][k]      = c43 * m.aux_w.x[1][j][k]      - c13 * m.aux_w.x[2][j][k];
                m.aux_w.x[m.im-1][j][k] = c43 * m.aux_w.x[m.im-2][j][k] - c13 * m.aux_w.x[m.im-3][j][k];
            }
        }


        // Grid-spacing reciprocals — constant for the entire grid
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const double inv_dr2   = 1.0 / (m.dr   * m.dr);
        const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
        const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);

        // Main Poisson solve — inline geometry per (i,j), rhs subtracted in divergence
        //
        // ===== THIS LOOP IS SERIAL, AND THAT IS THE FIX FOR THE THREADING DEFECT =====
        //
        // It is Gauss-Seidel IN PLACE: p_dyn[i][j][k] is written from p_dyn[i±1][j][k] and
        // p_dyn[i][j±1][k]. Those neighbours differ from (i,j) in exactly the two indices the
        // pragma that used to sit here handed to different threads:
        //
        //     #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        //
        // so cell (i,j,k) was read by the thread owning (i+1,j) or (i,j+1) while its owner was
        // writing it. schedule(dynamic) made it worse than a thread-count dependence: which
        // thread got which chunk varied with timing, so the SAME binary at the SAME thread count
        // gave different answers run to run.
        //
        // This is the defect PressureSolver.h records fixing in the SHARED solver, in the same
        // words — "writing p_dyn in place while reading p_dyn[i±1][j±1]" — cured there by
        // red-black colouring. This solver is the DEFAULT (ATURAN_PRESS_SOLVER=0) and never got
        // it, so the fix landed in the path almost nothing runs.
        //
        // MEASURED, nm=4, single run of this build. With the pragma: two 16-thread runs differ,
        // and 1 vs 16 threads differ in 7 of 7 output files. Serial: 1 thread and 16 threads
        // agree BIT-IDENTICALLY, two 16-thread runs agree bit-identically, and all of them agree
        // bit-identically with the 1-thread answer this loop gave BEFORE the change — one thread
        // ran collapse(2) in lexicographic order already, so serialising takes nothing back.
        // EVERY single-threaded measurement in README.md still stands unchanged.
        //
        // Serial rather than red-black because it costs almost nothing: computePressure is 0.003 s
        // of a 5.3 s step at 16 threads and 0.03 s when serial, ~0.5 % of a step. Red-black would
        // change the answer and strand those measurements to recover 30 milliseconds.
        // ATURAN_PRESS_SOLVER=1 still selects the shared red-black solver for anyone who wants
        // this loop parallel.
        //
        // The four boundary loops in this file KEEP their pragmas: each writes one i-, j- or
        // k-plane and reads only planes it does not write, so none of them has this problem.
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                const double rm             = m.rad.z[i];
                const double rm2            = rm * rm;
                const double sinthe         = sinthe_table[j];
                const double inv_rm         = 1.0 / rm;
                const double inv_rm2        = 1.0 / rm2;
                const double inv_rmsinthe   = 1.0 / (rm * sinthe);
                const double inv_rm2sinthe2 = inv_rm2 / (sinthe * sinthe);

                const double denom     = 2.0 * inv_dr2
                                       + 2.0 * inv_rm2 * inv_dthe2
                                       + 2.0 * inv_rm2sinthe2 * inv_dphi2;
                const double inv_denom = 1.0 / denom;
                const double num1      = inv_dr2;
                const double num2      = inv_rm2 * inv_dthe2;
                const double num3      = inv_rm2sinthe2 * inv_dphi2;

                for (int k = 1; k < m.km-1; k++) {
                    const double du_dr   = ((m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k])) * inv_2dr;
                    const double dv_dthe = ((m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k])) * inv_2dthe * inv_rm;
                    const double dw_dphi = ((m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1])) * inv_2dphi * inv_rmsinthe;

                    m.p_dyn.x[i][j][k] =
                        ((m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                       + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                       + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3
                       - (du_dr + dv_dthe + dw_dphi)) * inv_denom;
                }
            }
        }

        // p_dyn boundary conditions — r-direction
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[0][j][k]      = c43 * m.p_dyn.x[1][j][k]      - c13 * m.p_dyn.x[2][j][k];
                m.p_dyn.x[m.im-1][j][k] = c43 * m.p_dyn.x[m.im-2][j][k] - c13 * m.p_dyn.x[m.im-3][j][k];
            }
        }

        // p_dyn boundary conditions — theta-direction
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                m.p_dyn.x[i][0][k]      = c43 * m.p_dyn.x[i][1][k]      - c13 * m.p_dyn.x[i][2][k];
                m.p_dyn.x[i][m.jm-1][k] = c43 * m.p_dyn.x[i][m.jm-2][k] - c13 * m.p_dyn.x[i][m.jm-3][k];
            }
        }

        // p_dyn boundary conditions — phi-direction (periodic average)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                m.p_dyn.x[i][j][0]      = c43 * m.p_dyn.x[i][j][1]      - c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = c43 * m.p_dyn.x[i][j][m.km-2] - c13 * m.p_dyn.x[i][j][m.km-3];
                m.p_dyn.x[i][j][0]      = m.p_dyn.x[i][j][m.km-1]
                    = (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1]) / 2.0;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for computePressure\n", elapsed.count() * 1e-9);
        cout << "      ATURAN: computePressure ended" << endl;
    }

private:
    cUranusModel& m;
};

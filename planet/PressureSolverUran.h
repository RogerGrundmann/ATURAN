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
            if (sinthe_table[j] < 0.4) sinthe_table[j] = 0.4;
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
        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
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

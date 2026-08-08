/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * This file holds prepareProjectionBoundaries() and nothing else. It used to open with a third
 * pressure solver, cUranusModel::computePressure(), a 220-line in-place Gauss-Seidel that had no
 * caller anywhere in the repository — the run loop calls PressureSolverUran or the shared
 * PressureSolver, never this. It was deleted rather than left dormant: its Poisson loop carried
 * the SAME in-place stencil that made the default solver non-reproducible above one thread, and
 * it hoisted denom, num1..num3 and daux_udr/daux_vdthe/daux_wdphi to function scope, so anyone
 * who parallelised it to "speed it up" would have written both that race and a second one on the
 * shared scalars. ATSAT removed its inherited second solver for the same reason and says so in
 * PressureSolverSat.h. Deleting it changes no output: it never ran.
*/
#include "cUranusModel.h"

using namespace std;

/*
 * Radial-wall conditions for the pressure projection, mirrored from ATSAT's
 * cSaturnModel::prepareProjectionBoundaries. Required by the SHARED PressureSolver.h,
 * which asks the model to state its own wall treatment rather than assuming one.
 *
 * Reached only when ATURAN_PRESS_SOLVER=1 selects the shared solver; ATURAN's own
 * PressureSolverUran remains the default and does not call this.
 */
void cUranusModel::prepareProjectionBoundaries(bool rigid_lid){
        #pragma omp parallel for
        for(int j = 1; j < jm-1; j++){          // r-direction
            for(int k = 1; k < km-1; k++){
                if(rigid_lid){
                    aux_u.x[0][j][k]      = 0.0;
                    aux_u.x[im-1][j][k] = 0.0;
                } else {
                aux_u.x[0][j][k] = aux_u.x[3][j][k]
                    - 3.0 * aux_u.x[2][j][k] + 3.0 * aux_u.x[1][j][k];
                aux_u.x[im-1][j][k] = aux_u.x[im-4][j][k]
                    - 3.0 * aux_u.x[im-3][j][k] + 3.0 * aux_u.x[im-2][j][k];
                }

                aux_v.x[0][j][k] = aux_v.x[3][j][k]
                    - 3.0 * aux_v.x[2][j][k] + 3.0 * aux_v.x[1][j][k];
                aux_v.x[im-1][j][k] = aux_v.x[im-4][j][k]
                    - 3.0 * aux_v.x[im-3][j][k] + 3.0 * aux_v.x[im-2][j][k];

                aux_w.x[0][j][k] = aux_w.x[3][j][k]
                    - 3.0 * aux_w.x[2][j][k] + 3.0 * aux_w.x[1][j][k];
                aux_w.x[im-1][j][k] = aux_w.x[im-4][j][k]
                    - 3.0 * aux_w.x[im-3][j][k] + 3.0 * aux_w.x[im-2][j][k];

                rhs_u.x[0][j][k] = rhs_u.x[3][j][k]
                    - 3.0 * rhs_u.x[2][j][k] + 3.0 * rhs_u.x[1][j][k];
                rhs_u.x[im-1][j][k] = rhs_u.x[im-4][j][k]
                    - 3.0 * rhs_u.x[im-3][j][k] + 3.0 * rhs_u.x[im-2][j][k];

                rhs_v.x[0][j][k] = rhs_v.x[3][j][k]
                    - 3.0 * rhs_v.x[2][j][k] + 3.0 * rhs_v.x[1][j][k];
                rhs_v.x[im-1][j][k] = rhs_v.x[im-4][j][k]
                    - 3.0 * rhs_v.x[im-3][j][k] + 3.0 * rhs_v.x[im-2][j][k];

                rhs_w.x[0][j][k] = rhs_w.x[3][j][k]
                    - 3.0 * rhs_w.x[2][j][k] + 3.0 * rhs_w.x[1][j][k];
                rhs_w.x[im-1][j][k] = rhs_w.x[im-4][j][k]
                    - 3.0 * rhs_w.x[im-3][j][k] + 3.0 * rhs_w.x[im-2][j][k];
            }
        }

        #pragma omp parallel for
        for(int k = 1; k < km-1; k++){          // theta-direction
            for(int i = 1; i < im-1; i++){
                aux_u.x[i][0][k] = aux_u.x[i][3][k]
                    - 3.0 * aux_u.x[i][2][k] + 3.0 * aux_u.x[i][1][k];
                aux_v.x[i][0][k] = 0.0;
                aux_w.x[i][0][k] = 0.0;

                aux_u.x[i][jm-1][k] = aux_u.x[i][jm-4][k]
                    - 3.0 * aux_u.x[i][jm-3][k] + 3.0 * aux_u.x[i][jm-2][k];
                aux_v.x[i][jm-1][k] = 0.0;
                aux_w.x[i][jm-1][k] = 0.0;

                rhs_u.x[i][0][k] = rhs_u.x[i][3][k]
                    - 3.0 * rhs_u.x[i][2][k] + 3.0 * rhs_u.x[i][1][k];
                rhs_v.x[i][0][k] = 0.0;
                rhs_w.x[i][0][k] = 0.0;

                rhs_u.x[i][jm-1][k] = rhs_u.x[i][jm-4][k]
                    - 3.0 * rhs_u.x[i][jm-3][k] + 3.0 * rhs_u.x[i][jm-2][k];
                rhs_v.x[i][jm-1][k] = 0.0;
                rhs_w.x[i][jm-1][k] = 0.0;
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < im; i++){            // phi-direction
            for(int j = 0; j < jm; j++){
                aux_u.x[i][j][0] = c43 * aux_u.x[i][j][1] - c13 * aux_u.x[i][j][2];
                aux_u.x[i][j][km-1] = c43 * aux_u.x[i][j][km-2] - c13 * aux_u.x[i][j][km-3];
                aux_u.x[i][j][0] = aux_u.x[i][j][km-1] =
                    (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.0;

                aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
                aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
                aux_v.x[i][j][0] = aux_v.x[i][j][km-1] =
                    (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.0;

                aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
                aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
                aux_w.x[i][j][0] = aux_w.x[i][j][km-1] =
                    (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.0;

                rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
                rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
                rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] =
                    (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.0;

                rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
                rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
                rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] =
                    (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.0;

                rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
                rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
                rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] =
                    (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.0;
            }
        }

}


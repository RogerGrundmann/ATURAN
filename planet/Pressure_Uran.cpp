/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/
#include "cUranusModel.h"

using namespace std;

void cUranusModel::computePressure(){
    cout << endl << "      ATURAN: computePressure" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int j = 1; j < jm-1; j++){  // r-direction
        for(int k = 1; k < km-1; k++){

            aux_u.x[0][j][k] = aux_u.x[3][j][k] 
                - 3.0 * aux_u.x[2][j][k] + 3.0 * aux_u.x[1][j][k];  // extrapolation
            aux_u.x[im-1][j][k] = aux_u.x[im-4][j][k] 
                - 3.0 * aux_u.x[im-3][j][k] + 3.0 * aux_u.x[im-2][j][k];

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
    for(int k = 1; k < km-1; k++){  // theta-direction
        for(int i = 1; i < im-1; i++){

            aux_u.x[i][0][k] = aux_u.x[i][3][k] 
                - 3.0 * aux_u.x[i][2][k] + 3.0 * aux_u.x[i][1][k];  // extrapolation
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
    for(int i = 0; i < im; i++){  // phi-direction
        for(int j = 0; j < jm; j++){

            aux_u.x[i][j][0] = c43 * aux_u.x[i][j][1] - c13 * aux_u.x[i][j][2];  // von Neumann
            aux_u.x[i][j][km-1] = c43 * aux_u.x[i][j][km-2] - c13 * aux_u.x[i][j][km-3];
            aux_u.x[i][j][0] = aux_u.x[i][j][km-1] = (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.0;

            aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
            aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
            aux_v.x[i][j][0] = aux_v.x[i][j][km-1] = (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.0;

            aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
            aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
            aux_w.x[i][j][0] = aux_w.x[i][j][km-1] = (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.0;


            rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
            rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
            rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] = (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.0;

            rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
            rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
            rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] = (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.0;

            rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
            rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
            rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] = (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.0;
        }
    }


    double rm = 0.0;
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0.0;
    double rmsinthe = 0.0;
    double denom = 0.0;
    double num1 = 0.0;
    double num2 = 0.0;
    double num3 = 0.0;
    double daux_udr = 0.0;
    double daux_vdthe = 0.0;
    double daux_wdphi = 0.0;

    for(int i = 1; i < im-1; i++){
        rm = metricRadius(rad.z[i]);   // must match the integrator's metric, or the projection
                                       // solves a different geometry than the momentum equation

        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            if(sinthe == 0.0) sinthe = 1.0e-5;
            rmsinthe = rm * sinthe;
            denom = 2.0/dr2 + 2.0/(rm * rm * dthe2) 
                + 2.0/(rmsinthe * rmsinthe * dphi2);
            num1 = 1.0/dr2;
            num2 = 1.0/(rm * rm * dthe2);
            num3 = 1.0/(rmsinthe * rmsinthe * dphi2);

            for(int k = 1; k < km-1; k++){
                daux_udr = 
                     ((aux_u.x[i+1][j][k] - aux_u.x[i-1][j][k]) 
                    - (rhs_u.x[i+1][j][k] - rhs_u.x[i-1][j][k]))
                    /(2.0 * dr);
                daux_vdthe = 
                     ((aux_v.x[i][j+1][k] - aux_v.x[i][j-1][k]) 
                    - (rhs_v.x[i][j+1][k] - rhs_v.x[i][j-1][k]))
                    /(2.0 * dthe * rm);
                daux_wdphi = 
                     ((aux_w.x[i][j][k+1] - aux_w.x[i][j][k-1]) 
                    - (rhs_w.x[i][j][k+1] - rhs_w.x[i][j][k-1]))
                    /(2.0 * dphi * rmsinthe);


                p_dyn.x[i][j][k] = 
                     ((p_dyn.x[i+1][j][k] + p_dyn.x[i-1][j][k]) * num1 
                    + (p_dyn.x[i][j+1][k] + p_dyn.x[i][j-1][k]) * num2 
                    + (p_dyn.x[i][j][k+1] + p_dyn.x[i][j][k-1]) * num3 
                    - (daux_udr + daux_vdthe + daux_wdphi))/denom;

            }  // end k
        }  // end j
    }  // end i





    #pragma omp parallel for
    for(int k = 1; k < km-1; k++){  // r-direction
        for(int j = 1; j < jm-1; j++){

            p_dyn.x[0][j][k] = p_dyn.x[3][j][k] 
                - 3.0 * p_dyn.x[2][j][k] 
                + 3.0 * p_dyn.x[1][j][k];  // extrapolation
            p_dyn.x[im-1][j][k] = p_dyn.x[im-4][j][k] 
                - 3.0 * p_dyn.x[im-3][j][k] 
                + 3.0 * p_dyn.x[im-2][j][k];  // extrapolation

        }
    }

    #pragma omp parallel for
    for(int k = 1; k < km-1; k++){  // theta-direction
        for(int i = 1; i < im-1; i++){

            p_dyn.x[i][0][k] = p_dyn.x[i][3][k] 
                - 3.0 * p_dyn.x[i][2][k] + 3.0 * p_dyn.x[i][1][k];  // extrapolation
            p_dyn.x[i][jm-1][k] = p_dyn.x[i][jm-4][k] 
                - 3.0 * p_dyn.x[i][jm-3][k] + 3.0 * p_dyn.x[i][jm-2][k];  // extrapolation

        }
    }

    #pragma omp parallel for
    for(int i = 0; i < im; i++){  // phi-direction
        for(int j = 0; j < jm; j++){

            p_dyn.x[i][j][0] = c43 * p_dyn.x[i][j][1] - c13 * p_dyn.x[i][j][2];  // von Neumann
            p_dyn.x[i][j][km-1] = c43 * p_dyn.x[i][j][km-2] - c13 * p_dyn.x[i][j][km-3];
            p_dyn.x[i][j][0] = p_dyn.x[i][j][km-1] = (p_dyn.x[i][j][0] + p_dyn.x[i][j][km-1])/2.0;

        }
    }

/*
    #pragma omp parallel for
     for(int k = 0; k < km; k++){  // supresses whiggles in corners
        for(int j = 0; j <= 7; j++){
            for(int i = 0; i < im; i++){
                p_dyn.x[i][j][k] = 0.0;
            }
        }
        for(int j = jm-8; j < jm; j++){
            for(int i = 0; i < im; i++){
                p_dyn.x[i][j][k] = 0.0;
            }
        }
    }
*/
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for computePressure\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: computePressure ended" << endl;
    return;
}
/*
*
*/


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


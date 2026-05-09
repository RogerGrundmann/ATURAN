/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone velocity-initialiser class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * Header-only: all method bodies are inline.
 * Mirrors cUranusModel::UranusCellStructure() and its helpers.
*/

#pragma once

#include "cUranusModel.h"
#include "Utils.h"

#include <iostream>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class VelocityInitializerUran {
public:
    explicit VelocityInitializerUran(cUranusModel& model)
        : m(model)
    {}

    // ========================================================================
    // Main entry point — mirrors cUranusModel::UranusCellStructure()
    // ========================================================================
    void compute()
    {
        using namespace std;
        const int im = m.im, jm = m.jm, km = m.km;

        cout << endl << "      AGCM: init_velocities" << endl;

        // u — seed at pole, Hadley-cell boundary (45°), and equator;
        // form_diagonals produces the sign change exactly at 45° N/S
        init_u(m.u,   0);   // N. pole (descending)
        init_u(m.u,  45);   // 45°N: Hadley cell boundary (zero crossing)
        init_u(m.u,  90);   // equator (ascending)
        init_u(m.u, 135);   // 45°S: Hadley cell boundary (zero crossing)
        init_u(m.u, 180);   // S. pole (descending)

        // v — single Hadley cell per hemisphere: zero at pole and equator,
        //     maximum at mid-hemisphere (45°) with poleward flow aloft and
        //     equatorward return flow at the surface
        init_v_or_w(m.v,  90, 0.0,  0.0);   // equator: stagnation (ascending)

        // northern Hadley cell
//        init_v_or_w(m.v,   0, 0.0,  0.0);   // N. pole: stagnation (descending)
//        init_v_or_w(m.v,  45, -3.0, 2.5);   // mid-cell: poleward aloft, equatorward at surface
        init_v_or_w(m.v,   0, 0.0,  0.0);   // N. pole: stagnation (descending)
        init_v_or_w(m.v,  45, -0.3, 0.25);   // mid-cell: poleward aloft, equatorward at surface

        // southern Hadley cell (same magnitudes; sign flip applied below)
//        init_v_or_w(m.v, 180, 0.0,  0.0);   // S. pole: stagnation
//        init_v_or_w(m.v, 135, -3.0, 2.5);   // mid-cell
        init_v_or_w(m.v, 180, 0.0,  0.0);   // S. pole: stagnation
        init_v_or_w(m.v, 135, -0.3, 0.25);   // mid-cell

        // w — unchanged
//        init_v_or_w(m.w,  90, -0.65, -0.65);
        init_v_or_w(m.w,  90, -65.0, -65.0);

        init_v_or_w(m.w,   0, 0.0, 0.0);
        init_v_or_w(m.w, 180, 0.0, 0.0);

        init_v_or_w(m.w,  30, 60.0, 0.0);
        init_v_or_w(m.w, 150, 60.0, 0.0);

        init_v_or_w(m.w,  45, 185.0, 0.0);
        init_v_or_w(m.w, 135, 185.0, 0.0);

        init_v_or_w(m.w,  60, 220, 0.0);
        init_v_or_w(m.w, 120, 220.0, 0.0);


        // form diagonals — northern hemisphere
        form_diagonals(m.u,  0,  45);        // u: pole → 45°N (descending → zero)
        form_diagonals(m.u, 45,  90);        // u: 45°N → equator (zero → ascending)

        form_diagonals(m.v,  0,  45);        // v: pole → mid-cell
        form_diagonals(m.v, 45,  90);        // v: mid-cell → equator

        form_diagonals(m.w,  0,  30);
        form_diagonals(m.w, 30,  45);
        form_diagonals(m.w, 45,  60);
        form_diagonals(m.w, 60,  90);

        // form diagonals — southern hemisphere
        form_diagonals(m.u,  90, 135);       // u: equator → 45°S (ascending → zero)
        form_diagonals(m.u, 135, 180);       // u: 45°S → S. pole (zero → descending)

        form_diagonals(m.v,  90, 135);       // v: equator → mid-cell
        form_diagonals(m.v, 135, 180);       // v: mid-cell → S. pole

        form_diagonals(m.w,  90, 120);
        form_diagonals(m.w, 120, 135);
        form_diagonals(m.w, 135, 150);
        form_diagonals(m.w, 150, 180);

        // change direction for southern hemisphere v
        #pragma omp parallel for
        for(int i = 0; i < im; i++){
            for(int j = 91; j < jm; j++){
                for(int k = 0; k < km; k++){
                    m.v.x[i][j][k] = -m.v.x[i][j][k];
                }
            }
        }

        // non-dimensionalisation by u_0
        #pragma omp parallel for
        for(int i = 0; i < im; i++){
            for(int k = 0; k < km; k++){
                for(int j = 0; j < jm; j++){
                    m.u.x[i][j][k] = m.u.x[i][j][k]/m.u_0;
                    m.v.x[i][j][k] = m.v.x[i][j][k]/m.u_0;
                    m.w.x[i][j][k] = m.w.x[i][j][k]/m.u_0;
                }
            }
        }
        cout << "      AGCM: init_velocities ended" << endl;
        return;
    }

private:
    cUranusModel& m;

    void form_diagonals(Array& a, int start, int end)
    {
        const int im = m.im, km = m.km;
        #pragma omp parallel for
        for(int k = 0; k < km; k++){
            for(int j = start; j < end; j++){
                for(int i = 0; i < im; i++){
                    a.x[i][j][k] = (a.x[i][end][k] - a.x[i][start][k]) *
                        (j - start)/(double)(end - start) + a.x[i][start][k];
                }
            }
        }
        return;
    }

    void init_u(Array& u, int j)
    {
//        float ua_00 = 0.2894,
//              ua_90 = 0.12;   // ~mass-conserving: ua_00 * sin(45°)/(1-sin(45°))
        float ua_00 = 0.03,
              ua_90 = 0.01;   // ~mass-conserving: ua_00 * sin(45°)/(1-sin(45°))
        int tropopause_layer = m.get_tropopause_layer(j);
        float tropopause_height = m.get_layer_height(tropopause_layer);
        const int km = m.km;
        for(int k = 0; k < km; k++){
            for(int i = 0; i < tropopause_layer; i++){
                float layer_height = m.get_layer_height(i),
                    half_tropopause_height = tropopause_height/3.0;
                float ratio;
                if(layer_height < half_tropopause_height){
                    ratio = layer_height/(half_tropopause_height/2.2);
                    switch(j){
                        case 90:  u.x[i][90][k]  =  ua_00 * ratio; break;
                        case 45:  u.x[i][45][k]  =  0.0;           break;
                        case 135: u.x[i][135][k] =  0.0;           break;
                        case 0:   u.x[i][0][k]   = -ua_90 * ratio; break;
                        case 180: u.x[i][180][k] = -ua_90 * ratio; break;
                    }
                }else{
                    ratio = (tropopause_height-layer_height)/half_tropopause_height;
                    switch(j){
                        case 90:  u.x[i][90][k]  =  ua_00 * ratio; break;
                        case 45:  u.x[i][45][k]  =  0.0;           break;
                        case 135: u.x[i][135][k] =  0.0;           break;
                        case 0:   u.x[i][0][k]   = -ua_90 * ratio; break;
                        case 180: u.x[i][180][k] = -ua_90 * ratio; break;
                    }
                }
            }
        }
        return;
    }

    void init_v_or_w(Array& v_or_w, int j,
        double coeff_trop, double coeff_sl)
    {
        const int km = m.km;
        int tropopause_layer = m.get_tropopause_layer(j);
        double tropopause_height = m.get_layer_height(tropopause_layer);
        for(int k = 0; k < km; k++){
            if(is_ocean_surface(m.SeaMount, 20, j, k)){
                coeff_sl = v_or_w.x[0][j][k];
            }
            for(int i = 0; i < tropopause_layer; i++){
                v_or_w.x[i][j][k] = (coeff_trop - coeff_sl)
                    * m.get_layer_height(i)/tropopause_height + coeff_sl;
            }
        }
        init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);
        return;
    }

    void init_v_or_w_above_tropopause(Array& v_or_w, int j,
        double coeff)
    {
        const int im = m.im, km = m.km;
        int tropopause_layer = m.get_tropopause_layer(j);
        if(tropopause_layer >= im-1) return;
        double tropopause_height = m.get_layer_height(tropopause_layer);

        #pragma omp parallel for
        for(int k = 0; k < km; k++){
            for(int i = tropopause_layer; i < im; i++){
                v_or_w.x[i][j][k] = coeff
                    * (m.get_layer_height(im-1) - m.get_layer_height(i))
                    /(m.get_layer_height(im-1) - tropopause_height);
            }
        }
        return;
    }
};

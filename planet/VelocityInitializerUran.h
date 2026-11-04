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
#include <vector>
#include <cmath>
#include <algorithm>
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

        // ATURAN_V_INIT_SCALE=<factor> (2026-10-14, ported from ATNEPT; default 0 with
        // ATURAN_GEO_BALANCE, 1.0 = the cell of every run before): multiplies the prescribed
        // meridional overturning cell. Nothing balances it on the geostrophic start; on ATNEPT its
        // diffusion was the only term left in dv/dt at iteration 1.
        static const double v_init_scale = [](){
            const char* e = getenv("ATURAN_V_INIT_SCALE"); return e ? atof(e) : ((cUranusModel::geo_balance() && !cell_init()) ? 0.0 : 1.0); }();  // ATURAN_CELL_INIT=1 keeps the cell (see cell_init())
        if(v_init_scale != 1.0)
            printf("      ATURAN: ATURAN_V_INIT_SCALE = %.4g (prescribed meridional cell scaled)\n", v_init_scale);

        // non-dimensionalisation by u_0
        #pragma omp parallel for
        for(int i = 0; i < im; i++){
            for(int k = 0; k < km; k++){
                for(int j = 0; j < jm; j++){
                    m.u.x[i][j][k] = m.u.x[i][j][k]/m.u_0;
                    m.v.x[i][j][k] = (v_init_scale == 1.0) ? m.v.x[i][j][k]/m.u_0
                                                           : v_init_scale * m.v.x[i][j][k]/m.u_0;
                    m.w.x[i][j][k] = m.w.x[i][j][k]/m.u_0;
                }
            }
        }
        cout << "      AGCM: init_velocities ended" << endl;
        return;
    }

    // ATURAN_CELL_INIT (2026-11-04, ported from ATSAT_CELL_INIT, InitVelocity_Sat.cpp initCellStreamfunction;
    // user: "Neptune and Uranus need the cell structures initialized at the beginning too like Saturn").
    // DEFAULT 1 since 2026-11-04 (=0 restores u = 0 and the V_INIT_SCALE 0 start). With =1: (a) the prescribed meridional cell v is KEPT under ATURAN_GEO_BALANCE
    // (ATURAN_V_INIT_SCALE then defaults to 1 instead of 0; an explicit ATURAN_V_INIT_SCALE still wins), and
    // (b) after the initial damp_wiggles u is replaced by the radial branch that makes the cells divergence-free:
    // for an axisymmetric incompressible flow  d(rm^2 sin u)/dr + d(rm sin v)/dtheta = 0, so with a
    // streamfunction  rm^2 sin u = d psi/d theta,  rm sin v = - d psi/dr.  psi is integrated from v upward
    // (trapezoid on rad.z, psi = 0 at i = 0); a column whose v does not integrate to zero is closed by a
    // column-uniform return wind dv = psi(top) / (rm sin (r_top - r_0)), which also makes u = 0 at the lid.
    // u central in theta, pole rows copy their neighbours; built on k = 0, copied to every longitude; w untouched.
    // The TRUE sine is used, as on ATSAT; this model floors the metric sine in the RHS (ATURAN_SINTHE_FLOOR),
    // so poleward of the floor the cells are divergence-free in the true geometry only. Printed: max |u|,
    // max |dv| and the interior continuity residual (true sine, k = 1) before / after.
    static int cell_init(){
        static const int v = [](){ const char* e = getenv("ATURAN_CELL_INIT"); return e ? atoi(e) : 1; }();   // DEFAULT 1 since 2026-11-04 (user: "flip the ice giant defaults"; =0 restores)
        return v;
    }

    static double continuityResidualMax(cUranusModel& m, int& i_m, int& j_m){
        const int im = m.im, jm = m.jm;
        double mx = 0.0; i_m = j_m = 0;
        const int k = 1;
        for(int i = 1; i < im-1; i++) for(int j = 1; j < jm-1; j++){
            const double sinthe = sin(m.the.z[j]);
            double costhe = cos(m.the.z[j]);
            if(cUranusModel::costhe_abs() && j > 90) costhe = - costhe;
            const double rm = m.metricRadius(m.rad.z[i]);
            const double res = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k]) / (2. * m.dr) + 2. * m.u.x[i][j][k] / rm
                             + (m.v.x[i][j+1][k] - m.v.x[i][j-1][k]) / (2. * m.dthe) / rm + costhe / (rm * sinthe) * m.v.x[i][j][k];
            if(fabs(res) > mx){ mx = fabs(res); i_m = i; j_m = j; }
        }
        return mx;
    }

    void cellStreamfunction()
    {
        const int im = m.im, jm = m.jm, km = m.km;
        int i0_, j0_, i1_, j1_;
        const double res0 = continuityResidualMax(m, i0_, j0_);
        std::vector<double> psi((size_t)im * jm, 0.0), corr(jm, 0.0);
        const double r_span = m.rad.z[im-1] - m.rad.z[0];
        for(int j = 1; j < jm-1; j++){
            const double sj = sin(m.the.z[j]);
            double s_ = 0.0;
            for(int i = 1; i < im; i++){
                const double f0 = m.metricRadius(m.rad.z[i-1]) * sj * m.v.x[i-1][j][0], f1 = m.metricRadius(m.rad.z[i]) * sj * m.v.x[i][j][0];
                s_ -= 0.5 * (f0 + f1) * (m.rad.z[i] - m.rad.z[i-1]);
                psi[(size_t)i * jm + j] = s_;
            }
            corr[j] = psi[(size_t)(im-1) * jm + j];
            for(int i = 0; i < im; i++) psi[(size_t)i * jm + j] -= corr[j] * (m.rad.z[i] - m.rad.z[0]) / r_span;
        }
        double u_mx = 0.0, dv_mx = 0.0;
        std::vector<double> u_c((size_t)im * jm, 0.0), v_c((size_t)im * jm, 0.0);
        for(int i = 0; i < im; i++){
            const double rm = m.metricRadius(m.rad.z[i]);
            for(int j = 1; j < jm-1; j++){
                const double sj = sin(m.the.z[j]);
                u_c[(size_t)i * jm + j] = (psi[(size_t)i * jm + j+1] - psi[(size_t)i * jm + j-1]) / (2. * m.dthe) / (rm * rm * sj);
                const double dv = corr[j] / (rm * sj * r_span);
                v_c[(size_t)i * jm + j] = m.v.x[i][j][0] + dv;
                u_mx = std::max(u_mx, fabs(u_c[(size_t)i * jm + j])); dv_mx = std::max(dv_mx, fabs(dv));
            }
            u_c[(size_t)i * jm + 0] = u_c[(size_t)i * jm + 1];           u_c[(size_t)i * jm + jm-1] = u_c[(size_t)i * jm + jm-2];
            v_c[(size_t)i * jm + 0] = m.v.x[i][0][0];                    v_c[(size_t)i * jm + jm-1] = m.v.x[i][jm-1][0];
        }
        #pragma omp parallel for collapse(2)
        for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
            m.u.x[i][j][k] = u_c[(size_t)i * jm + j];
            m.v.x[i][j][k] = v_c[(size_t)i * jm + j];
        }
        const double res1 = continuityResidualMax(m, i1_, j1_);
        printf("      ATURAN: ATURAN_CELL_INIT: u from the streamfunction of v, max |u| %.4f m/s, column return wind max |dv| %.4f m/s;"
               " continuity residual (interior, k = 1, true sine) %.6f at (%d, %d) -> %.6f at (%d, %d)\n",
               u_mx * m.u_0, dv_mx * m.u_0, res0, i0_, j0_, res1, i1_, j1_);
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
        // ATURAN_U_INIT_SCALE (2026-10-12; DEFAULT 0 on the user's word, 1 restores): the amplitude
        // of the prescribed radial branch of the cells, as ATJUP_U_INIT_SCALE and
        // ATSAT_U_INIT_SCALE. On those two planets the branch peaked at 40 m/s and was the
        // model's velocity divergence. Here it is small (the amplitudes above are in m/s), and it
        // is switched off for the same reason and to treat the four models alike: it has no
        // matching horizontal branch. MEASUREMENT: see the README row of the knob.
        static const float u_init_scale = [](){ const char* e = getenv("ATURAN_U_INIT_SCALE"); return e ? (float)atof(e) : 0.0f; }();
        ua_00 *= u_init_scale; ua_90 *= u_init_scale;
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

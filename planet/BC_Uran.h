/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone boundary-condition class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * This is a header-only file: it contains the BC_Uran class, all inline method
 * bodies, and the inline cUranusModel delegation wrappers that forward each
 * cUranusModel::BC_xxx() call to the corresponding BC_Uran method.
 * BC_Uran.cpp is a minimal stub that provides a translation unit.
*/

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

class cUranusModel;
class Array;

using namespace std;

class BC_Uran {
public:

    explicit BC_Uran(cUranusModel& model) : m(model) {}

    void bcRadius();
    void bcTheta();
    void bcPhi();
    void initTropopauseLayers();

private:
    cUranusModel& m;
};


// -----------------------------------------------------------------------
// Inline implementation  (header-only, like ATOM's BC_Atm.h)
// -----------------------------------------------------------------------
#include "cUranusModel.h"
#include "Utils.h"

using namespace AtomUtils;


inline void BC_Uran::bcRadius()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    // All scalar/tracer fields: Neumann extrapolation at both radial boundaries.
    Array* fields[] = {
        &m.t, &m.u, &m.v, &m.w,
//        &m.t,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s, &m.h2s_cloud, &m.h2s_ice,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,  &m.j_nh4sh,
        &m.jT_h2s, &m.jT_nh3, &m.jT_nh4sh,
        &m.w_h2s,  &m.w_nh3,  &m.w_nh4sh,
        &m.massflux_h2s,  &m.massflux_nh3,  &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,   &m.difflux_nh3,   &m.difflux_nh4sh,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    // 2-point Neumann extrapolation: f[s] = (4/3)f[a] - (1/3)f[b].
    // The 3-point cubic (3f[a]-3f[b]+f[c]) amplifies alternating errors by 7x
    // per call and blows up near the SeaMount contour (same reason bcSolidGround
    // was switched to the 2-point formula; bcRadius had the same latent bug).
    #pragma omp parallel for schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[0][j][k]    = c43*F.x[1][j][k]    - c13*F.x[2][j][k];
                F.x[im-1][j][k] = c43*F.x[im-2][j][k] - c13*F.x[im-3][j][k];
            }
/*
            // Velocities: Neumann at inner surface, zero at outer boundary (top of atmosphere).
            m.u.x[0][j][k] = c43*m.u.x[1][j][k] - c13*m.u.x[2][j][k];
            m.v.x[0][j][k] = c43*m.v.x[1][j][k] - c13*m.v.x[2][j][k];
            m.w.x[0][j][k] = c43*m.w.x[1][j][k] - c13*m.w.x[2][j][k];
            m.u.x[im-1][j][k] = 0.0;
            m.v.x[im-1][j][k] = 0.0;
            m.w.x[im-1][j][k] = 0.0;
*/
//            m.t.x[im-1][j][k] = m.t_ref;
        }
    }
}


inline void BC_Uran::bcTheta()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    // All fields except v and w receive 2-point Neumann extrapolation at poles.
    Array* extrap_fields[] = {
        &m.t, &m.u,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s, &m.h2s_cloud, &m.h2s_ice,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,  &m.j_nh4sh,
        &m.jT_h2s, &m.jT_nh3, &m.jT_nh4sh,
        &m.w_h2s,  &m.w_nh3,  &m.w_nh4sh,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible
    };
    const int nf = (int)(sizeof(extrap_fields) / sizeof(extrap_fields[0]));

    // Flux fields that are singular-prone near poles: zero at pole boundary.
    Array* zero_at_poles[] = {
        &m.massflux_h2s, &m.massflux_nh3, &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,  &m.difflux_nh3,  &m.difflux_nh4sh,
    };
    const int nz = (int)(sizeof(zero_at_poles) / sizeof(zero_at_poles[0]));

    #pragma omp parallel for schedule(static)
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            m.v.x[i][0][k]    = 0.0;
            m.v.x[i][jm-1][k] = 0.0;
            m.w.x[i][0][k]    = 0.0;
            m.w.x[i][jm-1][k] = 0.0;

            for(int f = 0; f < nf; f++){
                Array& F = *extrap_fields[f];
                F.x[i][0][k]    = c43*F.x[i][1][k]    - c13*F.x[i][2][k];
                F.x[i][jm-1][k] = c43*F.x[i][jm-2][k] - c13*F.x[i][jm-3][k];
            }

            for(int f = 0; f < nz; f++){
                Array& F = *zero_at_poles[f];
                F.x[i][0][k]    = 0.0;
                F.x[i][jm-1][k] = 0.0;
            }
        }
    }
}


inline void BC_Uran::bcPhi()
{
    const int im = m.im, jm = m.jm, km = m.km;
    const double c43 = m.c43, c13 = m.c13;

    Array* fields[] = {
        &m.t, &m.u, &m.v, &m.w,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s, &m.h2s_cloud, &m.h2s_ice,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.nh4sh,
        &m.j_h2s,  &m.j_nh3,
        &m.jT_h2s, &m.jT_nh3,
        &m.w_h2s,  &m.w_nh3, &m.w_nh4sh,
        &m.massflux_h2s,  &m.massflux_nh3,  &m.massflux_nh4sh,
        &m.fluxlim_nh4sh,
        &m.difflux_h2s,   &m.difflux_nh3,
        &m.thermalmassflux,
        &m.CoriolisForce, &m.CentrifugalForce,
        &m.BuoyancyForce, &m.PresGradForce,
        &m.Q_Latent, &m.Q_Sensible
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                double lo = c43*F.x[i][j][1]    - c13*F.x[i][j][2];
                double hi = c43*F.x[i][j][km-2] - c13*F.x[i][j][km-3];
                F.x[i][j][0] = F.x[i][j][km-1] = 0.5*(lo + hi);
            }
        }
    }
}


inline void BC_Uran::initTropopauseLayers()
{
    const int jm = m.jm;

    m.tropopause_layers = std::vector<double>(jm, m.tropopause_pole);
    cout << endl << "      ATURAN: init_tropopause_layers" << endl;

    const int i_max  = m.im - 1;
    const int j_max  = jm - 1;
    const int j_half = j_max / 2;
    const double coeff_pole = 285.0;

    for(int j = j_half; j >= 0; j--){
        double x = coeff_pole * (1.0 - (double)(j_half - j) / (double)j_half);
        m.tropopause_layers[j] = AtomUtils::Agnesi(m.tropopause_equator, x);
        m.tropopause_layers[j] = std::round(m.tropopause_layers[j]
            / m.L_atm * (double)i_max);
        m.tropopause_layers[j] = m.tropopause_equator / m.L_atm * (double)i_max;
    }

    for(int j = j_max; j > j_half; j--)
        m.tropopause_layers[j] = m.tropopause_layers[j_max - j];

    cout << "      ATURAN: init_tropopause_layers ended" << endl;
}

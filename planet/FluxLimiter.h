/*
 * SHARED PHYSICS — TVD flux limiter for NH4SH advection. MUST BE BYTE-IDENTICAL IN EVERY MODEL
 * THAT USES IT; `make check-shared` verifies that against planet/SHARED.md5.
 *
 * Total Variation Diminishing (Superbee by default). Computes
 *
 *     fluxlim_nh4sh = transport_centered - transport_TVD
 *
 * which, added to rhs_nh4sh in the RHS, replaces the centred-difference advection with the
 * Superbee-limited upwind scheme and stops the spurious oscillations at the sharp NH4SH
 * cloud-formation boundary. To switch limiter, replace TVD::superbee_phi with TVD::van_leer_phi
 * at the three call sites below.
 *
 * ===== WHY THIS ONE IS SHARED WHEN THE REST OF THE CHEMISTRY IS NOT =====
 *
 * ChemistrySat and ChemistryJup hold four routines each. Three of them cannot be shared yet:
 * ChemMassRate is blocked on four unsettled disagreements (rate constants differing by 10x,
 * formation windows differing 3x, ATSAT's division by dt, and the diffusive sign), and
 * ThermalProperties and DiffMassFlux overlap only 31 % and 25 % line-for-line because the two
 * models compute genuinely different things — ATJUP's DiffMassFlux forms D*laplacian(c) directly
 * and carries Lewis-number groups, ATSAT's takes the divergence of a flux and has none.
 *
 * This routine is the exception: 97 % identical line for line, and the entire difference was
 *
 *     ATSAT   rm = m.metricRadius(m.rad.z[i])
 *     ATJUP   rm = m.rad.z[i]
 *
 * plus whitespace. metricRadius() is already the established hook for exactly this — PressureSolver
 * asks the model the same question — so the shared version asks it too.
 *
 * WHY THAT IS BIT-IDENTICAL IN BOTH, and it is not because either accessor does nothing. BOTH
 * models have a live metric radius; they apply it in different places, which is the whole reason
 * the hook exists:
 *
 *   ATJUP  ATJUP_METRIC_RADIUS is ON BY DEFAULT (unset means Jupiter's radius, cJupiterModel.cpp)
 *          and is applied by shifting rad.z ITSELF at initialisation, so rad.z[i] already carries
 *          it. metricRadius() is therefore the identity BY DESIGN, not by neglect, and
 *          m.metricRadius(m.rad.z[i]) is exactly the m.rad.z[i] the hand-written copy used.
 *   ATSAT  cannot do that — its rad.z is also the stretched coordinate the layer heights sit on,
 *          so shifting it would divide every radial derivative by ~116 — and instead shifts the
 *          metric FACTORS in the accessor. ATSAT_METRIC_RADIUS defaults to 0, i.e. off.
 *
 * So the two answer the same question from opposite ends and the shared call is right for both.
 * See cJupiterModel.h at metricRadius() and cSaturnModel.h at the same, which say this from each
 * model's side.
 *
 * The limiter is pure numerics — no rate constant, no temperature window, no species chemistry —
 * which is why it could be shared while the calibration questions around it stay open.
 */

#pragma once

#include "ATPhys.h"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class Array;


namespace TVD {
    inline double superbee_phi(double r){
        return std::max(0.0, std::max(std::min(2.0*r, 1.0), std::min(r, 2.0)));
    }
    inline double van_leer_phi(double r){
        return (r + std::abs(r)) / (1.0 + std::abs(r));
    }
}


template<class Planet>
class FluxLimiter {
public:
    explicit FluxLimiter(Planet& model) : m(model) {}

    void nh4sh();

private:
    Planet& m;
};


template<class Planet>
void FluxLimiter<Planet>::nh4sh()
{
    using namespace std;
    cout << endl << "      " << Planet::planet_tag() << ": FluxLimiterNH4SH" << endl;

    auto begin = chrono::high_resolution_clock::now();

    const int    im = m.im, jm = m.jm, km = m.km;
    const double dr   = m.dr;
    const double dthe = m.dthe;
    const double dphi = m.dphi;
    constexpr double eps = 1.0e-12;

    #pragma omp parallel for collapse(3) schedule(static)
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            for(int i = 1; i < im-1; i++){

                const double q    = m.nh4sh.x[i][j][k];
                const double q_rp = m.nh4sh.x[i+1][j][k];
                const double q_rm = m.nh4sh.x[i-1][j][k];
                const double q_tp = m.nh4sh.x[i][j+1][k];
                const double q_tm = m.nh4sh.x[i][j-1][k];
                const double q_pp = m.nh4sh.x[i][j][k+1];
                const double q_pm = m.nh4sh.x[i][j][k-1];

                const double u = m.u.x[i][j][k];
                const double v = m.v.x[i][j][k];
                const double w = m.w.x[i][j][k];

                const double rm           = m.metricRadius(m.rad.z[i]);
                const double sinthe       = max(ATPhys::polar_divisor_floor<Planet>(),
                                                abs(sin(m.the.z[j])));
                const double inv_rm       = 1.0 / rm;
                const double inv_rmsinthe = 1.0 / (rm * sinthe);

                double corr = 0.0;

                // antidiff = |u| * (q_{i+1} - 2q_i + q_{i-1}) / (2*dr)
                // correction = (1 - phi(r)) * antidiff

                // ---- r-direction ----
                {
                    const double df    = q_rp - q;
                    const double db    = q    - q_rm;
                    const double denom = df + (df >= 0.0 ? eps : -eps);
                    const double r = (u >= 0.0)
                        ? db / denom
                        : ((i+2 < im ? m.nh4sh.x[i+2][j][k] : q_rp) - q_rp) / denom;
                    corr += (1.0 - TVD::superbee_phi(r)) * abs(u) * (df - db) / (2.0 * dr);
                }

                // ---- theta-direction ----
                // Pole-symmetric handling: at j = 1 with v >= 0, q_tm = q[i][0][k]
                // is the Neumann-extrapolated boundary value (c43*q - c13*q_tp),
                // giving db = df/3 -> r = 1/3 (limiter partially active).
                // The mirror case at j = jm-2 with v < 0 wants q[i][jm][k], which
                // is off-grid; mirror the same Neumann extrapolation here so the
                // limiter behaves symmetrically across the equator instead of
                // collapsing to r = 0 (full antidiffusion) only at the south pole.
                {
                    const double df    = q_tp - q;
                    const double db    = q    - q_tm;
                    const double denom = df + (df >= 0.0 ? eps : -eps);
                    const double q_far = (j+2 < jm)
                        ? m.nh4sh.x[i][j+2][k]
                        : (m.c43 * q_tp - m.c13 * q);   // Neumann extrap of q[jm]
                    const double r = (v >= 0.0)
                        ? db / denom
                        : (q_far - q_tp) / denom;
                    corr += (1.0 - TVD::superbee_phi(r)) * abs(v) * inv_rm * (df - db) / (2.0 * dthe);
                }

                // ---- phi-direction ----
                {
                    const double df    = q_pp - q;
                    const double db    = q    - q_pm;
                    const double denom = df + (df >= 0.0 ? eps : -eps);
                    const double r = (w >= 0.0)
                        ? db / denom
                        : ((k+2 < km ? m.nh4sh.x[i][j][k+2] : q_pp) - q_pp) / denom;
                    corr += (1.0 - TVD::superbee_phi(r)) * abs(w) * inv_rmsinthe * (df - db) / (2.0 * dphi);
                }

                m.fluxlim_nh4sh.x[i][j][k] = corr;
            }
        }
    }

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for FluxLimiterNH4SH\n", elapsed.count() * 1e-9);
    cout << "      " << Planet::planet_tag() << ": FluxLimiterNH4SH ended" << endl;
}

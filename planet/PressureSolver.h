/*
 * SHARED PHYSICS — Poisson solver for the dynamic pressure, one implementation for every planet.
 * MUST BE BYTE-IDENTICAL IN EVERY MODEL THAT USES IT; `make check-shared` verifies that against
 * planet/SHARED.md5.
 *
 * It knows nothing about which planet it runs on. Beyond the usual fields it asks a model for:
 *
 *     Planet::planet_tag()          "ATJUP" / "ATSAT" — the log prefix, and the prefix its knobs
 *                                   are read under (ATJUP_PDYN_CAP, ATSAT_PRESS_SWEEPS, ...)
 *     Planet::sinthe_min()          the polar metric floor, shared with the momentum equations
 *     Planet::has_obstacle()        whether the model contains any solid body at all
 *     Planet::press_rigid_lid()     default for the rigid radial-wall condition on aux_u
 *     m.is_solid(i,j,k)             whether a cell is inside that body
 *     m.metricRadius(rm)            the metric radius for a cell (see THE METRIC RADIUS below)
 *     m.prepareProjectionBoundaries(rigid_lid)
 *                                   the model's own boundary conditions on the aux_ and rhs_
 *                                   fields before the divergence is taken (see THE BOUNDARY
 *                                   PREPARATION below)
 *
 * ===== WHAT THIS SOLVES, AND THE ONE THING THAT IS NOT SETTLED =====
 *
 * A projection method solves div(grad p) = div(u*), u* being the velocity tendency WITHOUT the
 * pressure gradient — which is what aux_u/aux_v/aux_w hold. That is the default source here
 * (<TAG>_PRESS_SRC=1) and what ATJUP, ATOM_turb and ATOM_Precipitation all use.
 *
 * <TAG>_PRESS_SRC=0 selects the legacy form div(aux - rhs), inherited by ATSAT from old
 * ATOM/atmosphere/Pressure_Atm.cpp. It is kept because it is what ATSAT's computePressure() did
 * for its whole history, and because switching it is the only way to attribute the difference —
 * but it should not be mistaken for an alternative discretisation of the same equation. aux - rhs
 * IS dpdr, identically, by the line that defined aux. So that source is the divergence of the
 * pressure gradient the previous stage already computed, and the equation being relaxed is
 *
 *     Laplacian_compact(p) = Laplacian_wide(p)
 *
 * two discretisations of one operator set equal to each other, which any smooth field satisfies to
 * truncation error. THE VELOCITY DIVERGENCE NEVER ENTERS. With PRESS_SRC=0 the sweep is a
 * checkerboard filter on p_dyn and nothing more, no matter how many sweeps it is given.
 *
 * ===== ONE GAUSS-SEIDEL SWEEP PER CALL IS NOT AN ELLIPTIC SOLVE =====
 *
 * run() is called once per physics iteration and does <TAG>_PRESS_SWEEPS sweeps, default 1. One
 * sweep moves information one cell, so the elliptic problem is never actually solved: p_dyn is a
 * local smoothing of the divergence, not the global pressure response of the flow. That matters
 * most where a body sits in the stream — the high pressure a body builds in front of itself is an
 * elliptic effect reaching several cells upstream, and one sweep per step cannot construct it.
 * Raising the count is the honest way to ask whether a pressure response is merely unconverged;
 * cost is linear in it.
 *
 * ===== THE RELAXATION IS RED-BLACK, AND WAS NOT ALWAYS =====
 *
 * Each sweep is two passes over a checkerboard colouring of (i+j+k): every cell of one colour has
 * all six of its stencil neighbours in the other, so within a pass nothing is read while it is
 * being written. This replaced `#pragma omp parallel for collapse(2)` over (i,j) writing p_dyn in
 * place while reading p_dyn[i±1][j±1] — a data race whose result depended on the thread count and
 * the scheduling.
 *
 * Red-black rather than Jacobi because of what that loop was reaching for: k ran serially inside a
 * thread, so k-1 was current and k+1 one sweep old, i.e. lexicographic Gauss-Seidel — correct in
 * serial, broken only by the (i,j) parallelism. Jacobi would have been the easier fix and would
 * have cost the convergence rate. The colour is selected with a `continue` rather than by striding
 * k, because the k loop carries a sliding window over the solid mask that assumes consecutive k.
 *
 * ===== THE POISSON METRIC, WHICH HAS NO KNOB ON PURPOSE =====
 *
 * The stencil weights are the LAPLACIAN's metric factors, 1/r^2 on theta and 1/(r^2 sin^2) on phi.
 * ATJUP carried ATJUP_POISSON_METRIC because its weights had once been the DIVERGENCE's (1/r and
 * 1/(r sin)); its default was already the corrected form and ATSAT never had the bug. The knob is
 * deliberately NOT carried into shared code: a switch whose only other position is a known-wrong
 * metric does not belong in a file every planet is required to run.
 *
 * Not included either, and also deliberately: the first-derivative parts of the spherical
 * Laplacian, (2/r) dp/dr and (cot/r^2) dp/dtheta. They are small on these grids (2 dr/r ~ 1e-4 of
 * the second derivative; cot(theta) dtheta ~ 0.017 away from the poles) and adding them would
 * break the symmetric seven-point form the max principle below relies on.
 *
 * ===== THE METRIC RADIUS =====
 *
 * Both models can run with the 1/r factors referred to the real planetary radius rather than to
 * rad.z's 1..2, and they implement it differently: ATJUP shifts rad.z itself at initialisation,
 * ATSAT leaves rad.z alone (it is also the stretched coordinate the layer heights sit on) and
 * shifts only the metric factors. So this file asks the model, through m.metricRadius(), and the
 * accessor is the identity on the model that has already done the shift.
 *
 * ===== THE BOUNDARY PREPARATION IS THE MODEL'S, NOT THIS FILE'S =====
 *
 * The conditions imposed on aux_* and rhs_* before the divergence is taken are the two models'
 * OWN and they genuinely differ — ATJUP sets one radial pass of 2-point Neumann on aux only;
 * ATSAT reproduces computePressure()'s three passes of 3-point cubic on aux AND rhs, with polar
 * zeroing and a phi average. Those are boundary conditions on the model's fields, so they stay
 * with the model, behind m.prepareProjectionBoundaries(). What is shared is the elliptic solve.
 * Forcing one form on both would have been the easy thing to write and would have silently
 * changed one planet's boundary conditions.
 */

#pragma once

#include "ATPhys.h"   // ATPhys::env_double / env_int — the per-planet knob prefix

#include <vector>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class Array;

template<class Planet>
class PressureSolver {
public:
    explicit PressureSolver(Planet& model) : m(model) {}

    void run();

private:
    Planet& m;
};

template<class Planet>
void PressureSolver<Planet>::run(){
    using namespace std;
    const char* TAG = Planet::planet_tag();
    cout << endl << "      " << TAG << ": PressureSolver" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // sin(theta) table — depends on j alone. The floor is the SAME one the momentum equations
    // use, so the two halves of the projection cannot disagree about the polar metric; that split
    // is a bug ATJUP had to find in its own pair of files. The exact-zero guard matters only where
    // sinthe_min() is 0, which is ATSAT's default.
    std::vector<double> sinthe_table(m.jm);
    for(int j = 0; j < m.jm; j++){
        sinthe_table[j] = sin(m.the.z[j]);
        if(sinthe_table[j] < Planet::sinthe_min())
            sinthe_table[j] = Planet::sinthe_min();
        if(sinthe_table[j] == 0.0) sinthe_table[j] = 1.0e-5;
    }

    // Divergence-source clamp for the p_dyn update (discrete max principle, from ATOM). Bounds
    // each cell's source contribution so a velocity spike — a steep obstacle flank, or a
    // low-density polar or top cell driven by the radiative coupling — cannot drive p_dyn
    // unbounded and NaN at the converging-meridian pole. A generous backstop, not a working limit.
    const double p_dyn_cap = ATPhys::env_double(TAG, "PDYN_CAP", 0.2);

    // Source term: 1 = div(aux), the projection source; 0 = div(aux - rhs), the legacy form.
    // See the header — the second is not an alternative discretisation of the same equation.
    const int src_projection = ATPhys::env_int(TAG, "PRESS_SRC", 1);

    const int n_sweeps_raw = ATPhys::env_int(TAG, "PRESS_SWEEPS", 1);
    const int n_sweeps = n_sweeps_raw > 0 ? n_sweeps_raw : 1;

    // ---- Solid-wall condition on p_dyn at the obstacle ----
    //
    // Without it the body is INVISIBLE to the pressure equation: the sweep runs over every
    // interior cell, solid ones included, and a fluid cell's stencil reads p_dyn straight out of
    // its solid neighbours, so p is relaxed THROUGH the body as if it were air. A body in a flow
    // works by building pressure in front of itself — that is what turns the stream aside instead
    // of letting it accelerate into the obstacle — and with no wall condition there is nothing to
    // build it against.
    //
    // The condition a rigid impermeable wall imposes on the pressure of a projection method is
    // dp/dn = 0. Discretely, in the finite-volume reading of this 7-point stencil, that is exactly
    // "no pressure flux through a solid face": the face is dropped from the numerator AND the
    // denominator rather than being fed a value from inside the body. Only FLUID cells are
    // treated; solid cells keep their one-sided source and their bcSolidGround values.
    //
    // It defaults to has_obstacle(), so a model with no solid body takes the plain six-face
    // branch. That is not a fudge for bit-identity: with nothing solid the two branches are the
    // same sum, and the plain one says so in one expression instead of six conditionals.
    const bool press_wall = ATPhys::env_int(TAG, "PRESS_WALL", Planet::has_obstacle() ? 1 : 0) != 0;

    // Rigid radial walls. aux_u is the wall-NORMAL intermediate velocity and it feeds du_dr in the
    // divergence source, so it carries whatever condition the radial walls are meant to impose.
    // Extrapolating it re-injects a wall-normal flux into the projection and leaves the column
    // mass budget open; zeroing it closes it, which is what a rigid lid and a rigid deep boundary
    // mean. The default is the model's, because the models disagree for a real reason: ATJUP's
    // i=0 is a floor, ATSAT's is the deep interior of a gas giant and whether a lid belongs there
    // at all is a modelling question no port gets to settle.
    const bool rigid_lid = ATPhys::env_int(TAG, "BC_RIGID_LID",
                                           Planet::press_rigid_lid() ? 1 : 0) != 0;

    // The model's own conditions on aux_* and rhs_*, before any divergence is taken.
    m.prepareProjectionBoundaries(rigid_lid);

    // Grid-spacing reciprocals — constant for the entire grid.
    const double inv_2dr   = 1.0 / (2.0 * m.dr);
    const double inv_2dthe = 1.0 / (2.0 * m.dthe);
    const double inv_2dphi = 1.0 / (2.0 * m.dphi);
    const double inv_dr2   = 1.0 / (m.dr   * m.dr);
    const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
    const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);
    const double inv_dthe  = 1.0 / m.dthe;
    const double inv_dphi  = 1.0 / m.dphi;

    for(int sweep = 0; sweep < n_sweeps; sweep++){

      for(int colour = 0; colour < 2; colour++){

        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < m.im-1; i++){
            for(int j = 1; j < m.jm-1; j++){

                // Geometry once per (i,j).
                typename Planet::CellGeometry geo;
                geo.rm       = m.metricRadius(m.rad.z[i]);
                geo.rm2      = geo.rm * geo.rm;
                geo.exp_rm   = m.coord_stretching ? 1.0 / (geo.rm + 1.0) : 1.0;
                geo.exp_2_rm = geo.exp_rm * geo.exp_rm;
                geo.sinthe   = sinthe_table[j];
                geo.sinthe2  = geo.sinthe * geo.sinthe;
                geo.costhe   = cos(m.the.z[j]);
                geo.cotanthe             = geo.costhe / geo.sinthe;
                geo.inv_rm               = 1.0 / geo.rm;
                geo.inv_rm2              = 1.0 / geo.rm2;
                geo.inv_rmsinthe         = 1.0 / (geo.rm * geo.sinthe);
                geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
                geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
                geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
                geo.inv_2dr   = inv_2dr;
                geo.inv_2dthe = inv_2dthe;
                geo.inv_2dphi = inv_2dphi;
                geo.inv_dr2   = inv_dr2;
                geo.inv_dthe2 = inv_dthe2;
                geo.inv_dphi2 = inv_dphi2;

                // Stencil weights — the Laplacian's own metric factors. See the header.
                const double num1 = geo.exp_2_rm       * inv_dr2;
                const double num2 = geo.inv_rm2        * inv_dthe2;
                const double num3 = geo.inv_rm2sinthe2 * inv_dphi2;
                const double denom = 2.0 * num1 + 2.0 * num2 + 2.0 * num3;

                const bool i_in_range = (i < m.im-2);
                const bool j_inner    = (j > 2) && (j < m.jm-2);

                // Sliding window for the k-direction solid status. It must advance on every k,
                // which is why the colour is skipped with a `continue` further down rather than
                // by striding this loop.
                bool sld_k0 = m.is_solid(i,j,0), sld_k1 = m.is_solid(i,j,1);

                for(int k = 1; k < m.km-1; k++){
                    const bool sld_ijk = sld_k1;
                    const bool sld_kp1 = m.is_solid(i,j,k+1);
                    const bool sld_km1 = sld_k0;
                    const bool sld_kp2 = (k < m.km-2) ? m.is_solid(i,j,k+2) : false;
                    const bool sld_km2 = (k > 2)      ? m.is_solid(i,j,k-2) : false;

                    sld_k0 = sld_k1;
                    sld_k1 = sld_kp1;

                    // Cells of the other colour are skipped AFTER the window bookkeeping above.
                    if(((i + j + k) & 1) != colour) continue;

                    // ---- Divergence source ----
                    double du_dr, dv_dthe, dw_dphi;

                    if(src_projection == 0){
                        // Legacy form: div(aux - rhs). Plain centred differences; the one-sided
                        // obstacle stencils below belong to the projection source and are not
                        // applied to a source that is not a velocity divergence at all.
                        du_dr   = ((m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k])
                                 - (m.rhs_u.x[i+1][j][k] - m.rhs_u.x[i-1][j][k])) * inv_2dr;
                        dv_dthe = ((m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k])
                                 - (m.rhs_v.x[i][j+1][k] - m.rhs_v.x[i][j-1][k])) * inv_2dthe;
                        dw_dphi = ((m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1])
                                 - (m.rhs_w.x[i][j][k+1] - m.rhs_w.x[i][j][k-1])) * inv_2dphi;
                    } else {
                        bool r_flag = false, the_flag = false, phi_flag = false;

                        // r direction
                        if(i_in_range && sld_ijk && !m.is_solid(i+1,j,k)){
                            du_dr = (-3.0 * m.aux_u.x[i][j][k] + 4.0 * m.aux_u.x[i+1][j][k]
                                     - m.aux_u.x[i+2][j][k]) * inv_2dr;
                            r_flag = true;
                        }

                        // theta direction
                        if(j_inner){
                            const bool air_jp1 = !m.is_solid(i,j+1,k);
                            const bool air_jm1 = !m.is_solid(i,j-1,k);

                            if(sld_ijk && air_jp1 && !m.is_solid(i,j+2,k)){
                                dv_dthe = (-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j+1][k]
                                           - m.aux_v.x[i][j+2][k]) * inv_2dthe;
                                the_flag = true;
                            }
                            if(sld_ijk && air_jm1 && !m.is_solid(i,j-2,k)){
                                dv_dthe = -(-3.0 * m.aux_v.x[i][j][k] + 4.0 * m.aux_v.x[i][j-1][k]
                                            - m.aux_v.x[i][j-2][k]) * inv_2dthe;
                                the_flag = true;
                            }
                            if((sld_ijk && air_jp1 && m.is_solid(i,j+2,k))
                               || (j == m.jm-2 && !sld_ijk && m.is_solid(i,j+1,k))){
                                dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                                the_flag = true;
                            }
                            if((sld_ijk && air_jm1 && m.is_solid(i,j-2,k))
                               || (j == 1 && sld_ijk && air_jm1)){
                                dv_dthe = (m.aux_v.x[i][j-1][k] - m.aux_v.x[i][j][k]) * inv_dthe;
                                the_flag = true;
                            }
                        }

                        // phi direction
                        const bool k_inner = (k > 2) && (k < m.km-2);
                        if(k_inner){
                            const bool air_kp1 = !sld_kp1;
                            const bool air_km1 = !sld_km1;

                            if(sld_ijk && air_kp1 && !sld_kp2){
                                dw_dphi = (-3.0 * m.aux_w.x[i][j][k] + 4.0 * m.aux_w.x[i][j][k+1]
                                           - m.aux_w.x[i][j][k+2]) * inv_2dphi;
                                phi_flag = true;
                            }
                            if(sld_ijk && air_km1 && !sld_km2){
                                dw_dphi = -(-3.0 * m.aux_w.x[i][j][k] + 4.0 * m.aux_w.x[i][j][k-1]
                                            - m.aux_w.x[i][j][k-2]) * inv_2dphi;
                                phi_flag = true;
                            }
                            if((sld_ijk && air_kp1 && sld_kp2)
                               || (k == m.km-2 && !sld_ijk && sld_kp1)){
                                dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k]) * inv_dphi;
                                phi_flag = true;
                            }
                            if((sld_ijk && air_km1 && sld_km2)
                               || (k == 1 && sld_ijk && air_km1)){
                                dw_dphi = (m.aux_w.x[i][j][k-1] - m.aux_w.x[i][j][k]) * inv_dphi;
                                phi_flag = true;
                            }
                        }

                        // Centred-difference fallbacks.
                        if(!r_flag)
                            du_dr   = (m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k]) * inv_2dr;
                        if(!the_flag)
                            dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k]) * inv_2dthe;
                        if(!phi_flag)
                            dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1]) * inv_2dphi;
                    }

                    double div_src = du_dr   * geo.exp_rm
                                   + dv_dthe * geo.inv_rm
                                   + dw_dphi * geo.inv_rmsinthe;

                    // Assemble the stencil, dropping any face that looks into the body. For a
                    // model with no obstacle press_wall is off and this is the plain six-face sum.
                    double acc = 0.0, den = denom;
                    if(press_wall && !sld_ijk){
                        acc = den = 0.0;
                        if(!m.is_solid(i+1,j,k)){ acc += num1 * m.p_dyn.x[i+1][j][k]; den += num1; }
                        if(!m.is_solid(i-1,j,k)){ acc += num1 * m.p_dyn.x[i-1][j][k]; den += num1; }
                        if(!m.is_solid(i,j+1,k)){ acc += num2 * m.p_dyn.x[i][j+1][k]; den += num2; }
                        if(!m.is_solid(i,j-1,k)){ acc += num2 * m.p_dyn.x[i][j-1][k]; den += num2; }
                        if(!sld_kp1)            { acc += num3 * m.p_dyn.x[i][j][k+1]; den += num3; }
                        if(!sld_km1)            { acc += num3 * m.p_dyn.x[i][j][k-1]; den += num3; }
                        // A fluid cell walled in on all six faces has no pressure equation left;
                        // leave it as it is rather than dividing by zero.
                        if(den <= 0.0) continue;
                    } else {
                        acc = (m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                            + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                            + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3;
                    }

                    // The clamp is a discrete max principle, so it uses the denominator this cell
                    // actually ends up with, not the full six-face one. The non-finite guard stops
                    // p_dyn running away to NaN at the pole.
                    const double src_max = den * p_dyn_cap;
                    if      (!std::isfinite(div_src)) div_src = 0.0;
                    else if (div_src >  src_max)      div_src =  src_max;
                    else if (div_src < -src_max)      div_src = -src_max;

                    m.p_dyn.x[i][j][k] = (acc - div_src) / den;
                } // k
            } // j
        } // i
      } // colour

        // ---- Outer boundary condition of the relaxation ----
        // 2-point Neumann, not the 3-point cubic: the cubic amplifies an alternating error 7x per
        // call, and p_dyn is not touched by either model's boundary-condition pass, so whatever
        // lands at i=0 here survives the whole following RK step and drives dpdr at i=1.
        //
        // These passes belong INSIDE the sweep loop: they are the outer boundary condition of the
        // relaxation, and leaving them outside would let the interior run away from its own edges
        // for n_sweeps-1 passes.
        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++){
            for(int j = 0; j < m.jm; j++){
                m.p_dyn.x[0][j][k]      = m.c43 * m.p_dyn.x[1][j][k]      - m.c13 * m.p_dyn.x[2][j][k];
                m.p_dyn.x[m.im-1][j][k] = m.c43 * m.p_dyn.x[m.im-2][j][k] - m.c13 * m.p_dyn.x[m.im-3][j][k];
            }
        }

        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++){
            for(int i = 0; i < m.im; i++){
                m.p_dyn.x[i][0][k]      = m.c43 * m.p_dyn.x[i][1][k]      - m.c13 * m.p_dyn.x[i][2][k];
                m.p_dyn.x[i][m.jm-1][k] = m.c43 * m.p_dyn.x[i][m.jm-2][k] - m.c13 * m.p_dyn.x[i][m.jm-3][k];
            }
        }

        #pragma omp parallel for collapse(2)
        for(int i = 0; i < m.im; i++){
            for(int j = 0; j < m.jm; j++){
                m.p_dyn.x[i][j][0]      = m.c43 * m.p_dyn.x[i][j][1]      - m.c13 * m.p_dyn.x[i][j][2];
                m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2] - m.c13 * m.p_dyn.x[i][j][m.km-3];
                m.p_dyn.x[i][j][0] = m.p_dyn.x[i][j][m.km-1]
                    = (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1]) / 2.0;
            }
        }

    } // sweep

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for PressureSolver\n", elapsed.count() * 1e-9);
    cout << "      " << TAG << ": PressureSolver ended" << endl;
}

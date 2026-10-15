/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone boundary-condition class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * This is a header-only file: it contains the BC_Uran class, all inline method
 * bodies, and the inline cUranusModel delegation wrappers that forward each
 * cUranusModel::BC_xxx() call to the corresponding BC_Uran method.
 * BC_Uran.cpp is NOT a stub, whatever this comment used to claim: it defines the four
 * bc_fields_* lists — radius, phi, theta_extrap and theta_zero — naming the ~25 arrays each
 * boundary pass applies to. That list is Uranus's own chemistry, and it is precisely the part
 * BoundaryConditions<Planet> cannot carry if it is to stay byte-identical across four planets
 * with different species. A file called a stub is a file someone eventually deletes.
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
#include "BoundaryConditions.h"
#include "Utils.h"

using namespace AtomUtils;


/*
 * The three boundary passes are the SHARED BoundaryConditions<Planet> now. What used to be here —
 * three field lists walked with a (4/3,-1/3) extrapolation — is the same algorithm the other three
 * models run; the lists moved to cUranusModel (BC_Uran.cpp) because they are Uranus's, and the two
 * places this model genuinely differs are named there as bc_margin() = 0 and
 * bc_default_form() = NEUMANN.
 *
 * BC_Uran stays as the name so no call site changes, and because initTropopauseLayers below is
 * ATURAN's own and has no counterpart in the shared header.
 */
// ---------------------------------------------------------------------------
// The two species floors of ATJUP, ported 2026-10-09. BOTH DEFAULT 1 (the user's word, same
// day); 0 restores each. Until then this model had neither: MEASURED at the defaults, 16 iterations (ATJUP/satchk/
// giants/URANa), several species are negative, most of them on the model-top plane -- which is
// the 2-point radial form f[s] = (4/3)f[a] - (1/3)f[b] going below zero wherever f[a] < f[b]/4 --
// and some in the interior, which is the transport.
//
//   ATURAN_BC_RADIUS_POSITIVE   1 = floor the radial boundary planes of the species at zero where
//                              they are written; 2 = fall back to the plain copy f[s] = f[a] there
//                              (still a zero-gradient wall, first order in those cells)
//   ATURAN_SPECIES_CLAMP        1 = set every negative species value to zero once per iteration,
//                              after the boundary conditions and before restoreVar, and report
//                              what that added per field at the end of the run
//
// 16 iterations, 3 threads (ATJUP/satchk/giants/URAN e, f, g): the boundary floor alone removes
// every top-plane negative and leaves the interior ones; with both no species is negative at
// any print. Against neither, the printed extrema that move do so in the sixth digit. NOT
// measured: anything longer than 16 iterations.
// ---------------------------------------------------------------------------
inline std::vector<Array*> cUranusModel::species_fields(){
    return { &h2o, &h2o_cloud, &h2o_ice,
             &h2s, &h2s_cloud, &h2s_ice,
             &nh3, &nh3_cloud, &nh3_ice,
             &ch4, &ch4_cloud, &ch4_ice,
             &nh4sh };
}

inline void cUranusModel::floorRadialSpecies(){
    static const int positive = [](){ const char* e = getenv("ATURAN_BC_RADIUS_POSITIVE"); return e ? atoi(e) : 1; }();
    if(positive == 0) return;
    std::vector<Array*> sp = species_fields();
    const int ns = (int)sp.size();
    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int f = 0; f < ns; f++){
                Array& F = *sp[f];
                if(F.x[0][j][k] < 0.0)
                    F.x[0][j][k]    = (positive == 2) ? std::max(0.0, F.x[1][j][k])    : 0.0;
                if(F.x[im-1][j][k] < 0.0)
                    F.x[im-1][j][k] = (positive == 2) ? std::max(0.0, F.x[im-2][j][k]) : 0.0;
            }
        }
    }
}

inline void cUranusModel::clampNegativeSpecies(){
    static const int on = [](){ const char* e = getenv("ATURAN_SPECIES_CLAMP"); return e ? atoi(e) : 1; }();
    if(on == 0) return;
    std::vector<Array*> sp = species_fields();
    if(clamp_added.size() != sp.size()){ clamp_added.assign(sp.size(), 0.0); clamp_cells.assign(sp.size(), 0); }
    for(size_t f = 0; f < sp.size(); f++){
        Array& F = *sp[f];
        double a = 0.0; long c = 0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:a,c)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++)
                    if(!(F.x[i][j][k] >= 0.0)){               // catches a NaN too
                        if(std::isfinite(F.x[i][j][k])){ a -= F.x[i][j][k]; c++; }
                        F.x[i][j][k] = 0.0;
                    }
        clamp_added[f] += a; clamp_cells[f] += c;
    }
}

inline void cUranusModel::clampNegativeReport(){
    if(clamp_added.empty()) return;
    static const char* const n[] = { "h2o", "h2o_cloud", "h2o_ice", "h2s", "h2s_cloud", "h2s_ice",
        "nh3", "nh3_cloud", "nh3_ice", "ch4", "ch4_cloud", "ch4_ice", "nh4sh" };
    printf("\n      ATURAN: negative-value clamp, cumulative since start (sum of the clipped cell values, field units)\n");
    for(size_t f = 0; f < clamp_added.size(); f++)
        if(clamp_cells[f] > 0)
            printf("        %-10s  gross %.4e over %10ld clippings\n", n[f], clamp_added[f], clamp_cells[f]);
}

inline void BC_Uran::bcRadius(){
    BoundaryConditions<cUranusModel>(m).bcRadius();
    m.floorRadialSpecies();
    // ATURAN_BC_LID_VW=<form> (2026-10-15, DEFAULT 2 WITH ATURAN_GEO_BALANCE, else 0 = the shared form, Neumann
    // w[im-1] = 4/3 w[im-2] - 1/3 w[im-3], i.e. dw/dr = 0): the lid value of the TANGENTIAL
    // velocities v, w only; t, the species and u keep theirs. 2 = linear extrapolation
    // 2 a - b (zero curvature), 3 = the cubic c - 3 b + 3 a. The balanced start carries the
    // thermal-wind shear of the isothermal top up to the lid; dw/dr = 0 there gives the cell
    // under the lid a curvature of -2/3 of its shear step, and the molecular diffusion of w
    // takes it out (79 N, level 39: -0.085 nondimensional at iteration 1, UDBGW48), spinning
    // the top down against the frozen p_baro. ATOM_Precipitation rejected the cubic lid for
    // v, w that decay to zero at the top (it overshoots through zero); here they do not.
    // MEASURED (224 it with PDYN_GEO_BC and DAMP_T_J=0; UGBC224 -> ULVW2_224): residuum 0.000464 ->
    // 0.000302 (old defaults UGA 0.000297), max |v| 0.088 -> 0.047, the lid jet held (max w
    // 214.895 -> 214.891 against 213.97 -> 213.11), T and OLR unchanged; 3 the same to 1 %.
    static const int lid_vw = [](){ const char* e = getenv("ATURAN_BC_LID_VW"); return e ? atoi(e) : (cUranusModel::geo_balance() ? 2 : 0); }();
    if(lid_vw == 2 || lid_vw == 3){
        const int a = m.im - 2, b = m.im - 3, c = m.im - 4, l = m.im - 1;
      #pragma omp parallel for
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                if(lid_vw == 2){
                    m.v.x[l][j][k] = 2.0 * m.v.x[a][j][k] - m.v.x[b][j][k];
                    m.w.x[l][j][k] = 2.0 * m.w.x[a][j][k] - m.w.x[b][j][k];
                }else{
                    m.v.x[l][j][k] = m.v.x[c][j][k] - 3.0 * m.v.x[b][j][k] + 3.0 * m.v.x[a][j][k];
                    m.w.x[l][j][k] = m.w.x[c][j][k] - 3.0 * m.w.x[b][j][k] + 3.0 * m.w.x[a][j][k];
                }
            }
        }
    }
}
inline void BC_Uran::bcTheta()  { BoundaryConditions<cUranusModel>(m).bcTheta();  }
inline void BC_Uran::bcPhi()    { BoundaryConditions<cUranusModel>(m).bcPhi();    }


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

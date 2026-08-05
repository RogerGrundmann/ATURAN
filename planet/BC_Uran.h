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
inline void BC_Uran::bcRadius() { BoundaryConditions<cUranusModel>(m).bcRadius(); }
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

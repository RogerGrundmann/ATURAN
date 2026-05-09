/*
 * cUranusModel delegation wrappers.
 * Each cUranusModel::Foo() method declared in cUranusModel.h but implemented
 * in a standalone friend class is forwarded here.
 * Keeping the includes in one translation unit avoids circular-include issues
 * in cUranusModel.cpp.
*/

#include "cUranusModel.h"
#include "BC_Uran.h"
#include "ChemistryUran.h"
#include "SaturationAdjustmentUran.h"

// ---- Boundary conditions --------------------------------------------------

void cUranusModel::BC_radius() { BC_Uran(*this).bcRadius(); }
void cUranusModel::BC_theta()  { BC_Uran(*this).bcTheta();  }
void cUranusModel::BC_phi()    { BC_Uran(*this).bcPhi();    }

void cUranusModel::init_tropopause_layers() { BC_Uran(*this).initTropopauseLayers(); }

// ---- Chemistry ------------------------------------------------------------

void cUranusModel::ChemMassRateUran()      { ChemistryUran(*this).ChemMassRateUran();      }
void cUranusModel::DiffMassFluxUran()      { ChemistryUran(*this).DiffMassFluxUran();      }
void cUranusModel::ThermalPropertiesUran() { ChemistryUran(*this).ThermalPropertiesUran(); }

// ---- Saturation adjustment ------------------------------------------------

void cUranusModel::Saturation_Adjustment(std::string gas,
    double &coeff_A,   double &coeff_B,
    double &coeff_A_i, double &coeff_B_i,
    double &t_0,       double &t_00,
    double &ep,        double &lv,  double &ls,
    double &cp,        double &r,
    double &C,         double &L0,  double &R,
    double &del_alf,   double &del_bet, double &m,
    Array &c, Array &cloud, Array &ice)
{
    SaturationAdjustmentUran(*this).run(gas,
        coeff_A, coeff_B, coeff_A_i, coeff_B_i,
        t_0, t_00, ep, lv, ls, cp, r,
        C, L0, R, del_alf, del_bet, m,
        c, cloud, ice);
}

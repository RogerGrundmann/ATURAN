/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone saturation-adjustment and microphysics class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * Reference algorithm:
 *   Tao, W.-K., Simpson, J., and McCumber, M.:
 *   "An Ice-Water Saturation Adjustment", AMS Notes and Correspondence, 1988.
*/

#pragma once

#include <cmath>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <string>
#include <iostream>

#include "SaturationAdjustment.h"

// cUranusModel.h is included by the translation unit that uses this header.
// We only need the forward declaration here to avoid a circular include.
class cUranusModel;
class Array;

using namespace std;

class SaturationAdjustmentUran {
public:

    explicit SaturationAdjustmentUran(cUranusModel& model) : m(model) {}

    // Selects between ATURAN's inherited routine and the SHARED SaturationAdjustment<Planet>
    // that ATSAT, ATJUP and ATNEPT instantiate. Default 0 = the inherited one, so the model is
    // byte-identical until this is set. ATURAN_SATADJ=1 selects the shared algorithm.
    //
    // NOT YET EVALUATED ON URANUS. On ATSAT the switch was measured to condense 16 % less peak
    // cloud water and move the deck a layer higher, and which of the two is right was not settled
    // there or on Neptune. The same caveat applies here and one more besides: the ice-phase
    // coefficient quadruple passed in the .cpp is the LIQUID one, because ATURAN's parameter set
    // has no ice pair for H2O, NH3 or CH4 — and Uranus is the coldest of the four, so the ice
    // branch is the one that matters most. Supplying real ice coefficients is what would make
    // this knob worth turning on.
    static int mirrored_enabled(){
        static const int v = [](){
            const char* e = getenv("ATURAN_SATADJ"); return e ? atoi(e) : 0; }();
        return v;
    }

    // Dispatch. The call-site signature is unchanged, so cUranusModel.cpp needs no edit.
    void run(const std::string& gas,
             double coeff_A,   double coeff_B,
             double coeff_A_i, double coeff_B_i,
             double t_0,       double t_00,
             double ep,        double lv,  double ls,
             double cp,        double r,
             double C,         double L0,  double R,
             double del_alf,   double del_bet,   double m_mol,
             Array& c,         Array& cloud,   Array& ice);

    // ATURAN's own algorithm, unchanged and still the default path.
    void run_legacy(const std::string& gas,
             double coeff_A,   double coeff_B,
             double coeff_A_i, double coeff_B_i,
             double t_0,       double t_00,
             double ep,        double lv,  double ls,
             double cp,        double r,
             double C,         double L0,  double R,
             double del_alf,   double del_bet,   double m_mol,
             Array& c,         Array& cloud,   Array& ice);

    // -----------------------------------------------------------------------
    // Static helper functions — usable without a SaturationAdjustmentUran instance
    // -----------------------------------------------------------------------
    static double clausius_clapeyron(double T_K, double A, double B){
        return std::exp(A / T_K + B);
    }

    static double saturation_vapour_pressure(double T_K,
            double C, double L0, double R, double del_alf, double del_bet){
        return std::exp(C
            + (-L0 / T_K + del_alf * std::log(T_K) + del_bet * T_K)
            / (1e-3 * R));
    }

    static double humility_critical(double x, double Hu_cr_max, double Hu_cr_mid){
        return (Hu_cr_max - Hu_cr_mid) * (x * x - 2.0 * x) + Hu_cr_max;
    }

private:
    cUranusModel& m;

    static constexpr int iter_prec_end = 30;
    static constexpr double q_diff_min = 1.0e-4;
};

// Implementation is in SaturationAdjustmentUran.cpp

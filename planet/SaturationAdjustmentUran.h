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

// cUranusModel.h is included by the translation unit that uses this header.
// We only need the forward declaration here to avoid a circular include.
class cUranusModel;
class Array;

using namespace std;

class SaturationAdjustmentUran {
public:

    explicit SaturationAdjustmentUran(cUranusModel& model) : m(model) {}

    void run(const std::string& gas,
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

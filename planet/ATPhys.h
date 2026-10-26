/*
 * SHARED PHYSICS — small helpers every planet model's physics needs. MUST BE BYTE-IDENTICAL IN
 * EVERY MODEL THAT USES IT; `make check-shared` verifies that against planet/SHARED.md5.
 *
 * Nothing here knows which planet it is running on. The environment helpers take the model's own
 * planet_tag() so one shared implementation reads each model's knobs under its own prefix, and
 * the saturation vapour pressure is the single expression that ATJUP, ATSAT and their saturation
 * adjustments were each computing separately — character for character the same formula in three
 * places, which is three places for it to drift.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

namespace ATPhys {


// Look up "<PLANET>_<name>", e.g. ATJUP_CONV_ADJ_LAPSE, so one shared implementation reads each
// model's own knobs under its own prefix. Not cached: the caller caches where it matters.
inline const char* env_for(const char* planet_tag, const char* name){
    std::string key(planet_tag);
    key += "_";
    key += name;
    return getenv(key.c_str());
}

inline double env_double(const char* planet_tag, const char* name, double fallback){
    const char* e = env_for(planet_tag, name);
    return e ? atof(e) : fallback;
}

inline int env_int(const char* planet_tag, const char* name, int fallback){
    const char* e = env_for(planet_tag, name);
    return e ? atoi(e) : fallback;
}


// Saturation vapour pressure [bar] from the Sanchez-Lavega / Clausius-Clapeyron coefficient set.
//
// This is the expression SaturationAdjustmentJup::saturation_vapour_pressure,
// SaturationAdjustmentSat::saturation_vapour_pressure and cSaturnModel::saturation_vapour_pressure
// all evaluated, in that order of discovery and with identical operand order. All four copies are
// gone now and every caller in both models comes here; R arrives in J/(kg K) and the 1e-3 takes it
// to the J/(g K) the coefficients are tabulated against.
inline double saturation_vapour_pressure(double T_K, double C, double L0, double R,
                                         double del_alf, double del_bet){
    return std::exp(C + (-L0 / T_K + del_alf * std::log(T_K) + del_bet * T_K) / (1e-3 * R));
}


// FITTED VAPOUR-PRESSURE SETS (2026-10-26), in saturation_vapour_pressure's own form, so a model swaps
// a set by assigning four members and every consumer follows. Properties of the substance, not of a
// planet: one copy here for all four models (<TAG>_SVP_FIT, default 0 in each). Sources and checks are
// written out at ATJUP's setSvpFit (InitValues_Jup.cpp); R is the model's own specific gas constant of the
// gas (R_h2o 461.52, R_nh3 488.21 J/(kg K)), against which the fits were made.
//   H2O ice     Murphy & Koop (2005), exact (the formula has this form)
//   NH3 ice     least-squares fit to Fray & Schmitt (2009), 110-195.4 K, 0.05 %
//   NH3 liquid  least-squares fit to the NIST Antoine equations, 195.4-330 K, 3.2 % (5 % at 350 K);
//               below 195.4 K it behaves like the shipped liquid set (within a factor 2 to 100 K)
struct SvpSet { double C, L0, del_alf, del_bet; };
inline SvpSet svp_fit_h2o_ice(){ return { -1.962499,  2641.4013,   1.629479,  -3.361398e-3 }; }
inline SvpSet svp_fit_nh3_ice(){ return { 16.086062,  1857.9898,   0.098768,  -1.179989e-3 }; }
inline SvpSet svp_fit_nh3_liq(){ return { 196.077451, 3586.5634, -15.978352,   2.837146e-2 }; }


// Two-coefficient Clausius-Clapeyron, E = exp(A/T + B). The short form used where a full
// Sanchez-Lavega coefficient set is not tabulated; Thermo_Jup.cpp is its only caller.
//
// It arrives here from SaturationAdjustmentJup and SaturationAdjustmentSat, which held it as
// character-for-character identical statics. ATSAT's copy had NO caller at all — it existed only
// because the class was written to mirror ATJUP's — so this is one live implementation replacing
// one live and one ornamental copy.
inline double clausius_clapeyron(double T_K, double A, double B){
    return std::exp(A / T_K + B);
}

// The floor applied to sin(theta) by terms that DIVIDE by (r*sin(theta))^2 — the TVD flux limiter
// and the spherical Laplacian in the chemistry. Without one they are infinite at the pole.
//
// THIS USED TO BE THE LITERAL 0.4, WRITTEN OUT AT FOUR SITES, and the comment beside two of them
// said "matches sinthe_min in RungeKutta". It did, once: ATJUP's metric floor WAS 0.4 before it was
// raised to 0.55 to stop a long-run polar blow-up, and the copies never followed. They are fossils
// of a single concept that drifted, not a second concept — which is why this is one accessor and
// not four constants.
//
//   <TAG>_SINTHE_TRACK=0  (default)  the historical literal 0.4, i.e. unchanged in both models
//   <TAG>_SINTHE_TRACK=1             follow the model's own metric floor, max(sinthe_min(), 0.4)
//
// WHAT THAT DOES TO EACH MODEL, which is not symmetric:
//   ATSAT  sinthe_min() defaults to 0.0, so max(0.0, 0.4) = 0.4 and the knob is INERT. It only
//          begins to matter if ATSAT adopts a metric floor, which is the open question — ATSAT's
//          RK integrates to j=2 where 1/sin^2 reaches 820, against the 3.3 ATJUP allows itself.
//   ATJUP  sinthe_min() is 0.55, so =1 moves these terms 0.4 -> 0.55 and makes the chemistry agree
//          with the momentum equations for the first time since the floor was raised. That CHANGES
//          ATJUP's results, which is why it is off by default rather than simply corrected.
//
// The 0.4 lower bound is kept under the max() deliberately: it guarantees the divisors stay finite
// even for a model that declares no metric floor at all, which is exactly ATSAT's position.
template<class Planet>
inline double polar_divisor_floor(){
    static const int    track = env_int(Planet::planet_tag(), "SINTHE_TRACK", 0);
    static const double v     = track ? std::max(Planet::sinthe_min(), 0.4) : 0.4;
    return v;
}


// ===== THE PHI SEAM: <TAG>_BC_PHI_PERIODIC =====
//
// The longitude grid has km points over 360 degrees with BOTH ends on it: k = 0 and k = km-1 are
// the same meridian. Every loop that differences in phi runs k = 1..km-2, so that meridian is
// integrated by nobody: bcPhi sets it to the mean of two one-sided extrapolations, the pressure
// solvers do the same to p_dyn, and the fields no list carries keep whatever was there.
//
// With the knob at 1 the meridian is an ordinary cell. The loops run k = 0..km-2 and take their
// phi neighbours through phi_wrap(), which closes the ring over the km-1 distinct meridians
// (the west neighbour of 0 is km-2, the east neighbour of km-2 is 0); k = km-1 is then only the
// copy of k = 0 that the output and the pointwise modules read, written by bcPhi, by the pressure
// solver and by the Shapiro filter. With the knob at 0 phi_wrap() is the identity and
// phi_first() is 1: the loops are the ones they were.
template<class Planet>
inline bool phi_periodic(){
    static const bool v = env_int(Planet::planet_tag(), "BC_PHI_PERIODIC",
                                  Planet::bc_default_phi_periodic()) != 0;
    return v;
}

// First k of a loop that differences in phi; the last one is km-2 in both cases.
inline int phi_first(bool periodic){ return periodic ? 0 : 1; }

// The index of a phi neighbour kk = k-2 .. k+2 of a cell k.
inline int phi_wrap(bool periodic, int kk, int km){
    if(!periodic) return kk;
    const int n = km - 1;
    return kk < 0 ? kk + n : (kk >= n ? kk - n : kk);
}

}  // namespace ATPhys

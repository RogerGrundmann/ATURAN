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

}  // namespace ATPhys

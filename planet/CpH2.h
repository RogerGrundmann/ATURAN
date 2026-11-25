/*
 * CpH2.h -- SHARED (byte-identical in ATJUP / ATSAT / ATNEPT / ATURAN, planet/SHARED.md5)
 *
 * <TAG>_CP_H2 = 0 | 1 | 2  and  <TAG>_CP_GAS = 0 | 1  (2026-11-24, row GIANTS-ICE-TP; user: "do both, start with the
 * T-dependent cp_h2", "yes, implement it this way with the 1-bar anchoring"). Knobs 0 / 0 = the constant cp_mix, every
 * written field as before. DEFAULTS since 2026-11-24 (user: "make the cp(T) knobs default on all four"; per model hooks
 * Planet::cp_h2_default() / cp_gas_default()): ATSAT CP_H2 1, ATURAN CP_GAS 1 + CP_H2 2, ATNEPT CP_GAS 1 + CP_H2 1,
 * ATJUP 0 / 0 (set-up not ported, the knobs are refused there).
 *
 * WHY. Every model's lapse rate is g / cp_mix with cp_h2 = 14.32e3 J/(kg K), the 298 K value of H2. Below ~250 K H2's
 * rotational heat capacity is smaller (normal H2, ortho:para frozen at 3:1: 10.6 kJ/(kg K) at 80 K, 12.5 at 150 K) or,
 * with ortho-para conversion (equilibrium H2), larger (15.7 at 80 K). And on the ice giants cp_mix is mass-weighted over
 * the SET-UP densities, ~40 % of which are low-cp condensables (CH4 0.19, H2S 0.3, NH3 0.1, H2O 0.08 of r_mix 1.677),
 * while the pressure behaves like H2-He. MEASURED (literature/voyager_rso/tp_compare.*, model vs Voyager 2 radio
 * occultation + Atreya's ECCM T-p): Saturn DL 5000 -8 / -20 / -16 K at 2 / 5 / 10 bar (slope R/cp 0.29 vs 0.30-0.33),
 * Uranus DL1 1000 +5 / +17 / +38 K (0.41 vs 0.31-0.36). PREDICTED offline (literature/cp_h2/predict.out, dlnT/dlnp =
 * R/cp(T) from the Voyager 1-bar T): Saturn with normal H2 within ~1 K of Atreya & Wong 2005 Fig. 3 to 10 bar; Uranus
 * with CP_GAS (+ equilibrium H2) within ~5 K; Neptune with CP_GAS + normal H2 within ~2 K of Voyager to 5 bar.
 *
 * WHAT. cp_h2(T) from the rotational level sum of H2 (theta_r = 87.57 K, theta_v = 6332 K): cp/R = 5/2 + c_rot + c_vib,
 * normal = 3/4 c_ortho + 1/4 c_para (each over its own levels), equilibrium = one sum with nuclear weights 1 (even J) /
 * 3 (odd J). Tabulated once at 0.5 K steps from 5 to 1500 K, linear in between, clamped outside; 298.15 K gives 28.69
 * J/(mol K) (tabulated 28.84; rigid rotor). The model's cp for lapse rates then is
 *   CP_GAS 0:  cp(T) = cp_mix + (r_h2 / r_mix) (cp_h2(T) - cp_h2)                 (only H2's share changes)
 *   CP_GAS 1:  cp(T) = (r_h2 cp_h2(T) + r_he cp_he + r_ch4 cp_ch4) / (r_h2 + r_he + r_ch4)   (the gas, no condensables)
 * with cp_h2(T) = cp_h2 (the model's constant) when CP_H2 is 0. enth(T) is its integral (any constant), so that energy
 * sums written as cp T become H(T). With both knobs 0 cp_lapse returns cp_mix and enth returns cp_mix * T exactly.
 * Used where the lapse rate is set: the set-up adiabat (per model, anchored at the column's 1-bar T), the RHS adiabatic
 * compression term (per model), the moist convective adjustment (ConvectiveAdjustment.h: h_L = H(T) + g z - sum L q) and
 * the buoyancy sink of the turbulence (Turbulence.h: Gam = g / cp(T)). NOT in: the dry CONV_ADJ (default off; warns),
 * the latent-heat divisors (cp_mix rho) of the saturation adjustment, precipitation and chemistry.
 * A model that has not ported the set-up returns false from Planet::cp_lapse_supported(); the knobs are then refused
 * with a message and stay 0 (2026-11-24: ATJUP; ported: ATSAT, ATURAN, ATNEPT). MEASURED 11-24 (from scratch to 1000,
 * <repo>/satchk/cph2_2026-11-24): Saturn CP_H2 1 -1..-8 K vs A&W Fig. 3 at 1-20 bar (default -6..-19); Uranus CP_GAS 1 +
 * CP_H2 2 -4..+8 K vs Atreya 2020 Fig. 5 at 1-10 bar (default +5..+38) and the H2S deck on Fig. 6's; Neptune CP_GAS 1 +
 * CP_H2 1 within ~2 K of Voyager to 5 bar (as the default). Planet::cp_lapse_in_dynamics() false keeps cp_mix in the RHS term,
 * the MCA and the turbulence (ATSAT, design (A): its set-up lapse is the CONFIG gam = 0.7 K/km with R_ref 3.89 J/(g K),
 * cp_eff = g/gam = 14286, not g/cp_mix; its set-up uses cp_eff cp_h2(T)/cp_h2 instead, see InitVariables_Sat.cpp).
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace CpH2 {

    constexpr double R_GAS = 8.314462618, M_H2 = 2.01588e-3, TH_R = 87.57, TH_V = 6332.0;
    constexpr double T_LO = 5.0, T_HI = 1500.0, DT = 0.5;

    // cp/R of H2 at T: mode 1 normal (3:1 frozen), 2 equilibrium
    inline double cpR_exact(double T, int mode){
        auto moments = [T](int J0, double wgt, double& Z, double& E1, double& E2){
            for(int J = J0; J < 60; J += 2){
                const double e = TH_R * J * (J + 1.0), b = wgt * (2.0 * J + 1.0) * std::exp(-e / T);
                Z += b; E1 += b * e; E2 += b * e * e;
            }
        };
        auto crot = [T](double Z, double E1, double E2){ const double m1 = E1 / Z, m2 = E2 / Z; return (m2 - m1 * m1) / (T * T); };
        double c;
        if(mode == 1){
            double Zp = 0, Ep = 0, Fp = 0, Zo = 0, Eo = 0, Fo = 0;
            moments(0, 1.0, Zp, Ep, Fp); moments(1, 1.0, Zo, Eo, Fo);
            c = 0.25 * crot(Zp, Ep, Fp) + 0.75 * crot(Zo, Eo, Fo);
        } else {
            double Z = 0, E1 = 0, E2 = 0;
            moments(0, 1.0, Z, E1, E2); moments(1, 3.0, Z, E1, E2);
            c = crot(Z, E1, E2);
        }
        const double x = TH_V / T, cv = (x < 700.0) ? x * x * std::exp(x) / ((std::exp(x) - 1.0) * (std::exp(x) - 1.0)) : 0.0;
        return 2.5 + c + cv;
    }

    struct Table { std::vector<double> cp, H; };     // J/(kg K) and its integral from T_LO (J/kg)
    inline const Table& table(int mode){
        static const Table t1 = [](){ Table t; const int n = int((T_HI - T_LO) / DT) + 1; t.cp.resize(n); t.H.resize(n);
            for(int i = 0; i < n; i++) t.cp[i] = cpR_exact(T_LO + i * DT, 1) * R_GAS / M_H2;
            t.H[0] = t.cp[0] * T_LO; for(int i = 1; i < n; i++) t.H[i] = t.H[i-1] + 0.5 * DT * (t.cp[i] + t.cp[i-1]); return t; }();
        static const Table t2 = [](){ Table t; const int n = int((T_HI - T_LO) / DT) + 1; t.cp.resize(n); t.H.resize(n);
            for(int i = 0; i < n; i++) t.cp[i] = cpR_exact(T_LO + i * DT, 2) * R_GAS / M_H2;
            t.H[0] = t.cp[0] * T_LO; for(int i = 1; i < n; i++) t.H[i] = t.H[i-1] + 0.5 * DT * (t.cp[i] + t.cp[i-1]); return t; }();
        return mode == 1 ? t1 : t2;
    }
    // H2 cp [J/(kg K)] and enthalpy [J/kg] at T; mode 0 = the model's constant
    inline double cp_h2(double T, int mode, double cp_const){
        if(mode <= 0) return cp_const;
        const Table& t = table(mode); const int n = (int)t.cp.size();
        const double x = (std::min(std::max(T, T_LO), T_HI) - T_LO) / DT; const int i = std::min(n - 2, (int)x); const double f = x - i;
        return t.cp[i] + f * (t.cp[i+1] - t.cp[i]);
    }
    inline double enth_h2(double T, int mode, double cp_const){
        if(mode <= 0) return cp_const * T;
        const Table& t = table(mode); const int n = (int)t.cp.size();
        if(T <= T_LO) return t.cp[0] * T;
        if(T >= T_HI) return t.H[n-1] + t.cp[n-1] * (T - T_HI);
        const double x = (T - T_LO) / DT; const int i = std::min(n - 2, (int)x); const double f = x - i;
        return t.H[i] + DT * f * (t.cp[i] + 0.5 * f * (t.cp[i+1] - t.cp[i]));   // exact for linear cp in the step
    }

    // unset = the model's default (Planet::cp_h2_default() / cp_gas_default(), 2026-11-24, user: "make the cp(T) knobs
    // default on all four": ATSAT 1 / 0, ATURAN 2 / 1, ATNEPT 1 / 1, ATJUP 0 / 0 -- its set-up is not ported); =0 restores.
    template<class Planet> int env_knob(const char* name, int dflt){
        const char* e = std::getenv((std::string(Planet::planet_tag()) + "_" + name).c_str()); return e ? std::atoi(e) : dflt; }
    template<class Planet> int mode(){
        static const int v = [](){ const int k = env_knob<Planet>("CP_H2", Planet::cp_h2_default());
            if(k != 0 && !Planet::cp_lapse_supported()){
                std::printf("      %s: %s_CP_H2=%d REFUSED -- this model has not ported the cp(T) set-up adiabat (CpH2.h); stays 0\n",
                            Planet::planet_tag(), Planet::planet_tag(), k); return 0; }
            return (k == 1 || k == 2) ? k : 0; }();
        return v; }
    template<class Planet> bool gas(){
        static const bool v = [](){ const int k = env_knob<Planet>("CP_GAS", Planet::cp_gas_default() ? 1 : 0);
            if(k != 0 && !Planet::cp_lapse_supported()){
                std::printf("      %s: %s_CP_GAS=%d REFUSED -- this model has not ported the cp(T) set-up adiabat (CpH2.h); stays 0\n",
                            Planet::planet_tag(), Planet::planet_tag(), k); return false; }
            return k != 0; }();
        return v; }
    template<class Planet> bool on(){ return mode<Planet>() > 0 || gas<Planet>(); }
    // in the dynamics (RHS compression term, MCA, turbulence) only where the model says so: ATSAT keeps cp_mix there
    // (design (A), 2026-11-24: its set-up lapse is the config gam, not g/cp_mix) and uses cp(T) in the set-up alone.
    template<class Planet> bool dyn_on(){ return on<Planet>() && Planet::cp_lapse_in_dynamics(); }

    // the model's cp for lapse rates [J/(kg K)] and its enthalpy [J/kg] at T [K]. A struct template so that each model
    // can make it a friend (template<class M> friend struct CpH2::Lapse;) -- it reads the private gas constants.
    template<class M> struct Lapse {
        static double cp(const M& m, double T){
            if(!dyn_on<M>()) return m.cp_mix;
            const double ch = cp_h2(T, mode<M>(), m.cp_h2);
            if(gas<M>()) return (m.r_h2 * ch + m.r_he * m.cp_he + m.r_ch4 * m.cp_ch4) / (m.r_h2 + m.r_he + m.r_ch4);
            return m.cp_mix + (m.r_h2 / m.r_mix) * (ch - m.cp_h2);
        }
        static double enth(const M& m, double T){
            if(!dyn_on<M>()) return m.cp_mix * T;
            const double hh = enth_h2(T, mode<M>(), m.cp_h2);
            if(gas<M>()) return (m.r_h2 * hh + (m.r_he * m.cp_he + m.r_ch4 * m.cp_ch4) * T) / (m.r_h2 + m.r_he + m.r_ch4);
            return (m.cp_mix - (m.r_h2 / m.r_mix) * m.cp_h2) * T + (m.r_h2 / m.r_mix) * hh;
        }
        static double cp_mix(const M& m){ return m.cp_mix; }
    };
    template<class M> double cp_lapse(const M& m, double T){ return Lapse<M>::cp(m, T); }
    template<class M> double enth(const M& m, double T){ return Lapse<M>::enth(m, T); }
    // T with enth(T) = h (Newton from T0; enth rises monotonically)
    template<class M> double enth_inv(const M& m, double h, double T0){
        if(!dyn_on<M>()) return h / Lapse<M>::cp_mix(m);
        double T = (T0 > 0.0) ? T0 : h / Lapse<M>::cp_mix(m);
        for(int it = 0; it < 30; it++){
            const double f = enth(m, T) - h, c = cp_lapse(m, T); const double dT = f / c;
            T -= dT; if(T < 1.0) T = 1.0;
            if(std::fabs(dT) < 1.0e-10 * T) break;
        }
        return T;
    }
}

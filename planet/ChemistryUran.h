/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone chemistry class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * Header-only: all method bodies are inline.
*/

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Array.h"
#include "Array_1D.h"
#include "cUranusModel.h"
#include "FluxLimiter.h"

class ChemistryUran {
    friend class cUranusModel;

public:
    explicit ChemistryUran(cUranusModel& model) : m(model) {}
    ~ChemistryUran() = default;
    ChemistryUran(const ChemistryUran&) = delete;
    ChemistryUran& operator=(const ChemistryUran&) = delete;

    // -----------------------------------------------------------------------
    void ChemMassRateUran()
    {
        using namespace std;
        cout << endl << "      ATURAN: ChemMassRateUran" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;
        // ATURAN_BC_PHI_PERIODIC (ATPhys.h): with 1 the seam meridian k = 0 is computed like any other.
        const int k_first = ATPhys::phi_first(ATPhys::phi_periodic<cUranusModel>());

        // reaction constants are the same for all cells
        const double keq = m.m_nh4sh / (m.m_nh3 * m.m_h2s);
        const double A   = 1500.0;
        const double B   = -0.3;
        const double T_d = 3020.0;

        // ATURAN_NH4SH_EQ=1: the rate law contributes nothing (w_* = 0); the diffusion term stays.
        const bool eq_on = nh4sh_eq();
        // ATURAN_CHEM_DIFFLUX_PLUS (2026-11-23, row GIANTS-CHEM-DIFFLUX; DEFAULT 0 = the minus below, unchanged):
        // the sign of the diffusive term in massflux_*. ATSAT's ATSAT_CHEM_DIFFLUX_PLUS (SAT-CHEM-DIVERGE item 4,
        // ChemistrySat.h): Fick's tendency is +D lap(c), so massflux = w MINUS difflux anti-diffuses. Belongs
        // with ATURAN_CHEM_DIFFLUX_LAPLACIAN=1 (DiffMassFluxUran); either alone is not the physical pair.
        // dsign is +1.0 or -1.0, so the default path is the subtraction it always was, exactly.
        static const double dsign = [](){
            const char* e = getenv("ATURAN_CHEM_DIFFLUX_PLUS"); return (e && atoi(e) != 0) ? 1.0 : -1.0; }();

        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = k_first; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    double t_u = m.t.x[i][j][k] * m.t_ref;
                    double kf  = react_rate_const(t_u, T_d, A, B);
                    double kb  = kf / keq;

                    if(eq_on){

                        // ATURAN_NH4SH_EQ: EquilibriumNH4SH forms / decomposes the solid; no rate law.

                        m.w_nh3.x[i][j][k]   = 0.0;

                        m.w_h2s.x[i][j][k]   = 0.0;

                        m.w_nh4sh.x[i][j][k] = 0.0;

                    }

                    else if((t_u <= m.t_0_nh4sh) && (t_u >= m.t_00_nh4sh)){
                        const double c_nh3   = m.nh3.x[i][j][k]   / m.m_nh3;
                        const double c_h2s   = m.h2s.x[i][j][k]   / m.m_h2s;
                        const double c_nh4sh = m.nh4sh.x[i][j][k] / m.m_nh4sh;

                        const double R_diff  = kf * c_nh3 * c_h2s - kb * c_nh4sh;

                        m.w_nh3.x[i][j][k]   = -m.m_nh3   * R_diff;
                        m.w_h2s.x[i][j][k]   = -m.m_h2s   * R_diff;
                        m.w_nh4sh.x[i][j][k] =  m.m_nh4sh * R_diff;
                    }

                    m.massflux_h2s.x[i][j][k]   = m.w_h2s.x[i][j][k]   + dsign * m.difflux_h2s.x[i][j][k];
                    m.massflux_nh3.x[i][j][k]   = m.w_nh3.x[i][j][k]   + dsign * m.difflux_nh3.x[i][j][k];
                    m.massflux_nh4sh.x[i][j][k] = m.w_nh4sh.x[i][j][k] + dsign * m.difflux_nh4sh.x[i][j][k];
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ChemMassRateUran\n", elapsed.count() * 1e-9);

        cout << "      ATURAN: ChemMassRateUran ended" << endl;
        return;
    }

    // ===================================================================================================
    // ATURAN_NH4SH_EQ (2026-11-22, port of ATSAT_NH4SH_EQ / ATJUP_NH4SH_EQ; default 0 = the rate law, output unchanged):
    // NH4SH as an EQUILIBRIUM ADJUSTMENT, replacing the rate law of ChemMassRateUran (whose w_* are then zero).
    //
    // WHY. The rate law has no source for its numbers (A / B / T_d = 1500 / -0.3 / 3020, keq = a ratio of molar
    // masses, not an equilibrium constant, no heat released). On Uranus it makes only ~0.001 g/m3 NH4SH, so the
    // set-up NH3 survives as a 10-20 g/m3 NH3 cloud, although H2S outnumbers NH3 (r_h2s / m_h2s : r_nh3 / m_nh3
    // = 1.5) and in equilibrium the NH4SH deck should take ALL the NH3 (de Pater & Lissauer: no NH3 cloud on the
    // ice giants; the surplus H2S forms the H2S ice cloud above).
    //
    // WHAT. NH3(g) + H2S(g) <-> NH4SH(s) is in equilibrium where p_NH3 p_H2S = K_p(T),
    //     log10 K_p [atm^2] = 14.82 - 4705 / T          (Lewis 1969; Weidenschilling & Lewis 1973),
    // with p_x = (c_x / m_x) R* T [Pa], c_x the species density [kg/m3]. Per fluid cell x [kmol/m3] of NH4SH
    // forms (x > 0) or decomposes (x < 0, at most all of it), heating the cell by L x m_nh4sh / (cp_mix rho),
    // L = 4705 ln(10) R* / m_nh4sh = 1.762e6 J/kg (the van't Hoff enthalpy of the same K_p). T and K_p are solved
    // together (bisection). N, S and mass conserved by construction. No window. rho = p_stat / (R_mix T), the
    // local density the saturation adjustment reads (and DiffMassFluxUran writes into rho_mix). Called after
    // the four saturation adjustments (set-up: parcel sweep; loop: in place); NH4SH joins ATURAN_RK_SCALAR_SYNC
    // (H2S is synced already); the MCA carries the enthalpy in NH4SH's tracer L.
    // ===================================================================================================
    static double nh4sh_formation_heat(){ return 4705.0 * std::log(10.0) * 8314.46 / 51.1114; }

    static bool nh4sh_eq(){
        static const bool v = [](){ const char* e = getenv("ATURAN_NH4SH_EQ"); return !e || atoi(e) != 0; }();   // default 1 since 2026-11-22
        return v;
    }

    // One cell to equilibrium; returns x [kmol/m3] formed (> 0) or decomposed (< 0) and applies it.
    double equilibrateNH4SHCell(int i, int j, int k, bool& all_decomposed, double& sat_entry, double& dT)
    {
        constexpr double R_univ = 8314.46, atm = 101325.0;
        static const double L_f = nh4sh_formation_heat();
        auto Kp = [&](double T){ return std::pow(10.0, 14.82 - 4705.0 / T) * atm * atm; };   // Pa^2
        all_decomposed = false; sat_entry = 0.0; dT = 0.0;
        const double T0 = m.t.x[i][j][k] * m.t_ref;
        if (!(T0 > 0.0)) return 0.0;
        double rho = m.p_stat.x[i][j][k] * 1e5 / (m.R_mix * T0);
        if (!(rho > 0.0) || !std::isfinite(rho)) rho = m.r_mix;
        const double a = std::max(0.0, m.nh3.x[i][j][k])   / m.m_nh3;          // kmol/m3
        const double b = std::max(0.0, m.h2s.x[i][j][k])   / m.m_h2s;
        const double s = std::max(0.0, m.nh4sh.x[i][j][k]) / m.m_nh4sh;
        const double heat = L_f * m.m_nh4sh / (m.cp_mix * rho);                // K per kmol/m3 formed
        auto f = [&](double x){
            const double T = T0 + heat * x, RT = R_univ * T;
            return (a - x) * (b - x) * RT * RT - Kp(T); };
        const double f0 = f(0.0);
        double lo, hi;
        if (f0 > 0.0)         { lo = 0.0; hi = std::min(a, b); sat_entry = f0 / Kp(T0); }
        else if (s > 0.0) {
            if (f(-s) <= 0.0) { lo = hi = -s; all_decomposed = true; }
            else              { lo = -s; hi = 0.0; }
        } else return 0.0;
        if (hi > lo) for (int it = 0; it < 80; it++) {
            const double xm = 0.5 * (lo + hi);
            if (f(xm) > 0.0) lo = xm; else hi = xm;
            if (hi - lo <= 1.0e-15 * std::max(a, std::max(b, s))) break;
        }
        const double x = (f0 > 0.0) ? lo : hi;                                 // never past the root
        if (x == 0.0) return 0.0;
        m.nh3.x[i][j][k]   = std::max(0.0, m.nh3.x[i][j][k]   - m.m_nh3   * x);
        m.h2s.x[i][j][k]   = std::max(0.0, m.h2s.x[i][j][k]   - m.m_h2s   * x);
        m.nh4sh.x[i][j][k] = std::max(0.0, m.nh4sh.x[i][j][k] + m.m_nh4sh * x);
        dT = heat * x;
        m.t.x[i][j][k] = (T0 + dT) / m.t_ref;
        return x;
    }

    // ATURAN_NH4SH_EQ_SWEEP=1 (2026-11-22, row GIANTS-NH4SH-EQ; default 0 = the sweep below unchanged): ECCM-style
    // start of the set-up parcel sweep. The ice-giant domains lie wholly ABOVE the NH4SH cloud base (Atreya et al.:
    // 30-40 bar; Neptune's floor 31.5, Uranus's 17.3 bar), so the NH4SH that the deep NH3 / H2S form on the way up
    // belongs to layers below the domain. With 0 the first fluid cell of each column equilibrates its own full
    // deep abundances IN PLACE, keeping all that NH4SH and its heat (Uranus it 8: 130 / 99 g/m3 in the two bottom
    // rows, max dT 28.7 K; it then relaxes, -16 % column N by it 700). With 1 that cell is the cloud base of the
    // parcel instead: the parcel is brought to equilibrium at the cell's set-up T (closed form, no heat), the cell
    // keeps the equilibrium vapour, NO NH4SH is deposited, and every cell above gets only what forms between it
    // and the level below -- the Weidenschilling & Lewis (1973) / Atreya & Wong (2005) construction.
    static bool nh4sh_eq_sweep(){
        static const bool v = [](){ const char* e = getenv("ATURAN_NH4SH_EQ_SWEEP"); return !e || atoi(e) != 0; }();   // default 1 since 2026-11-22
        return v;
    }

    // Equilibrium of a parcel at FIXED T: x [kmol/m3] with (a-x)(b-x) (R* T)^2 = K_p(T), 0 <= x <= min(a, b).
    // Smaller root of x^2 - (a+b) x + ab - K = 0, K = K_p / (R* T)^2, in the cancellation-free form.
    static double parcelNH4SHFixedT(double a, double b, double T){
        constexpr double R_univ = 8314.46, atm = 101325.0;
        if (!(T > 0.0) || !(a > 0.0) || !(b > 0.0)) return 0.0;
        const double RT = R_univ * T;
        const double K = std::pow(10.0, 14.82 - 4705.0 / T) * atm * atm / (RT * RT);
        if (a * b <= K) return 0.0;
        const double d = a - b;
        return 2.0 * (a * b - K) / ((a + b) + std::sqrt(d * d + 4.0 * K));
    }

    // where == "init": the set-up's column construction (Weidenschilling & Lewis 1973): each column is swept
    // bottom-up with a rising parcel; a cell takes NH3 and H2S as the parcel's remaining mass fractions (capped
    // by its own set-up value), is equilibrated, and passes what is left upward. Not conservative by design (it
    // builds a state). Any other where: every fluid cell in place, conservative. Columns are independent, so
    // the result does not depend on the thread count.
    void EquilibriumNH4SH(const char* where)
    {
        using namespace std;
        auto begin = chrono::high_resolution_clock::now();
        const bool sweep = std::string(where) == "init";
        const bool eccm = sweep && nh4sh_eq_sweep();                               // NH4SH_EQ_SWEEP, see above
        long n_form = 0, n_decomp = 0, n_all = 0, n_base = 0;
        double dm_form = 0.0, dm_decomp = 0.0, dT_max = 0.0, s_max = 0.0, dm_base = 0.0;
        #pragma omp parallel for collapse(2) schedule(dynamic, 8) \
            reduction(+:n_form, n_decomp, n_all, dm_form, dm_decomp, n_base, dm_base) reduction(max:dT_max, s_max)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                double q_nh3 = 1.0e300, q_h2s = 1.0e300;                            // the parcel, mass fractions
                bool base_ = eccm;                                                    // NH4SH_EQ_SWEEP: first fluid cell = parcel base
                for (int i = 0; i < m.im; i++) {
                    if (m.SeaMount.x[i][j][k] == 1.0) continue;
                    const double T_ = m.t.x[i][j][k] * m.t_ref;
                    double rho = (T_ > 0.0) ? m.p_stat.x[i][j][k] * 1e5 / (m.R_mix * T_) : 0.0;
                    if (!(rho > 0.0) || !isfinite(rho)) rho = m.r_mix;
                    if (sweep) {
                        m.nh3.x[i][j][k] = min(m.nh3.x[i][j][k], q_nh3 * rho);
                        m.h2s.x[i][j][k] = min(m.h2s.x[i][j][k], q_h2s * rho);
                    }
                    if (base_) {
                        // NH4SH_EQ_SWEEP: the parcel's NH4SH below this level is not in the domain -- equilibrate at the
                        // set-up T, keep the vapour, deposit nothing, release no heat.
                        base_ = false;
                        const double a_ = max(0.0, m.nh3.x[i][j][k]) / m.m_nh3, b_ = max(0.0, m.h2s.x[i][j][k]) / m.m_h2s;
                        const double xb = parcelNH4SHFixedT(a_, b_, T_);
                        m.nh3.x[i][j][k] = max(0.0, m.nh3.x[i][j][k] - m.m_nh3 * xb);
                        m.h2s.x[i][j][k] = max(0.0, m.h2s.x[i][j][k] - m.m_h2s * xb);
                        q_nh3 = m.nh3.x[i][j][k] / rho; q_h2s = m.h2s.x[i][j][k] / rho;
                        if (xb > 0.0) { n_base++; dm_base += m.m_nh4sh * xb; }
                        continue;
                    }
                    bool all_ = false; double sat_ = 0.0, dT_ = 0.0;
                    const double x = equilibrateNH4SHCell(i, j, k, all_, sat_, dT_);
                    if (sweep) { q_nh3 = m.nh3.x[i][j][k] / rho; q_h2s = m.h2s.x[i][j][k] / rho; }
                    s_max = max(s_max, sat_); dT_max = max(dT_max, fabs(dT_));
                    if (all_) n_all++;
                    if (x > 0.0)      { n_form++;   dm_form   += m.m_nh4sh * x; }
                    else if (x < 0.0) { n_decomp++; dm_decomp -= m.m_nh4sh * x; }
                }
            }
        }
        auto end = chrono::high_resolution_clock::now();
        printf("      ATURAN: EquilibriumNH4SH (%s%s, Lewis K_p, L %.4e J/kg): formed in %ld cells (sum %.4e kg/m3),"
               " decomposed in %ld (sum %.4e kg/m3, %ld completely), max p_NH3 p_H2S / K_p - 1 at entry %.3e,"
               " max dT %.4f K; %.3f s\n", where, sweep ? ", parcel sweep" : "",
               nh4sh_formation_heat(), n_form, dm_form, n_decomp, dm_decomp, n_all, s_max, dT_max,
               chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1e-9);
        if (eccm)
            printf("      ATURAN: EquilibriumNH4SH NH4SH_EQ_SWEEP: parcel base in the first fluid cell of %ld columns, NH4SH below"
                   " the domain (not deposited) sum %.4e kg/m3\n", n_base, dm_base);
    }

    // -----------------------------------------------------------------------
    // The TVD limiter is the SHARED FluxLimiter<Planet> template — ATURAN's hand-written copy
    // was the same code as ATSAT's, ATJUP's and ATNEPT's, differing only where each model answers
    // two questions the shared version now asks through hooks: rm = m.rad.z[i] against
    // m.metricRadius(m.rad.z[i]), and a local `constexpr sinthe_min = 0.4` against
    // ATPhys::polar_divisor_floor<Planet>(). Both reproduce ATURAN's previous values exactly by
    // default, which is why swapping it in is byte-identical. See FluxLimiter.h for why this one
    // routine could be shared while the rest of the chemistry cannot.
    void FluxLimiterNH4SH(){ FluxLimiter<cUranusModel>(m).nh4sh(); }

    // -----------------------------------------------------------------------
    void DiffMassFluxUran()
    {
        using namespace std;
        cout << endl << "      ATURAN: DiffMassFluxUran" << endl;

        // ATURAN_CHEM_DIFFLUX_LAPLACIAN (2026-11-23, row GIANTS-CHEM-DIFFLUX; DEFAULT 0 = the two-pass scalar below,
        // unchanged): ATSAT's ATSAT_CHEM_DIFFLUX_LAPLACIAN (SAT-CHEM-DIVERGE item 4, default there since e7f77de).
        // Pass 2 below is NOT a divergence: it sums the components of d(j)/dx as scalars, with |.| on the
        // theta part and no metric factors, and j_* carries the thermal-diffusion part jT ~ DT/T dT/dr, which
        // is nonzero where the species is ZERO -- so difflux_* creates NH3 and H2S out of nothing wherever
        // dT/dr varies. MEASURED (ATURAN/satchk/nh3ice_2026-11-23, cold-region CWB, 40 it): above 1.1 bar
        // +41 g/m2 NH3 and +26.6 H2S from 'chemistry source' = -difflux_*, all of it deposited as NH3 ice at
        // ~55 K (the ~0.03 g/m3 NH3 ice at 0.1-1 bar, growing linearly, on Uranus and Neptune).
        // =1: difflux_x = D_x * laplacian_spherical(c_x), the diffusive tendency (zero where c is uniform);
        // pair it with ATURAN_CHEM_DIFFLUX_PLUS=1. j_* and thermalmassflux are computed as before either way.
        static const int difflux_laplacian = [](){
            const char* e = getenv("ATURAN_CHEM_DIFFLUX_LAPLACIAN"); return e ? atoi(e) : 0; }();

        auto begin = std::chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;
        // ATURAN_BC_PHI_PERIODIC (ATPhys.h): with 1 the seam meridian k = 0 is computed like any other.
        const int k_first = ATPhys::phi_first(ATPhys::phi_periodic<cUranusModel>());

        // Diffusion coefficients for gas-phase species only.
        // NH4SH is a solid precipitate — no gas-phase diffusion.
        {
            const double nu_mix = m.mue_mix / m.r_mix;
            m.D_nh3  = nu_mix / m.sc_nh3;
            m.D_h2s  = nu_mix / m.sc_h2s;
            m.DT_nh3 = nu_mix / m.sc_nh3;
            m.DT_h2s = nu_mix / m.sc_h2s;
        }

        // Update local mixture density from ideal gas law before computing fluxes.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                for(int k = 0; k < km; k++){
                    m.rho_mix.x[i][j][k] = m.p_stat.x[i][j][k] * 1e5
                        / (m.R_mix * m.t.x[i][j][k] * m.t_ref);
                }
            }
        }

        // Pass 1: compute jT_* and j_* (ordinary + thermal diffusion fluxes).
        // All writes go to jT_* and j_* at [i][j][k]; reads are from
        // neighbouring cells of t, h2s, nh3, nh4sh — no write-write races.
        const double sinthe_min = m.sinthe_min();   // the model's polar metric floor

        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = k_first; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double rm       = m.rad.z[i];
                    const double sinthe   = std::max(sinthe_min, std::abs(sin(m.the.z[j])));
                    const double rmsinthe = rm * sinthe;

                    double dtdr     = 0.0,   dtdthe     = 0.0,   dtdphi     = 0.0;
                    double dh2sdr   = 0.0, dh2sdthe   = 0.0, dh2sdphi   = 0.0;
                    double dnh3dr   = 0.0, dnh3dthe   = 0.0, dnh3dphi   = 0.0;
                    double dnh4shdr = 0.0, dnh4shdthe = 0.0, dnh4shdphi = 0.0;

                    derivative_1_order(i, j, k, dtdr,     dtdthe,     dtdphi,     m.t);
                    derivative_1_order(i, j, k, dh2sdr,   dh2sdthe,   dh2sdphi,   m.h2s);
                    derivative_1_order(i, j, k, dnh3dr,   dnh3dthe,   dnh3dphi,   m.nh3);
                    derivative_1_order(i, j, k, dnh4shdr, dnh4shdthe, dnh4shdphi, m.nh4sh);

                    m.jT_nh3.x[i][j][k]   = m.DT_nh3   / m.t.x[i][j][k] * dtdr;
                    m.jT_h2s.x[i][j][k]   = m.DT_h2s   / m.t.x[i][j][k] * dtdr;
                    m.jT_nh4sh.x[i][j][k] = m.DT_nh4sh / m.t.x[i][j][k] * dtdr;

                    // Pole-symmetric sum of gradient components: dXdr and dXdphi are
                    // pole-symmetric for symmetric inputs, dXdthe is pole-antisymmetric.
                    // std::abs() on the theta term symmetrizes the resulting scalar.
                    const double dnh3   = dnh3dr   + std::abs(dnh3dthe)/rm   + dnh3dphi/rmsinthe;
                    const double dh2s   = dh2sdr   + std::abs(dh2sdthe)/rm   + dh2sdphi/rmsinthe;
                    const double dnh4sh = dnh4shdr + std::abs(dnh4shdthe)/rm + dnh4shdphi/rmsinthe;

                    m.j_nh3.x[i][j][k]   = m.c_mix/m.r_mix * (m.m_nh3   * m.D_nh3   * dnh3
                                          - m.jT_nh3.x[i][j][k]);
                    m.j_h2s.x[i][j][k]   = m.c_mix/m.r_mix * (m.m_h2s   * m.D_h2s   * dh2s
                                          - m.jT_h2s.x[i][j][k]);
                    m.j_nh4sh.x[i][j][k] = m.c_mix/m.r_mix * (m.m_nh4sh * m.D_nh4sh * dnh4sh
                                          - m.jT_nh4sh.x[i][j][k]);
                }
            }
        }

        // With the periodic seam k = km-1 is the copy of k = 0. bcPhi does that for the fields of
        // its list; j_nh4sh and jT_nh4sh are not on it, so they are copied here.
        if(ATPhys::phi_periodic<cUranusModel>()){
            #pragma omp parallel for collapse(2) schedule(static)
            for(int i = 0; i < im; i++){
                for(int j = 0; j < jm; j++){
                    m.j_nh4sh.x[i][j][km-1]  = m.j_nh4sh.x[i][j][0];
                    m.jT_nh4sh.x[i][j][km-1] = m.jT_nh4sh.x[i][j][0];
                }
            }
        }

        // Pass 2: compute difflux_* (divergence of j_*) and thermalmassflux.
        // Reads j_* from neighbouring cells — all written by pass 1.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = k_first; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double rm       = m.rad.z[i];
                    const double sinthe   = std::max(sinthe_min, std::abs(sin(m.the.z[j])));
                    const double rmsinthe = rm * sinthe;

                    double dj_h2sdr   = 0.0, dj_h2sdthe   = 0.0, dj_h2sdphi   = 0.0;
                    double dj_nh3dr   = 0.0, dj_nh3dthe   = 0.0, dj_nh3dphi   = 0.0;
                    double dj_nh4shdr = 0.0, dj_nh4shdthe = 0.0, dj_nh4shdphi = 0.0;
                    double dtdr       = 0.0, dtdthe       = 0.0, dtdphi       = 0.0;

                    derivative_1_order(i, j, k, dj_h2sdr,   dj_h2sdthe,   dj_h2sdphi,   m.j_h2s);
                    derivative_1_order(i, j, k, dj_nh3dr,   dj_nh3dthe,   dj_nh3dphi,   m.j_nh3);
                    derivative_1_order(i, j, k, dj_nh4shdr, dj_nh4shdthe, dj_nh4shdphi, m.j_nh4sh);
                    derivative_1_order(i, j, k, dtdr,       dtdthe,       dtdphi,       m.t);

                    if(difflux_laplacian){
                        // A genuine D * grad^2(c), as ATSAT and ATJUP compute it. See the note at the top of this routine.
                        m.difflux_h2s.x[i][j][k]   = m.D_h2s   * laplacian_spherical(i, j, k, m.h2s);
                        m.difflux_nh3.x[i][j][k]   = m.D_nh3   * laplacian_spherical(i, j, k, m.nh3);
                        m.difflux_nh4sh.x[i][j][k] = m.D_nh4sh * laplacian_spherical(i, j, k, m.nh4sh);
                    }else{
                        m.difflux_h2s.x[i][j][k]   = dj_h2sdr   + std::abs(dj_h2sdthe)/rm   + dj_h2sdphi/rmsinthe;
                        m.difflux_nh3.x[i][j][k]   = dj_nh3dr   + std::abs(dj_nh3dthe)/rm   + dj_nh3dphi/rmsinthe;
                        m.difflux_nh4sh.x[i][j][k] = dj_nh4shdr + std::abs(dj_nh4shdthe)/rm + dj_nh4shdphi/rmsinthe;
                    }

                    m.thermalmassflux.x[i][j][k] =
                          (m.j_nh3.x[i][j][k]   * m.cp_nh3
                         + m.j_h2s.x[i][j][k]   * m.cp_h2s
                         + m.j_nh4sh.x[i][j][k] * m.cp_nh4sh)
                         * (dtdr + std::abs(dtdthe)/rm + dtdphi/rmsinthe)
                         + m.t.x[i][j][k] * (m.w_nh3.x[i][j][k]   * m.k_nh3
                                           + m.w_h2s.x[i][j][k]   * m.k_h2s
                                           + m.w_nh4sh.x[i][j][k] * m.k_nh4sh);
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for DiffMassFluxUran\n", elapsed.count() * 1e-9);

        cout << "      ATURAN: DiffMassFluxUran ended" << endl;
        return;
    }

    // -----------------------------------------------------------------------
    void ThermalPropertiesUran()
    {
        using namespace std;
        cout << endl << "      ATURAN: ThermalPropertiesUran" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        double r_mixture = 3.61;

        m.rho_cond_mix = m.rho_cond_h2 + m.rho_cond_he + m.rho_cond_h2s
            + m.rho_cond_nh3 + m.rho_cond_nh4sh + m.rho_cond_h2o + m.rho_cond_ch4;
        m.r_mix  = m.r_h2  + m.r_he  + m.r_h2s  + m.r_nh3  + m.r_nh4sh  + m.r_h2o  + m.r_ch4;
        m.c_mix  = m.c_h2  + m.c_he  + m.c_h2s  + m.c_nh3  + m.c_nh4sh  + m.c_h2o  + m.c_ch4;

        double M_mix = m.r_mix/m.c_mix;

        // MASS-FRACTION WEIGHTING, which is ATJUP's form (ChemistryJup.h:441) and was not ATURAN's.
        //
        // These three are per-unit-MASS properties, so a mixture value is sum(mass fraction * x),
        // i.e. sum(r_x * X_x) / r_mix. What stood here divided each species by its OWN MOLAR MASS
        // instead of by the mixture density: sum(r_x/m_x * X_x). Those weights are mole densities
        // and they sum to c_mix = 0.4934 kmol/m3, not to 1, so the result was neither a mass
        // average nor a molar one — and every species heavier than H2 was suppressed by the ratio
        // m_x/r_mix, which for CH4 is a factor of 9.6.
        //
        // THIS IS A DEFECT AND NOT A MODELLING CHOICE, on the routine's own evidence. R_mix below
        // is already normalised by the density sum, so the two lines disagreed with each other;
        // and for an ideal gas cp - cv = R, which makes the old cp_mix say
        //     gamma = cp/(cp - R_mix) = 6380.32/(6380.32 - 2447.84) = 1.6224
        // for an atmosphere that is 51.5 % H2 by mass — essentially the monatomic 5/3. A
        // H2-dominated atmosphere should be near 1.4; the corrected value gives 1.3923.
        //
        // ON URANUS THIS MOVES THE INITIAL PROFILE, because gam is not a config constant here:
        // cUranusModel.cpp:142 computes gam = g*1e3/cp_mix immediately before init_temperature,
        // and init_temperature builds T_bottom = T_top + gam*L_atm. Uranus shares Neptune's
        // mixture exactly (r_mix = 1.677), so cp_mix is the same number; only g = 8.69 and
        // L_atm = 360 km differ:
        //     cp_mix        6380.3192059532 -> 8683.24... J/(kg K)   +36.1 %
        //     gam           1.3620 -> 1.0008 K/km
        //     deep equator   551.0 ->  421.0 K   (T_top 60.65 K + gam*360 km)
        //
        // ATURAN IS OTHERWISE UNTOUCHED and still sits at its initial fork: none of the eleven
        // shared physics headers have reached it. This fix is applied here only because the
        // defect is in ATURAN's own chemistry file and is independent of that port.
        // ATSAT and ATNEPT carry the same defect and are fixed alongside.
        m.cp_mix  = (m.r_h2 * m.cp_h2 + m.r_h2s * m.cp_h2s
                  +  m.r_he * m.cp_he + m.r_nh3 * m.cp_nh3
                  +  m.r_nh4sh * m.cp_nh4sh + m.r_h2o * m.cp_h2o
                  +  m.r_ch4 * m.cp_ch4) / m.r_mix;
        m.mue_mix = (m.r_h2 * m.mue_h2 + m.r_h2s * m.mue_h2s
                  +  m.r_he * m.mue_he + m.r_nh3 * m.mue_nh3
                  +  m.r_nh4sh * m.mue_nh4sh + m.r_h2o * m.mue_h2o
                  +  m.r_ch4 * m.mue_ch4) / m.r_mix;
        m.k_mix   = (m.r_h2 * m.k_h2 + m.r_h2s * m.k_h2s
                  +  m.r_he * m.k_he + m.r_nh3 * m.k_nh3
                  +  m.r_nh4sh * m.k_nh4sh + m.r_h2o * m.k_h2o
                  +  m.r_ch4 * m.k_ch4) / m.r_mix;
        m.R_mix   = (m.r_h2 * m.R_h2 + m.r_he * m.R_he + m.r_h2o * m.R_h2o
            + m.r_h2s * m.R_h2s + m.r_nh3 * m.R_nh3 + m.r_nh4sh * m.R_nh4sh
            + m.r_ch4 * m.R_ch4)
            /(m.r_h2 + m.r_he + m.r_h2o + m.r_h2s + m.r_nh3 + m.r_nh4sh + m.r_ch4);

        cout.precision(10);
        cout.setf(ios::fixed);
        cout << endl
            << "     r_mixture[kg/m³] = " << r_mixture << endl << endl

            << "     rho_cond_h2[kg/m³] = "    << m.rho_cond_h2    << endl
            << "     rho_cond_he[kg/m³] = "    << m.rho_cond_he    << endl
            << "     rho_cond_nh3[kg/m³] = "   << m.rho_cond_nh3   << endl
            << "     rho_cond_h2s[kg/m³] = "   << m.rho_cond_h2s   << endl
            << "     rho_cond_nh4sh[kg/m³] = " << m.rho_cond_nh4sh << endl
            << "     rho_cond_h2o[kg/m³] = "   << m.rho_cond_h2o   << endl
            << "     rho_cond_ch4[kg/m³] = "   << m.rho_cond_ch4   << endl
            << "     rho_cond_mix[kg/m³] = "   << m.rho_cond_mix   << endl << endl

            << "     r_h2[kg/m³] = "     << m.r_h2    << endl
            << "     r_he[kg/m³] = "     << m.r_he    << endl
            << "     r_nh3[kg/m³] = "    << m.r_nh3   << endl
            << "     r_h2s[kg/m³] = "    << m.r_h2s   << endl
            << "     r_nh4sh[kg/m³] = "  << m.r_nh4sh << endl
            << "     r_h2o[kg/m³] = "    << m.r_h2o   << endl
            << "     r_ch4[kg/m³] = "    << m.r_ch4   << endl
            << "     r_mix[kg/m³] = "    << m.r_mix   << endl << endl

            << "     c_h2[kmol/m³] = "   << m.c_h2    << endl
            << "     c_he[kmol/m³] = "   << m.c_he    << endl
            << "     c_nh3[kmol/m³] = "  << m.c_nh3   << endl
            << "     c_h2s[kmol/m³] = "  << m.c_h2s   << endl
            << "     c_nh4sh[kmol/m³] = "<< m.c_nh4sh << endl
            << "     c_h2o[kmol/m³] = "  << m.c_h2o   << endl
            << "     c_ch4[kmol/m³] = "  << m.c_ch4   << endl
            << "     c_mix[kmol/m³] = "  << m.c_mix   << endl << endl

            << "     ep_h2[kg/m³] = "    << m.ep_h2   << endl
            << "     ep_he[kg/m³] = "    << m.ep_he   << endl
            << "     ep_nh3[kg/m³] = "   << m.ep_nh3  << endl
            << "     ep_h2s[kg/m³] = "   << m.ep_h2s  << endl
            << "     ep_nh4sh[kg/m³] = " << m.ep_nh4sh << endl
            << "     ep_ch4[kg/m³] = "   << m.ep_ch4  << endl
            << "     ep_h2o[kg/m³] = "   << m.ep_h2o  << endl << endl

            << "     X_h2[/] = "    << m.X_h2    << endl
            << "     X_he[/] = "    << m.X_he    << endl
            << "     X_nh3[/] = "   << m.X_nh3   << endl
            << "     X_h2s[/] = "   << m.X_h2s   << endl
            << "     X_nh4sh[/] = " << m.X_nh4sh << endl
            << "     X_ch4[/] = "   << m.X_ch4   << endl
            << "     X_h2o[/] = "   << m.X_h2o   << endl << endl

            << "     gam[/] = "         << m.gam    << endl
            << "     M_mix[kg/Kmol] = " << M_mix    << endl
            << "     r_mix[kg/m³] = "   << m.r_mix  << endl
            << "     c_mix[kmol/m³] = " << m.c_mix  << endl << endl

            << "     cp_mix[J/(kg*K)] = "  << m.cp_mix  << endl
            << "     mue_mix[Ns/m²] = "   << m.mue_mix << endl
            << "     k_mix[W/(m*K)] = "   << m.k_mix   << endl
            << "     R_mix[J/(kg*K)] = "  << m.R_mix   << endl << endl;

        // Initial rho_mix from ideal gas law using the reference state.
        // DiffMassFluxUran refreshes this every call with the current p_stat and t.
        #pragma omp parallel for schedule(static) collapse(2)
        for(int i = 0; i < m.im; i++){
            for(int j = 0; j < m.jm; j++){
                for(int k = 0; k < m.km; k++){
                    m.rho_mix.x[i][j][k] = m.p_stat.x[i][j][k] * 1e5
                        / (m.R_mix * m.t.x[i][j][k] * m.t_ref);
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThermalPropertiesUran\n", elapsed.count() * 1e-9);

        cout << "      ATURAN: ThermalPropertiesUran ended" << endl;
        return;
    }

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // Spherical Laplacian for ATURAN_CHEM_DIFFLUX_LAPLACIAN, ported from ATSAT's ChemistrySat::laplacian_spherical
    // (itself ATJUP's ChemistryJup::laplacian_spherical): rm through metricRadius(), the phi seam through
    // ATPhys::phi_wrap, the 1/(rm sinthe) divisor floored by ATPhys::polar_divisor_floor (0.4, as in the
    // shared FluxLimiter.h). coord_stretching is false here, so exp_rm folds to 1.
    double laplacian_spherical(int i, int j, int k, Array& c)
    {
        const double rm       = m.metricRadius(m.rad.z[i]);
        const double exp_rm   = m.coord_stretching ? 1.0 / (rm + 1.0) : 1.0;
        const double exp_2_rm = exp_rm * exp_rm;
        const double sinthe   = sin(m.the.z[j]);
        const double costhe   = cos(m.the.z[j]);
        const double rmsinthe = rm * std::max(sinthe, ATPhys::polar_divisor_floor<cUranusModel>());

        const double d2cdr2   = (c.x[i+1][j][k] - 2.0*c.x[i][j][k] + c.x[i-1][j][k]) / (m.dr   * m.dr) * exp_2_rm;
        const double d2cdthe2 = (c.x[i][j+1][k] - 2.0*c.x[i][j][k] + c.x[i][j-1][k]) / (m.dthe * m.dthe);
        const bool phi_per = ATPhys::phi_periodic<cUranusModel>();
        const int kw = ATPhys::phi_wrap(phi_per, k-1, m.km), ke = ATPhys::phi_wrap(phi_per, k+1, m.km);
        const double d2cdphi2 = (c.x[i][j][ke] - 2.0*c.x[i][j][k] + c.x[i][j][kw]) / (m.dphi * m.dphi);

        const double dcdr   = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0 * m.dr) * exp_rm;
        const double dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0 * m.dthe);

        return d2cdr2
             + 2.0 / rm * dcdr
             + d2cdthe2 / (rm * rm)
             + costhe / (rm * rmsinthe) * dcdthe
             + d2cdphi2 / (rmsinthe * rmsinthe);
    }

    void derivative_1_order(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        const double dr = m.dr, dthe = m.dthe, dphi = m.dphi;

        // The phi neighbours: k-1, k+1, k-2, k+2, taken round the ring at the seam when
        // ATURAN_BC_PHI_PERIODIC=1 (ATPhys.h).
        const bool phi_per = ATPhys::phi_periodic<cUranusModel>();
        const int kw  = ATPhys::phi_wrap(phi_per, k-1, m.km), ke  = ATPhys::phi_wrap(phi_per, k+1, m.km);
        const int kw2 = ATPhys::phi_wrap(phi_per, k-2, m.km), ke2 = ATPhys::phi_wrap(phi_per, k+2, m.km);

        m.SeaMount.x[i][j][k] = 0.0;                                    // no SeaMount introduced

        // radial gradients
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i+1][j][k] == 0.0))
            &&(m.SeaMount.x[i+2][j][k] == 0.0)){
            c.x[i][j][k] = c.x[i+3][j][k]
                - 3.0 * c.x[i+2][j][k] + 3.0 * c.x[i+1][j][k];
            dcdr = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i+1][j][k]
                - c.x[i+2][j][k])/(2.0 * dr);
        }
        if(((m.SeaMount.x[i-1][j][k] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i+1][j][k] == 0.0)){
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k])/(2.0 * dr);
        }
        if(((m.SeaMount.x[i-1][j][k] == 0.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i+1][j][k] == 0.0)){
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k])/(2.0 * dr);
        }

        // south gradients
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i][j+1][k] == 0.0))
            &&(m.SeaMount.x[i][j+2][k] == 0.0)){
            c.x[i][j][k] = m.c43 * c.x[i][j+1][k] - m.c13 * c.x[i][j+2][k];
            dcdthe = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i][j+1][k]
                - c.x[i][j+2][k])/(2.0 * dthe);
        }
        if(((m.SeaMount.x[i][j-1][k] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j+1][k] == 0.0)){
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k])/(2.0 * dthe);
        }

        // north gradients
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i][j-1][k] == 0.0))
            &&(m.SeaMount.x[i][j-2][k] == 0.0)){
            c.x[i][j][k] = m.c43 * c.x[i][j-1][k] - m.c13 * c.x[i][j-2][k];
            dcdthe = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i][j-1][k]
                - c.x[i][j-2][k])/(2.0 * dthe);
        }
        if(((m.SeaMount.x[i][j-1][k] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j+1][k] == 0.0)){
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k])/(2.0 * dthe);
        }
        if(((m.SeaMount.x[i][j-1][k] == 0.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j+1][k] == 0.0)){
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k])/(2.0 * dthe);
        }

        // east gradients
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i][j][ke] == 0.0))
            &&(m.SeaMount.x[i][j][ke2] == 0.0)){
            c.x[i][j][k] = m.c43 * c.x[i][j][ke] - m.c13 * c.x[i][j][ke2];
            dcdphi = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i][j][ke]
                - c.x[i][j][ke2])/(2.0 * dphi);
        }
        if(((m.SeaMount.x[i][j][kw] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j][ke] == 0.0)){
            dcdphi = (c.x[i][j][ke] - c.x[i][j][kw])/(2.0 * dphi);
        }

        // west gradients
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i][j][kw] == 0.0))
            &&(m.SeaMount.x[i][j][kw2] == 0.0)){
            c.x[i][j][k] = m.c43 * c.x[i][j][kw] - m.c13 * c.x[i][j][kw2];
            dcdphi = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i][j][kw]
                - c.x[i][j][kw2])/(2.0 * dphi);
        }
        if(((m.SeaMount.x[i][j][kw] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j][ke] == 0.0)){
            dcdphi = (c.x[i][j][ke] - c.x[i][j][kw])/(2.0 * dphi);
        }
        if(((m.SeaMount.x[i][j][kw] == 0.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j][ke] == 0.0)){
            dcdphi = (c.x[i][j][ke] - c.x[i][j][kw])/(2.0 * dphi);
        }
        return;
    }

    // -----------------------------------------------------------------------
    void derivative_1_order_boundary(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        const int im = m.im, jm = m.jm, km = m.km;
        const double dr = m.dr, dthe = m.dthe, dphi = m.dphi;

        if(i == 0)    dcdr   = (- 3.0 * c.x[0][j][k]    + 4.0 * c.x[1][j][k]
                            - c.x[2][j][k])/(2.0 * dr);
        if(i == im-1) dcdr   = (- 3.0 * c.x[im-1][j][k] + 4.0 * c.x[im-2][j][k]
                            - c.x[im-3][j][k])/(2.0 * dr);

        if(j == 0)    dcdthe = (- 3.0 * c.x[i][0][k]    + 4.0 * c.x[i][1][k]
                            - c.x[i][2][k])/(2.0 * dthe);
        if(j == jm-1) dcdthe = (- 3.0 * c.x[i][jm-1][k] + 4.0 * c.x[i][jm-2][k]
                            - c.x[i][jm-3][k])/(2.0 * dthe);

        if(k == 0)    dcdphi = (- 3.0 * c.x[i][j][0]    + 4.0 * c.x[i][j][1]
                            - c.x[i][j][2])/(2.0 * dphi);
        if(k == km-1) dcdphi = (- 3.0 * c.x[i][j][km-1] + 4.0 * c.x[i][j][km-2]
                            - c.x[i][j][km-3])/(2.0 * dphi);
        return;
    }

    // -----------------------------------------------------------------------
    double react_rate_const(double &T_K, const double &T_d,
        const double &A, const double &B)
    {
        return A * pow(T_K, B) * exp(-T_d/T_K);
    }

private:
    cUranusModel& m;

    // Superbee: most compressive TVD limiter — best for sharp cloud fronts.
    static double superbee_phi(double r){
        return std::max(0.0, std::max(std::min(2.0*r, 1.0), std::min(r, 2.0)));
    }

    // Van Leer: smooth, differentiable — good general-purpose alternative.
    static double van_leer_phi(double r){
        return (r + std::abs(r)) / (1.0 + std::abs(r));
    }
};

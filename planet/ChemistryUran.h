/*
 * Atmosphere General Circulation Modell (ATURAN)
 * Standalone chemistry class for the Uranus model.
 * Declared as friend of cUranusModel so it may access all private members
 * through the stored reference.
 *
 * Header-only: all method bodies are inline.
*/

#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>

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

        // reaction constants are the same for all cells
        const double keq = m.m_nh4sh / (m.m_nh3 * m.m_h2s);
        const double A   = 1500.0;
        const double B   = -0.3;
        const double T_d = 3020.0;

        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    double t_u = m.t.x[i][j][k] * m.t_ref;
                    double kf  = react_rate_const(t_u, T_d, A, B);
                    double kb  = kf / keq;

                    if((t_u <= m.t_0_nh4sh) && (t_u >= m.t_00_nh4sh)){
                        const double c_nh3   = m.nh3.x[i][j][k]   / m.m_nh3;
                        const double c_h2s   = m.h2s.x[i][j][k]   / m.m_h2s;
                        const double c_nh4sh = m.nh4sh.x[i][j][k] / m.m_nh4sh;

                        const double R_diff  = kf * c_nh3 * c_h2s - kb * c_nh4sh;

                        m.w_nh3.x[i][j][k]   = -m.m_nh3   * R_diff;
                        m.w_h2s.x[i][j][k]   = -m.m_h2s   * R_diff;
                        m.w_nh4sh.x[i][j][k] =  m.m_nh4sh * R_diff;
                    }

                    m.massflux_h2s.x[i][j][k]   = m.w_h2s.x[i][j][k]   - m.difflux_h2s.x[i][j][k];
                    m.massflux_nh3.x[i][j][k]   = m.w_nh3.x[i][j][k]   - m.difflux_nh3.x[i][j][k];
                    m.massflux_nh4sh.x[i][j][k] = m.w_nh4sh.x[i][j][k] - m.difflux_nh4sh.x[i][j][k];
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ChemMassRateUran\n", elapsed.count() * 1e-9);

        cout << "      ATURAN: ChemMassRateUran ended" << endl;
        return;
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

        auto begin = std::chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;

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
        for(int k = 1; k < km-1; k++){
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

        // Pass 2: compute difflux_* (divergence of j_*) and thermalmassflux.
        // Reads j_* from neighbouring cells — all written by pass 1.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = 1; k < km-1; k++){
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

                    m.difflux_h2s.x[i][j][k]   = dj_h2sdr   + std::abs(dj_h2sdthe)/rm   + dj_h2sdphi/rmsinthe;
                    m.difflux_nh3.x[i][j][k]   = dj_nh3dr   + std::abs(dj_nh3dthe)/rm   + dj_nh3dphi/rmsinthe;
                    m.difflux_nh4sh.x[i][j][k] = dj_nh4shdr + std::abs(dj_nh4shdthe)/rm + dj_nh4shdphi/rmsinthe;

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
    void derivative_1_order(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        const double dr = m.dr, dthe = m.dthe, dphi = m.dphi;

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
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i][j][k+1] == 0.0))
            &&(m.SeaMount.x[i][j][k+2] == 0.0)){
            c.x[i][j][k] = m.c43 * c.x[i][j][k+1] - m.c13 * c.x[i][j][k+2];
            dcdphi = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i][j][k+1]
                - c.x[i][j][k+2])/(2.0 * dphi);
        }
        if(((m.SeaMount.x[i][j][k-1] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j][k+1] == 0.0)){
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1])/(2.0 * dphi);
        }

        // west gradients
        if(((m.SeaMount.x[i][j][k] == 1.0)&&(m.SeaMount.x[i][j][k-1] == 0.0))
            &&(m.SeaMount.x[i][j][k-2] == 0.0)){
            c.x[i][j][k] = m.c43 * c.x[i][j][k-1] - m.c13 * c.x[i][j][k-2];
            dcdphi = (- 3.0 * c.x[i][j][k] + 4.0 * c.x[i][j][k-1]
                - c.x[i][j][k-2])/(2.0 * dphi);
        }
        if(((m.SeaMount.x[i][j][k-1] == 1.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j][k+1] == 0.0)){
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1])/(2.0 * dphi);
        }
        if(((m.SeaMount.x[i][j][k-1] == 0.0)&&(m.SeaMount.x[i][j][k] == 0.0))
            &&(m.SeaMount.x[i][j][k+1] == 0.0)){
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1])/(2.0 * dphi);
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

/*
 * SHARED PHYSICS — multi-layer grey-body radiation, one implementation for every planet.
 * MUST BE BYTE-IDENTICAL IN EVERY MODEL THAT USES IT; `make check-shared` verifies that against
 * planet/SHARED.md5.
 *
 * It knows nothing about which planet it runs on. Everything planet-specific arrives through the
 * Planet template parameter:
 *
 *     Planet::planet_tag()       "ATJUP" / "ATSAT" — the log prefix, and the prefix its knobs are
 *                                read under (ATJUP_CIA_STRENGTH, ATSAT_SOLAR, ...)
 *     m.surface_index(j,k)       the first fluid level of a column (the same accessor the
 *                                turbulence closure and the microphysics use)
 *     Planet::rad_F_int()        intrinsic heat flux                      [W/m2]
 *     Planet::rad_S_solar()      solar constant at this planet's orbit    [W/m2]
 *     Planet::rad_albedo_bond()  Bond albedo, so S*(1-A) is ABSORBED
 *     Planet::rad_x_H2()         H2 mole fraction
 *     Planet::rad_x_He()         He mole fraction
 *
 * WHICH NUMBERS LIVE WHERE, AND WHY THE LINE IS DRAWN THERE. The five above are MEASURED
 * PROPERTIES OF A PLANET: they have different values on Jupiter and Saturn because the planets
 * differ, and a model that got one wrong would be wrong about that planet. They belong in the
 * planet's own header, where they can be read next to the rest of its parameters.
 *
 * Everything below — C_cia, he_ratio, the four kappa values, tau_cloud_cap, opac_cal, k_sw — is a
 * GREY-OPACITY CALIBRATION, not a measurement. One set of numbers is in force on all four planets
 * today, and they were tuned on Jupiter: C_cia so that CIA alone puts Jupiter's thermal
 * photosphere near 0.5 bar, opac_cal so that the combined opacity brings it to 0.25-0.35 bar.
 * Both of those are now measured and printed every run, by the diagnostic below, and as of the
 * recalibration recorded there Jupiter hits both.
 *
 * THAT TUNING DOES NOT TRANSFER, and it is now measured on all four rather than suspected on one.
 * Every model run single-threaded to nm = 224 with <TAG>_RADIATION=1, all other knobs off. Values
 * are the LAST sample; dT/diter is the drift of the photosphere temperature per ITERATION over the
 * final 32, which is the column that says how much the rest of the row can be trusted:
 *
 *      planet      tau=1  T(tau=1)  T_eff(in)      OLR   OLR/B    OLR/in   dT/diter
 *      ATJUP      0.3214    117.29     124.66   11.463   1.068     0.837    -0.0084
 *      ATSAT      0.0670     63.38      94.12    1.234   1.349     0.277    +0.0012
 *      ATURAN     0.1713    136.45      59.04   20.734   1.055    30.093    +0.1603
 *      ATNEPT     0.0576    205.74      59.28  100.385   0.988   143.407    +0.3359
 *
 * THREE OF THESE ROWS MOVED WHEN A CONSTANT WAS CORRECTED, AND THE PREVIOUS ONES WERE BETTER-LOOKING
 * AND LESS TRUE. ATURAN read 7.591 here and ATNEPT 2.860. Both were produced by TWO ERRORS PARTLY
 * CANCELLING. mue_ch4 carried methane's viscosity in CENTIPOISE as if it were Pa*s and mue_h2o
 * carried liquid water's, so the mass-weighted mue_mix came out ~150x too large — liquid-water
 * viscosity in a hydrogen atmosphere. On the ice giants mue_mix sets the species diffusivities and
 * hence the diffusive-enthalpy sink in rhs_t, so that sink was ~150x overweighted, cold enough to
 * drive a front into the t_min floor and, in the mean, to hold the column DOWN against a heating
 * excess nothing else was opposing. Correcting the constant removes the overweight and unmasks the
 * excess: it did not create it. A number that looks better because two mistakes disagree is the
 * failure mode this table exists to expose, and it caught it here on its own rows.
 *
 * THE SLOPE COLUMN EXISTS BECAUSE THE FIRST VERSION OF THIS TABLE WAS MEASURED AT TWO ITERATIONS
 * AND WAS WRONG. Not wrong about what it measured — wrong to be read as a steady state. A row
 * without a slope beside it cannot tell a converged planet from one mid-excursion, and the two gas
 * giants barely moving is exactly why the short run looked trustworthy at the time.
 *
 * THE SPLIT IS NOW GAS GIANT AGAINST ICE GIANT, AND THE PREVIOUS VERSION OF THIS PARAGRAPH SAID IT
 * WAS NOT. That statement was correct about the numbers in front of it — with the bad viscosity
 * ATJUP drifted -0.0084 K/iter and ATURAN -0.0088, indistinguishable — and it is wrong now. With the
 * constants right the gas giants are quiescent (-0.0084, +0.0012) and BOTH ice giants diverge
 * (+0.1603, +0.3359), an order of magnitude clear of either. The line is where it first looked, but
 * it was not visible while an overweighted sink was pinning Uranus flat. Two readings of the same
 * split, both honest on their own evidence; this one rests on constants that are not known to be
 * wrong.
 *
 * ATNEPT'S EARLIER ROW IS WORTH KEEPING IN VIEW, as the clearest demonstration on this table of what
 * the photosphere diagnostic is for. It read 17.2287 bar / 116.10 K / 14.651, and none of that was
 * opacity: init_PressureStatic anchored a p_bottom taken at T_bottom to a base of (T/t_ref) rather
 * than (T/T_bottom), inflating Neptune's whole pressure field by (T_bottom/t_ref)^n = 607x and
 * putting the TOP of the domain at 15 bar. tau = 1 was reached in the first layer or two below the
 * ceiling, so the scheme reported the top of the grid and faithfully radiated the 116 K it found
 * there. A calibration read of that row would have been a search for an opacity error that did not
 * exist.
 *
 * READ THE LAST THREE COLUMNS AGAINST EACH OTHER, because that is what the photosphere temperature
 * was added to make possible. OLR/B is the SCHEME's excess over a blackbody at the level it calls
 * the photosphere; OLR/in is the whole error. Over 224 iterations and four planets the scheme term
 * stays between 0.99 and 1.35 while the total error now spans a factor of FIVE HUNDRED, 0.277 to
 * 143.407. That is the strongest form this argument has taken: the error moved by more than two
 * orders of magnitude on ATNEPT alone when a viscosity was corrected, and the scheme term moved from
 * 0.923 to 0.988. The two-stream sum is not the variable and never has been. The temperature of the
 * column it is handed is — and on Neptune it was once not even the temperature so much as WHICH
 * LEVEL sat at that pressure.
 *
 * THE SIGN IS STILL NOT COMMON, so "the grey scheme over-emits" remains false as a general statement,
 * but the margins are no longer comparable. Jupiter and Saturn UNDER-emit, photospheres 7.4 K and
 * 30.7 K too COLD, both quiescent. Uranus and Neptune OVER-emit by 77.4 K and 146.5 K and are
 * climbing. The discriminator is still where tau reaches 1: too opaque puts the photosphere high and
 * cold (Saturn, 0.067 bar), too transparent puts it deep and hot. All four sit between 0.058 and
 * 0.32 bar. Jupiter, the planet the opacity was tuned on, is the only one whose photosphere
 * temperature is close to right, and it is the only one that has never needed a constant corrected.
 *
 * WHICH ICE GIANT IS WORSE HAS NOW REVERSED TWICE, which is a reason to distrust any single reading
 * of this table. At two iterations Neptune looked catastrophic (14.651) and Uranus mild (2.508). Run
 * to 224 with the bad viscosity, Uranus was the worse (7.591, quiescent) and Neptune the better
 * (2.860, moving). With the viscosity corrected, Neptune is far the worse again at 143.407 against
 * Uranus's 30.093, and both diverge. Nothing about the radiation scheme changed across any of those
 * three readings.
 *
 * WHAT IS LEFT, stated so the next reader does not mistake the current rows for a calibration
 * result: an unopposed HEATING excess on both ice giants, previously masked by a sink that was
 * ~150x too strong. Until that is found, the opacity constants below cannot be judged against these
 * two rows at all — a photosphere 77 K and 146 K too warm says nothing about kappa.
 *
 * They stay here, as one shared calibration, precisely so that stays ONE question rather than four
 * independently drifting answers. The per-planet lever is the runtime knob
 * (<TAG>_CIA_STRENGTH, <TAG>_OPACITY_STRENGTH), not a second copy of the constant. The day any
 * planet is genuinely calibrated on its own measurements, its coefficient moves up into its model
 * header alongside F_int and stops being shared — and that move should be a commit that says so.
 *
 * WHAT THE TWO COPIES ACTUALLY DIFFERED BY, before this: RadiationJup.h and RadiationSat.h were
 * 293 and 276 lines. Of the real differences, five were the planetary constants above, one was
 * the column base (i_topography against a SeaMount scan — now the shared accessor), and one was
 * the density (a stored rho_mix against a locally formed one — already reconciled in its own
 * commit, since it changed a number and had to be measured rather than merged). The rest was the
 * pair of column diagnostics ATSAT grew and ATJUP never had, which are now in both.
 */
/*
 * Multi-layer grey-body radiation — DIAGNOSTIC SCAFFOLD.
 *
 * PROVENANCE. Only the ARCHITECTURE comes from ATOM's MultiLayerRadiation: per-column grey layers,
 * Stefan-Boltzmann emission sigma*T^4, layer emissivity eps = 1 - exp(-tau), a two-stream up/down
 * net-flux solve, and the flux divergence as a heating rate. The Earth physics of the original
 * (surface albedo feedback, ocean/land split, Bignami and Atwater-Ball laws, the surface-BC-
 * entangled tridiagonal Thomas solve) applies to no giant planet and is not here.
 *
 * VERTICAL ORIENTATION. i = surface_index(j,k) is the deep bottom — the intrinsic-heat boundary —
 * and i = im-1 is the top, radiating to space.
 *
 * SCOPE. The result goes to DIAGNOSTIC arrays only (radiation, epsilon, Q_rad); it does not touch
 * t or any rhs. Wiring the heating into the temperature equation is a separate step behind its own
 * second knob (<TAG>_RAD_COUPLING) in both models. In radiative equilibrium every interface net
 * flux is uniform (= F_int) and Q_rad -> 0; a column out of equilibrium shows the tendency toward
 * it, and that is the intended self-test.
 *
 * Gated by <TAG>_RADIATION (default 0 = off, bit-identical).
*/

#pragma once

#include "ATPhys.h"   // ATPhys::env_double / env_int — the per-planet knob prefix

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

class Array;

template<class Planet>
class Radiation {
public:

    explicit Radiation(Planet& model) : m(model) {}

    static constexpr double sigma = 5.670374419e-8;  // Stefan-Boltzmann [W/m2/K4]

    // --- Absorbed shortwave (solar) ---
    // Injecting only F_int leaves a column radiating more than entered it, so it can never be in
    // radiative balance — matching an observed OLR that way would be right for the wrong reason.
    //
    // Solar is a SEPARATE SHORTWAVE CHANNEL, not a downward thermal flux at the top. The thermal
    // opacity here is CIA-dominated (tau ~ P^2/T), which is the wrong absorber for sunlight —
    // sunlight is absorbed by CH4 bands and haze high up — and mixing the two would also corrupt
    // the net-flux diagnostic, netting solar against thermal emission in Fu - Fd.
    //
    // Diurnally averaged insolation for a fast rotator is S*(1-A)*cos(lat)/pi, whose area-weighted
    // mean is exactly S*(1-A)/4:
    //   mean = S(1-A)/pi * integral(cos^2) / integral(cos) = S(1-A)/pi * (pi/2)/2 = S(1-A)/4.
    // This is an annual mean. On a planet with small obliquity (Jupiter, 3.1 deg) that is nearly
    // exact; on one with large obliquity (Saturn, 26.7 deg) a real seasonal cycle is being averaged
    // away, and that approximation is not accounted for anywhere else.
    static constexpr double k_sw = 1.0;   // shortwave optical depth per bar, measured from the top

    // --- H2/He collision-induced absorption (CIA), grey parametrization ---
    // CIA is the dominant thermal-IR opacity on the giant planets. It needs two collision partners,
    // each with number density n ~ P/T, so the volume absorption goes as (P/T)^2; over a
    // hydrostatic mass path dm = dP/g this gives a LAYER optical depth
    //     dtau = C_cia * comp * (P/T) * dP / g
    // i.e. the characteristic tau ~ P^2/T of CIA. C_cia is a grey coefficient lumping the
    // frequency-integrated Borysow binary coefficient / k_B / mean molecular mass; frequency-
    // resolved tables are a later refinement. comp = x_H2^2 + he_ratio*x_H2*x_He weights the H2-H2
    // and H2-He collision-pair probabilities. See the calibration note in the banner above.
    static constexpr double P0_phys  = 1.0e5;   // physical pressure [Pa] at the p_stat = 1 level
    static constexpr double he_ratio = 0.6;     // H2-He / H2-H2 grey binary-coefficient ratio
    static constexpr double C_cia    = 3.5e-6;  // grey CIA coefficient (Jupiter photosphere ~0.5 bar)

    // --- CH4/NH3 gas bands + cloud/ice continuum, additive to the CIA tau ---
    // Each absorber adds dtau = kappa * q * dP/g, where q is the mass mixing ratio [kg/kg] and
    // dP/g the hydrostatic mass path [kg/m2]: CH4 (7.7 um band) and NH3 (rotational + v2) grey
    // band mass-absorption coefficients, and suspended cloud liquid and ice grey continuum after
    // Stephens-type values.
    //
    // Both models over-condense (as ATOM did: column condensate >> observed), which would make
    // every cloudy layer a blackbody and crush the OLR, so each layer's CLOUD tau is capped at
    // tau_cloud_cap.
    //
    // What is calibrated is the grey OPTICAL DEPTH kappa*q, never kappa alone — so these are only
    // as good as the species fields they multiply. The CH4 field in particular is large (q_ch4 ~
    // 0.246 kg/kg on Jupiter, ~100x the real Jovian value), which is why kappa_ch4 is
    // correspondingly small.
    //
    // opac_cal WAS 0.25, AND THAT VALUE WAS STALE BY EXACTLY THE FACTOR IT LOOKS LIKE. It was
    // tuned when the species fields were divided by the wrong density, so the mixing ratios it was
    // fitted against were far too small; once that was fixed the same coefficient put Jupiter's
    // photosphere at 0.031 bar instead of the 0.25-0.35 it was chosen for. Recalibrated on
    // Jupiter, 30 iterations from cold, by sweeping ATJUP_OPACITY_STRENGTH:
    //
    //     strength    opac_cal    photosphere    mean OLR   (against 13.695 W/m2 in)
    //       1.0        0.25         0.0307 bar     8.319
    //       0.20       0.05         0.1636         10.222
    //       0.12       0.03         0.2886         11.576
    //       0.10       0.025        0.3259         12.070      <-- adopted
    //       0.08       0.02         0.3555         12.665
    //       0.05       0.0125       0.4059         13.817
    //       0.0        0            0.5108         18.091      <-- CIA alone
    //
    // THE TWO TARGETS DO NOT COINCIDE, and the photosphere is the one being hit. Radiative balance
    // (OLR = 13.7 in) wants ~0.0125, which puts the photosphere at 0.406 bar, outside the range;
    // 0.025 puts the photosphere mid-range and leaves the column radiating 12% less than enters
    // it. That residual is not opac_cal's to fix — it is 30 iterations from cold with a model top
    // this file elsewhere notes is anomalously cold, and radiation does not yet feed T by default.
    //
    // The bottom row is the floor: with the gas and cloud opacity switched off entirely, CIA alone
    // puts the photosphere at 0.51 bar. That is C_cia's own calibration target, reproduced, and it
    // is why the fix belonged in opac_cal and not in C_cia.
    static constexpr double kappa_ch4     = 0.003;  // CH4 grey band mass opacity [m2/kg]
    static constexpr double kappa_nh3     = 1.5;    // NH3 grey band mass opacity [m2/kg]
    static constexpr double kappa_cloud   = 25.0;   // liquid-cloud grey mass opacity [m2/kg]
    static constexpr double kappa_ice     = 12.0;   // ice-cloud grey mass opacity [m2/kg]
    static constexpr double tau_cloud_cap = 0.4;    // per-layer cap on the cloud optical depth
    static constexpr double opac_cal      = 0.025;  // calibration on the whole gas+cloud opacity

    // Grey two-stream multi-layer solve. Fills radiation / epsilon / Q_rad.
    void run();

private:
    Planet& m;
};

template<class Planet>
void Radiation<Planet>::run(){
    const char* TAG = Planet::planet_tag();
    std::cout << std::endl << "      " << TAG << ": Radiation (grey CIA + CH4/NH3 + clouds)" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    const int im = m.im, jm = m.jm, km = m.km;
    const double t_ref = m.t_ref;
    const double g     = m.g;

    // The planet's own numbers. Copied into locals here rather than used as Planet:: expressions
    // below, so that the OpenMP region and the printf both see plain doubles.
    const double F_int       = Planet::rad_F_int();
    const double S_solar     = Planet::rad_S_solar();
    const double albedo_bond = Planet::rad_albedo_bond();
    const double x_H2        = Planet::rad_x_H2();
    const double x_He        = Planet::rad_x_He();

    // Composition-weighted collision-pair factor and the runtime CIA strength multiplier.
    const double comp = x_H2 * x_H2 + he_ratio * x_H2 * x_He;
    const double cia_mult  = ATPhys::env_double(TAG, "CIA_STRENGTH",     1.0);
    const double opac_mult = ATPhys::env_double(TAG, "OPACITY_STRENGTH", 1.0);
    const double cia_coeff = C_cia * cia_mult * comp / g;

    // Shortwave knobs. <TAG>_SOLAR=0 restores the pure-internal (F_int only) scaffold for A/B;
    // <TAG>_SOLAR_STRENGTH scales the absorbed flux; <TAG>_SW_TAU_PER_BAR moves the absorption
    // level (larger = absorbed higher up).
    const int    solar_on   = ATPhys::env_int   (TAG, "SOLAR",          1);
    const double solar_mult = ATPhys::env_double(TAG, "SOLAR_STRENGTH", 1.0);
    const double sw_tau_bar = ATPhys::env_double(TAG, "SW_TAU_PER_BAR", k_sw);
    const double pi_ = 3.14159265358979323846;

    // Column diagnostics: the mean outgoing longwave, and where tau = 1 measured down from the
    // top falls. That second number is the photosphere, and it is the calibration question this
    // scaffold exists to ask — see the banner.
    // tau1_t_sum carries the TEMPERATURE at that same level. OLR is sigma*T^4 evaluated where the
    // column becomes opaque, so the photosphere temperature is the one number that says whether an
    // over-emitting scheme is emitting wrongly or emitting a column that is too warm — and until
    // this was added, only the column MAXIMUM had ever been looked at, which answers neither.
    double olr_sum = 0.0, wsum = 0.0, tau1_p_sum = 0.0, tau1_t_sum = 0.0;
    long long n_tau1 = 0;

    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:olr_sum, wsum, tau1_p_sum, tau1_t_sum, n_tau1)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){

            const int i_base = std::max(m.surface_index(j, k), 0);  // deep boundary / first fluid cell
            const int i_top  = im - 1;                              // top boundary (to space)
            if(i_base >= i_top) continue;                           // degenerate column (all solid)

            // Thread-local column scratch.
            std::vector<double> B(im, 0.0);        // layer grey-body emission sigma*T^4 [W/m2]
            std::vector<double> eps(im, 0.0);      // layer emissivity
            std::vector<double> tau_l(im, 0.0);    // layer optical depth, kept for the tau=1 scan
            std::vector<double> Fd(im + 1, 0.0);   // downward flux at bottom-interface of layer i (index i); top = im
            std::vector<double> Fu(im + 1, 0.0);   // upward   flux at bottom-interface of layer i (index i); top = im

            // Layer emission and total (CIA + gas + cloud) emissivity.
            for(int i = i_base; i <= i_top; i++){
                const double T = m.t.x[i][j][k] * t_ref;
                B[i] = sigma * T * T * T * T;

                // Physical pressures [Pa] from the model's dimensionless p_stat (~bars).
                const double P_lo = m.p_stat.x[i][j][k] * P0_phys;                     // layer bottom (higher P)
                const double P_hi = (i < i_top) ? m.p_stat.x[i + 1][j][k] * P0_phys    // layer top (lower P)
                                                : 0.0;                                 // top layer: vacuum above
                double dP = P_lo - P_hi;                                               // layer pressure thickness
                if(dP < 0.0) dP = 0.0;
                const double P_mean = 0.5 * (P_lo + P_hi);
                const double dm = dP / g;                                              // hydrostatic mass path [kg/m2]

                // CIA layer optical depth  tau_cia = cia_coeff * (P/T) * dP   (tau ~ P^2/T).
                const double tau_cia = (T > 0.0) ? cia_coeff * (P_mean / T) * dP : 0.0;

                // The gas and cloud terms need q as a DIMENSIONLESS mass mixing ratio, but the
                // species fields are mass DENSITIES in kg/m3, so they are divided by the local
                // mixture density. rho_mix is read DIRECTLY, not through rho_at(): this is one of
                // the places that must always see the local field, never the constant fallback.
                const double rho_c = m.rho_mix.x[i][j][k];
                const double inv_rho = (rho_c > 0.0 && std::isfinite(rho_c)) ? 1.0 / rho_c : 0.0;
                const double q_ch4 = std::max(0.0, m.ch4.x[i][j][k]) * inv_rho;
                const double q_nh3 = std::max(0.0, m.nh3.x[i][j][k]) * inv_rho;
                const double tau_gas = opac_mult * opac_cal * (kappa_ch4 * q_ch4 + kappa_nh3 * q_nh3) * dm;

                // Cloud/ice continuum optical depth (all three condensing species), capped.
                const double q_liq = (std::max(0.0, m.h2o_cloud.x[i][j][k])
                                    + std::max(0.0, m.nh3_cloud.x[i][j][k])
                                    + std::max(0.0, m.ch4_cloud.x[i][j][k])) * inv_rho;
                const double q_ice = (std::max(0.0, m.h2o_ice.x[i][j][k])
                                    + std::max(0.0, m.nh3_ice.x[i][j][k])
                                    + std::max(0.0, m.ch4_ice.x[i][j][k])) * inv_rho;
                double tau_cloud = opac_mult * opac_cal * (kappa_cloud * q_liq + kappa_ice * q_ice) * dm;
                if(tau_cloud > tau_cloud_cap) tau_cloud = tau_cloud_cap;

                tau_l[i] = tau_cia + tau_gas + tau_cloud;
                eps[i]   = 1.0 - std::exp(-tau_l[i]);
            }

            // Where does tau = 1 from the top fall? A column that never reaches it is not counted,
            // which is why the printf reports the population as well as the mean.
            {
                double cum = 0.0;
                for(int i = i_top; i >= i_base; i--){
                    cum += tau_l[i];
                    if(cum >= 1.0){
                        tau1_p_sum += m.p_stat.x[i][j][k];
                        tau1_t_sum += m.t.x[i][j][k] * t_ref;
                        n_tau1++;
                        break;
                    }
                }
            }

            // --- Shortwave (solar) channel, independent of the thermal sweeps. ---
            // Absorbed insolation at the top of this column, then Beer-Lambert attenuation with an
            // optical depth measured DOWN FROM THE TOP so that Fsw at the top interface is exactly
            // F_sun_toa. Energy is then conserved by construction: everything absorbed in the
            // layers plus the residual reaching the deep boundary sums to F_sun_toa exactly.
            std::vector<double> Fsw(im + 2, 0.0);   // downward shortwave at bottom-interface of layer i
            double F_sun_toa = 0.0;
            if(solar_on != 0){
                // cos(latitude) = sin(colatitude); colatitude = j*pi/(jm-1), so j=0/jm-1 are the
                // poles (cos_lat = 0) and j=(jm-1)/2 the equator (cos_lat = 1).
                const double colat   = pi_ * (double)j / (double)(jm - 1);
                const double cos_lat = std::max(0.0, std::sin(colat));
                F_sun_toa = solar_mult * S_solar * (1.0 - albedo_bond) * cos_lat / pi_;

                // p_stat is ~bars; the interface above layer i_top is vacuum (p = 0), which is
                // also the zero point of the shortwave optical depth.
                for(int i = i_base; i <= i_top + 1; i++){
                    const double p_bar = (i <= i_top) ? m.p_stat.x[i][j][k] : 0.0;
                    Fsw[i] = F_sun_toa * std::exp(-sw_tau_bar * std::max(0.0, p_bar));
                }
            }

            // Downward sweep (top -> bottom). Fd[i] is the downward flux leaving the bottom of
            // layer i; Fd[i_top+1] is the incoming THERMAL flux at the top boundary, which is zero
            // — space is cold, and the solar term is the separate shortwave channel above.
            Fd[i_top + 1] = 0.0;
            for(int i = i_top; i >= i_base; i--)
                Fd[i] = Fd[i + 1] * (1.0 - eps[i]) + eps[i] * B[i];

            // Bottom boundary: net flux up through the deep interface must carry the intrinsic flux
            // PLUS whatever shortwave survived to depth (absorbed there and re-emitted in the
            // thermal channel), so the deep boundary emits Fd_bottom + F_int + Fsw_bottom.
            Fu[i_base] = Fd[i_base] + F_int + Fsw[i_base];

            // Upward sweep (bottom -> top). Fu[i+1] is the upward flux leaving the top of layer i
            // (= entering the bottom of layer i+1).
            for(int i = i_base; i <= i_top; i++)
                Fu[i + 1] = Fu[i] * (1.0 - eps[i]) + eps[i] * B[i];

            // Outgoing longwave at the top, area-weighted, for the balance check.
            {
                const double colat = pi_ * (double)j / (double)(jm - 1);
                const double wgt   = std::sin(colat);
                olr_sum += wgt * (Fu[i_top + 1] - Fd[i_top + 1]);
                wsum    += wgt;
            }

            // Interface net flux (up - down) and per-layer heating from its divergence. In
            // radiative equilibrium every interface net flux equals F_int and the divergence
            // (hence Q_rad) vanishes.
            for(int i = i_base; i <= i_top; i++){
                const double net_bot = Fu[i]     - Fd[i];       // net flux at bottom interface of layer i
                const double net_top = Fu[i + 1] - Fd[i + 1];   // net flux at top interface of layer i

                // Layer thickness in METRES. The models specify L_atm in km, so the raw layer
                // height difference would make Q_rad 1000x too large.
                double dz = (i < i_top) ? m.layer_thickness_m(i) : m.layer_thickness_m(i - 1);
                if(dz <= 0.0) dz = 1.0;   // guard degenerate layer heights

                // Shortwave heating: convergence of the downward solar beam in this layer. Fsw
                // decreases downward, so (Fsw[i+1] - Fsw[i]) >= 0 is absorbed energy and the
                // contribution to Q_rad is a HEATING term.
                const double sw_absorbed = Fsw[i + 1] - Fsw[i];            // [W/m2] absorbed in layer i
                const double Q_sw        = sw_absorbed / dz;               // [W/m3]

                m.epsilon.x[i][j][k]   = eps[i];
                m.radiation.x[i][j][k] = 0.5 * (net_bot + net_top);        // layer-centre net THERMAL flux [W/m2]
                m.Q_rad.x[i][j][k]     = (net_bot - net_top) / dz + Q_sw;  // total radiative heating [W/m3]
            }

            // Copy the deepest fluid value into any solid cells below so plots/BCs see a filled
            // column.
            for(int i = 0; i < i_base; i++){
                m.epsilon.x[i][j][k]   = m.epsilon.x[i_base][j][k];
                m.radiation.x[i][j][k] = m.radiation.x[i_base][j][k];
                m.Q_rad.x[i][j][k]     = 0.0;
            }
        }
    }

    const double olr  = (wsum > 0.0) ? olr_sum / wsum : 0.0;
    const double tau1 = (n_tau1 > 0) ? tau1_p_sum / double(n_tau1) : 0.0;
    const double sw_  = solar_mult * S_solar * (1.0 - albedo_bond) / 4.0;
    const double in_  = F_int + sw_;
    printf("      %s: radiation — mean OLR %.3f W/m2 against %.3f in (F_int %.2f + solar %.2f);"
           " thermal photosphere (tau=1) at %.4f bar in %lld of %d columns\n",
           TAG, olr, in_, F_int, sw_, tau1, n_tau1, jm * km);

    // The photosphere temperature, and the two things it is worth comparing against. B_tau1 is what
    // a blackbody AT that temperature emits: if it tracks the OLR, the scheme is faithfully
    // radiating the column it was handed and any excess is in the TEMPERATURE PROFILE, not here. If
    // it does not, the fault is in the two-stream sum. The two temperature gaps are DIFFERENT
    // questions and are easy to confuse: T_tau1 - T_eff(OLR) is the scheme's own excess in kelvin,
    // while T_tau1 - T_eff(in) is how far the COLUMN sits from the temperature the planet's energy
    // budget calls for. The second is only a defect on a planet actually near radiative
    // equilibrium; on one with a large internal flux the photosphere belongs above T_eff(in).
    const double tau1_T   = (n_tau1 > 0) ? tau1_t_sum / double(n_tau1) : 0.0;
    const double B_tau1   = sigma * tau1_T * tau1_T * tau1_T * tau1_T;
    const double T_eff_in = std::pow(in_ / sigma, 0.25);
    const double T_eff_ol = (olr > 0.0) ? std::pow(olr / sigma, 0.25) : 0.0;
    printf("      %s: radiation — photosphere T %.2f K emits %.3f W/m2 as a blackbody;"
           " T_eff(OLR) %.2f K vs T_eff(in) %.2f K\n",
           TAG, tau1_T, B_tau1, T_eff_ol, T_eff_in);

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Radiation\n", elapsed.count() * 1e-9);
    std::cout << "      " << TAG << ": Radiation ended" << std::endl;
}

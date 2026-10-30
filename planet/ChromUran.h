/*
 * ATURAN_CHROM=1 (2026-10-30, default 0) -- the CHROMOPHORE TRACER, ported from ATJUP's ChromJup.h (ATJUP_CHROM,
 * 2026-10-24) and ATSAT's ChromSat.h on the user's word ("can you put them also in the rest giants"), with the
 * ice giant's SOURCE: CH4 (the user's pick, 2026-10-30).
 *
 * On Jupiter and Saturn the colours come from NH3 broken up by UV high in the atmosphere. On Uranus NH3 is
 * frozen out far below the UV-lit top (NH3 vapour 0 above layer ~10 in the 1000-iteration runs of 10-30),
 * so a literal port would be empty at the top and unshielded at depth. What UV makes there instead is the
 * hydrocarbon haze from CH4 photolysis, and CH4 vapour reaches the lid. So the same equation, with CH4:
 *
 *     d chrom/dt |_chem = J * c_src * F_sun(theta) * A_uv  -  chrom / tau
 *
 *     c_src    ATURAN_CHROM_SRC=1 (default): CH4 vapour; 0: NH3 vapour + cloud + ice (the Jupiter recipe)
 *     F_sun    sin(colatitude): the diurnal-mean insolation shape, no season (1 at the equator).
 *              Uranus' 97.8 deg obliquity (its poles get more sun than its equator over a year) is NOT in it.
 *     A_uv     exp(-N_above / N_uv): the UV left after the source's VAPOUR column N_above [kg/m2] above the
 *              cell (from its upper face to the lid) has absorbed it -- the source shields itself.
 *              N_uv = ATURAN_CHROM_NUV [kg/m2], default 1.5e-5 with CH4: tau = 1 at 5.6e20 molecules/m2,
 *              the CH4 cross-section ~1.8e-17 cm2 at Lyman alpha; 3e-4 with NH3 (ATJUP's). No haze absorption.
 *     J        ATURAN_CHROM_J, effective production rate [1/s], default 1e-6 (ATJUP's)
 *     tau      ATURAN_CHROM_TAU [days of 86400 s], default 30 (ATJUP's)
 *
 * The real CH4 photolysis happens at microbar pressures, far above this model's lid; here the UV meets the
 * first CH4 below the lid, so with the default N_uv the tracer is made in a thin sheet at the top of the CH4
 * vapour and spreads only by transport and diffusion (no haze settling). Diffused with the source's Schmidt
 * number. A_uv is built once per iteration (chromUV(), before RungeKuttaUran) into chrom_uv.
 *
 * UNCALIBRATED ON PURPOSE, as on ATJUP: read chrom as a relative intensity. It takes nothing from CH4 and
 * feeds nothing back (no radiation, no density). Uranus writes no restart file, so there is no restart part.
 *
 * Off: chrom is never touched and stays 0 and no output carries it, so a default run is byte-identical.
 */

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace UranChrom {
    struct Par { double J, N_uv, tau_s; int src; };
    inline const Par& par(){
        static const Par p = [](){
            auto d = [](const char* n, double v){ const char* e = getenv(n); return e ? atof(e) : v; };
            const int src = (int)d("ATURAN_CHROM_SRC", 1.0);
            return Par{ d("ATURAN_CHROM_J", 1.0e-6), d("ATURAN_CHROM_NUV", src ? 1.5e-5 : 3.0e-4),
                        d("ATURAN_CHROM_TAU", 30.0) * 86400.0, src };
        }();
        return p;
    }
}

// chrom_uv = exp(-N_above / N_uv), N_above the source's vapour column from the upper face of the cell to the
// lid (the lid layer im-1 has no thickness, as in every column sum here).
inline void cUranusModel::chromUV(){
    if(!chrom_on()) return;
    const double n_uv = UranChrom::par().N_uv;
    const Array& gas = UranChrom::par().src ? ch4 : nh3;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            double n_above = 0.0;
            chrom_uv.x[im-1][j][k] = 1.0;
            for(int i = im - 2; i >= 0; i--){
                chrom_uv.x[i][j][k] = n_uv > 0.0 ? std::exp(-n_above / n_uv) : 1.0;
                n_above += std::max(0.0, gas.x[i][j][k]) * layer_thickness_m(i);
            }
        }
    }
}

inline double cUranusModel::chromSource(int i, int j, int k){
    const UranChrom::Par& p = UranChrom::par();
    const double c = p.src ? ch4.x[i][j][k] : nh3.x[i][j][k] + nh3_cloud.x[i][j][k] + nh3_ice.x[i][j][k];
    const double f_sun = std::max(0.0, std::sin(the.z[j]));
    const double rate = p.J * std::max(0.0, c) * f_sun * chrom_uv.x[i][j][k]
                      - (p.tau_s > 0.0 ? chrom.x[i][j][k] / p.tau_s : 0.0);     // kg/m3/s
    return rate * (L_atm * 1.0e3 / u_0);                                          // per nondimensional time
}

inline void cUranusModel::chromInfo(){
    if(!chrom_on()) return;
    const UranChrom::Par& p = UranChrom::par();
    printf("      ATURAN: ATURAN_CHROM on -- chromophore tracer: J %.3e 1/s, N_uv %.3e kg/m2 (%s vapour self-shielding), tau %.1f d, source %s\n",
           p.J, p.N_uv, p.src ? "CH4" : "NH3", p.tau_s / 86400.0, p.src ? "CH4 vapour" : "NH3 vapour + cloud + ice");
}

#include "cUranusModel.h"

using namespace std;

/*
 * Area-weighted horizontal mean, at each level, of the expression the buoyancy takes the anomaly
 * of. Mirrors ATSAT's computeBuoyancyRefLevel(), which mirrors ATJUP's.
 *
 * The mean has to be the mean OF the quantity whose anomaly is taken, or the anomaly no longer has
 * zero mean at that height — which is the whole point: with it subtracted, only horizontal density
 * contrasts drive vertical motion and hydrostatic balance is left to carry the mean. sin(theta) is
 * the spherical area weight; cells with no positive temperature are skipped, which also catches
 * NaN, and a non-finite contribution is dropped rather than poisoning the level.
 */
void cUranusModel::computeBuoyancyRefLevel(){
    if((int)buoy_ref_level.size() != im) buoy_ref_level.assign(im, 0.0);

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double sum = 0.0, wsum = 0.0;
        for(int j = 0; j < jm; j++){
            const double wgt = sin(the.z[j]);              // spherical area weight
            for(int k = 0; k < km; k++){
                if(!(t.x[i][j][k] > 0.0)) continue;        // also catches NaN
                const double b = g * p_stat.x[i][j][k]
                               / (r_mix * R_mix * t.x[i][j][k] * t_ref);
                if(!std::isfinite(b)) continue;
                sum  += wgt * b;
                wsum += wgt;
            }
        }
        buoy_ref_level[i] = (wsum > 0.0) ? sum / wsum : 0.0;
    }
}

/*
 * Hydrostatic pressure perturbation: p_hydro(r) = integral of the buoyancy ANOMALY from the base,
 * so d(p_hydro)/dr = buoyancy by construction. Ported from ATJUP's computeHydrostaticPressure().
 *
 * THE PROBLEM IT SOLVES, measured on this model rather than assumed. The buoyancy term in rhs_u is
 * ~1e5 too small because p_stat is in bar where the ideal-gas density needs pascals. Restoring that
 * factor was measured over 224 iterations and moved the answer 0.5 % — from OLR/in 30.093 to
 * 30.247 — because the pressure projection absorbs a larger radial body force and returns a
 * matching dpdr. The balance is enforced numerically either way, so nothing is left over to drive
 * vertical motion, and this model's radial velocity sits at ~1e-4 where ATSAT's is ~1e-1. With no
 * overturning, thermal diffusion flattens the column unopposed.
 *
 * Splitting it off analytically is what changes that: the radial force balances exactly and by
 * construction, and what enters the momentum equation is the HORIZONTAL gradient of p_hydro —
 * smaller by the 1/r that the horizontal derivative carries, and the part that physically drives a
 * circulation. p_dyn is left with the barotropic and non-hydrostatic remainder, which is the part
 * it can actually represent.
 *
 * The base is the DEEP boundary: the reference isobaric surface, where the gas is densest and
 * horizontal pressure contrasts are hardest to sustain, while the model top is an arbitrary cut
 * through a continuing atmosphere. ATURAN_HYDRO_REF=1 integrates downward from the top instead so
 * the choice can be measured; it is not the intended configuration.
 */
void cUranusModel::computeHydrostaticPressure(){
    static const bool ref_top = [](){
        const char* e = getenv("ATURAN_HYDRO_REF"); return e && atoi(e) != 0; }();
    // The same factor rhs_u puts on the buoyancy, read the same way, so the two cannot drift apart.
    static const double buoy_scale_local = [](){
        const char* e = getenv("ATURAN_BUOY_SCALE"); return e ? atof(e) : 1.0; }();

    // L_atm is in KILOMETRES (see the config: "extension of the troposhere in km"), and the length
    // scale this nondimensionalisation needs is in metres. ATJUP writes the factor
    // 1.0e5 * (L_atm * 1.0e3) / (u_0 * u_0); the port that brought the split here dropped the
    // 1.0e3, so p_hydro came out 1000x too small and the split it feeds did essentially nothing.
    // ATURAN_HYDRO_ND_KM=1 restores the dropped-factor version so the two can be attributed apart.
    static const bool nd_km = [](){
        const char* e = getenv("ATURAN_HYDRO_ND_KM"); return e && atoi(e) != 0; }();
    const double nd = 1.0e5 * (L_atm * (nd_km ? 1.0 : 1.0e3)) / (u_0 * u_0);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            auto buoy = [&](int i)->double{
                const double T = t.x[i][j][k];
                if(!(T > 0.0)) return 0.0;
                const double b = g * p_stat.x[i][j][k] / (r_mix * R_mix * T * t_ref)
                               - buoy_ref_level[i];
                const double f = -nd * buoy_scale_local * buoyancy * b;
                return std::isfinite(f) ? f : 0.0;
            };
            if(!ref_top){
                p_hydro.x[0][j][k] = 0.0;
                for(int i = 1; i < im; i++)
                    p_hydro.x[i][j][k] = p_hydro.x[i-1][j][k]
                                       + 0.5 * (buoy(i-1) + buoy(i)) * dr;
            } else {
                p_hydro.x[im-1][j][k] = 0.0;
                for(int i = im-2; i >= 0; i--)
                    p_hydro.x[i][j][k] = p_hydro.x[i+1][j][k]
                                       - 0.5 * (buoy(i+1) + buoy(i)) * dr;
            }
        }
    }

}

void cUranusModel::RungeKuttaUran(){
    cout << endl << "      ATURAN: RungeKuttaUran" << endl;


    auto begin = std::chrono::high_resolution_clock::now();

    // The buoyancy base state and the hydrostatic pressure built from it, refreshed once per RK4
    // step before any stage reads them. Both are inert unless ATURAN_HYDRO_SPLIT is set.
    computeBuoyancyRefLevel();
    computeHydrostaticPressure();

    // Precompute sin/cos tables — depend only on j.
    // sinthe is clamped to a minimum to prevent 1/sin²θ blow-up near the poles.
    const double sinthe_floor = sinthe_min();   // the model's polar metric floor, cUranusModel.h
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = std::max(sinthe_floor, std::abs(sin(the.z[j])));
        costhe_tbl[j] = cos(the.z[j]);
    }

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr   * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    // Turbulence clamps, used only when the closure is on. k* has a physical ceiling; dis* has a
    // FLOOR because it appears in denominators throughout the closure (nue = k/dis among them) and
    // a zero there is an infinity one step later.
    const bool turb_on_rk = turb_active;
    const double tke_max_nd = 1000.0 / (u_0 * u_0);   // 1000 m2/s2
    constexpr double dis_min_nd = 1.0e-10;            // matches the closure's dis_min

    // Physical bounds, divided by t_ref here rather than written as nondimensional literals.
    // See cUranusModel.h for why, and for what the model actually reaches.
    const double t_min = t_min_K() / t_ref;
    const double t_max = t_max_K() / t_ref;

    /*
     * FOUR BARRIER-SEPARATED STAGES, which is the whole point of this file.
     *
     * What was here ran all four RK4 stages inside ONE parallel loop over cells, holding k1..k4 as
     * per-thread scalars and overwriting t,u,v,w at the cell as it went — while RHSUran
     * differentiates those same live fields at i+-1, j+-1, k+-1, cells other threads are writing
     * at the same moment. Two runs of the same binary at 24 threads differed in 4 of 7 output
     * files; that race is why.
     *
     * Each stage is now two passes with a barrier between them, and the running sum
     * k1 + 2k2 + 2k3 + k4 lives in the acc_* arrays rather than in registers:
     *
     *   pass A   evaluate RHSUran over every cell, reading a state no one is writing
     *   pass B   fold rhs_* into acc_*, then form the NEXT stage's input from the n-copies
     *
     * c_in is the offset of the next stage's input from y_n, wgt this stage's weight in the sum.
     * Stage 3 forms no input, which is why the write is guarded. This is ATSAT's 71082e7 applied
     * to Uranus, and the acceptance test is its: reproducibility at any thread count.
     */
    for(int stage = 0; stage < 4; stage++){

        const double c_in = (stage == 0 || stage == 1) ? 0.5 * dt : (stage == 2 ? dt : 0.0);
        const double wgt  = (stage == 0 || stage == 3) ? 1.0 : 2.0;

        // Tells RHSUran which stage it is inside, so item 1's budget instrument can record the
        // decomposition at stage 0 only — the evaluation at y_n. Read-only in the loop below.
        tbud_stage = stage;

        // ---- pass A: evaluate the right-hand sides ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 1; j < jm-1; j++){
            CellGeometry geo;
            // metricRadius(): identity unless ATURAN_METRIC_RADIUS is set, so this is inert by
            // default. See the note at metricRadius() in the model header — this model has no
            // metric radius, i.e. rad.z runs 1..2 and the metric puts the planet's surface
            // L_atm from its centre instead of R, making every HORIZONTAL derivative
            // R/L_atm times too large. ATJUP measured the same defect on itself as 499x and
            // fixed it by shifting rad.z at initialisation; ATSAT routes it through this
            // accessor, as here.
            geo.rm      = metricRadius(rad.z[i]);
            geo.rm2     = geo.rm * geo.rm;
            geo.sinthe  = sinthe_tbl[j];
            geo.sinthe2 = geo.sinthe * geo.sinthe;
            geo.costhe  = costhe_tbl[j];
            geo.cotanthe             = geo.costhe / geo.sinthe;
            geo.inv_rm               = 1.0 / geo.rm;
            geo.inv_rm2              = 1.0 / geo.rm2;
            geo.inv_rmsinthe         = 1.0 / (geo.rm * geo.sinthe);
            geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
            geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
            geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
            geo.inv_2dr   = inv_2dr;
            geo.inv_2dthe = inv_2dthe;
            geo.inv_2dphi = inv_2dphi;
            geo.inv_dr2   = inv_dr2;
            geo.inv_dthe2 = inv_dthe2;
            geo.inv_dphi2 = inv_dphi2;

                for(int k = 1; k < km-1; k++)
                    cUranusModel::RHSUran(i, j, k, geo);
            }
        }

        // ---- pass B: fold into the accumulator, then form the next stage's input ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 1; j < jm-1; j++){
                for(int k = 1; k < km-1; k++){
                    acc_t.x[i][j][k] = (stage == 0) ? wgt * rhs_t.x[i][j][k]
                                       : acc_t.x[i][j][k] + wgt * rhs_t.x[i][j][k];
                    acc_u.x[i][j][k] = (stage == 0) ? wgt * rhs_u.x[i][j][k]
                                       : acc_u.x[i][j][k] + wgt * rhs_u.x[i][j][k];
                    acc_v.x[i][j][k] = (stage == 0) ? wgt * rhs_v.x[i][j][k]
                                       : acc_v.x[i][j][k] + wgt * rhs_v.x[i][j][k];
                    acc_w.x[i][j][k] = (stage == 0) ? wgt * rhs_w.x[i][j][k]
                                       : acc_w.x[i][j][k] + wgt * rhs_w.x[i][j][k];
                    acc_ch4.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4.x[i][j][k]
                                       : acc_ch4.x[i][j][k] + wgt * rhs_ch4.x[i][j][k];
                    acc_ch4_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_cloud.x[i][j][k]
                                       : acc_ch4_cloud.x[i][j][k] + wgt * rhs_ch4_cloud.x[i][j][k];
                    acc_ch4_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_ice.x[i][j][k]
                                       : acc_ch4_ice.x[i][j][k] + wgt * rhs_ch4_ice.x[i][j][k];
                    acc_h2o.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o.x[i][j][k]
                                       : acc_h2o.x[i][j][k] + wgt * rhs_h2o.x[i][j][k];
                    acc_h2o_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_cloud.x[i][j][k]
                                       : acc_h2o_cloud.x[i][j][k] + wgt * rhs_h2o_cloud.x[i][j][k];
                    acc_h2o_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_ice.x[i][j][k]
                                       : acc_h2o_ice.x[i][j][k] + wgt * rhs_h2o_ice.x[i][j][k];
                    acc_h2s.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s.x[i][j][k]
                                       : acc_h2s.x[i][j][k] + wgt * rhs_h2s.x[i][j][k];
                    acc_h2s_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s_cloud.x[i][j][k]
                                       : acc_h2s_cloud.x[i][j][k] + wgt * rhs_h2s_cloud.x[i][j][k];
                    acc_h2s_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s_ice.x[i][j][k]
                                       : acc_h2s_ice.x[i][j][k] + wgt * rhs_h2s_ice.x[i][j][k];
                    acc_nh3.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3.x[i][j][k]
                                       : acc_nh3.x[i][j][k] + wgt * rhs_nh3.x[i][j][k];
                    acc_nh3_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_cloud.x[i][j][k]
                                       : acc_nh3_cloud.x[i][j][k] + wgt * rhs_nh3_cloud.x[i][j][k];
                    acc_nh3_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_ice.x[i][j][k]
                                       : acc_nh3_ice.x[i][j][k] + wgt * rhs_nh3_ice.x[i][j][k];
                    acc_nh4sh.x[i][j][k] = (stage == 0) ? wgt * rhs_nh4sh.x[i][j][k]
                                       : acc_nh4sh.x[i][j][k] + wgt * rhs_nh4sh.x[i][j][k];
                    if(turb_on_rk){
                        acc_tke.x[i][j][k] = (stage == 0) ? wgt * rhs_tke.x[i][j][k]
                                           : acc_tke.x[i][j][k] + wgt * rhs_tke.x[i][j][k];
                        acc_dis.x[i][j][k] = (stage == 0) ? wgt * rhs_dis.x[i][j][k]
                                           : acc_dis.x[i][j][k] + wgt * rhs_dis.x[i][j][k];
                    }

                    if(stage < 3){
                        t.x[i][j][k] = tn.x[i][j][k] + c_in * rhs_t.x[i][j][k];
                        u.x[i][j][k] = un.x[i][j][k] + c_in * rhs_u.x[i][j][k];
                        v.x[i][j][k] = vn.x[i][j][k] + c_in * rhs_v.x[i][j][k];
                        w.x[i][j][k] = wn.x[i][j][k] + c_in * rhs_w.x[i][j][k];
                        ch4.x[i][j][k] = ch4n.x[i][j][k] + c_in * rhs_ch4.x[i][j][k];
                        ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + c_in * rhs_ch4_cloud.x[i][j][k];
                        ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + c_in * rhs_ch4_ice.x[i][j][k];
                        h2o.x[i][j][k] = h2on.x[i][j][k] + c_in * rhs_h2o.x[i][j][k];
                        h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + c_in * rhs_h2o_cloud.x[i][j][k];
                        h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + c_in * rhs_h2o_ice.x[i][j][k];
                        h2s.x[i][j][k] = h2sn.x[i][j][k] + c_in * rhs_h2s.x[i][j][k];
                        h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + c_in * rhs_h2s_cloud.x[i][j][k];
                        h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + c_in * rhs_h2s_ice.x[i][j][k];
                        nh3.x[i][j][k] = nh3n.x[i][j][k] + c_in * rhs_nh3.x[i][j][k];
                        nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + c_in * rhs_nh3_cloud.x[i][j][k];
                        nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + c_in * rhs_nh3_ice.x[i][j][k];
                        nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + c_in * rhs_nh4sh.x[i][j][k];
                        if(turb_on_rk){
                            tke.x[i][j][k] = AtomUtils::clamp(tken.x[i][j][k]
                                + c_in * rhs_tke.x[i][j][k], 0.0, tke_max_nd);
                            dis.x[i][j][k] = std::max(disn.x[i][j][k]
                                + c_in * rhs_dis.x[i][j][k], dis_min_nd);
                        }
                    }
                }
            }
        }
    }

    // ---- the step itself: y_{n+1} = y_n + dt/6 * (k1 + 2k2 + 2k3 + k4) ----
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            for(int k = 1; k < km-1; k++){
                t.x[i][j][k] = AtomUtils::clamp(
                    tn.x[i][j][k] + dt * acc_t.x[i][j][k] / 6.0, t_min, t_max);
                u.x[i][j][k] = un.x[i][j][k] + dt * acc_u.x[i][j][k] / 6.0;
                v.x[i][j][k] = vn.x[i][j][k] + dt * acc_v.x[i][j][k] / 6.0;
                w.x[i][j][k] = wn.x[i][j][k] + dt * acc_w.x[i][j][k] / 6.0;
                ch4.x[i][j][k] = ch4n.x[i][j][k] + dt * acc_ch4.x[i][j][k] / 6.0;
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + dt * acc_ch4_cloud.x[i][j][k] / 6.0;
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + dt * acc_ch4_ice.x[i][j][k] / 6.0;
                h2o.x[i][j][k] = h2on.x[i][j][k] + dt * acc_h2o.x[i][j][k] / 6.0;
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + dt * acc_h2o_cloud.x[i][j][k] / 6.0;
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + dt * acc_h2o_ice.x[i][j][k] / 6.0;
                h2s.x[i][j][k] = h2sn.x[i][j][k] + dt * acc_h2s.x[i][j][k] / 6.0;
                h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + dt * acc_h2s_cloud.x[i][j][k] / 6.0;
                h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + dt * acc_h2s_ice.x[i][j][k] / 6.0;
                nh3.x[i][j][k] = nh3n.x[i][j][k] + dt * acc_nh3.x[i][j][k] / 6.0;
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + dt * acc_nh3_cloud.x[i][j][k] / 6.0;
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + dt * acc_nh3_ice.x[i][j][k] / 6.0;
                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + dt * acc_nh4sh.x[i][j][k] / 6.0;
                if(turb_on_rk){
                    tke.x[i][j][k] = AtomUtils::clamp(std::max(tken.x[i][j][k]
                        + dt * acc_tke.x[i][j][k] / 6.0, 0.0), 0.0, tke_max_nd);
                    dis.x[i][j][k] = std::max(disn.x[i][j][k]
                        + dt * acc_dis.x[i][j][k] / 6.0, dis_min_nd);
                }
            }
        }
    }
}

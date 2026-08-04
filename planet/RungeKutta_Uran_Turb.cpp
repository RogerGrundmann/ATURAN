#include "cUranusModel.h"

using namespace std;

void cUranusModel::RungeKuttaUran(){
    cout << endl << "      ATURAN: RungeKuttaUran" << endl;


    auto begin = std::chrono::high_resolution_clock::now();

    // Precompute sin/cos tables — depend only on j.
    // sinthe is clamped to a minimum to prevent 1/sin²θ blow-up near the poles.
    constexpr double sinthe_min = 0.4;
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = std::max(sinthe_min, std::abs(sin(the.z[j])));
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

        // ---- pass A: evaluate the right-hand sides ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 1; j < jm-1; j++){
            CellGeometry geo;
            geo.rm      = rad.z[i];
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

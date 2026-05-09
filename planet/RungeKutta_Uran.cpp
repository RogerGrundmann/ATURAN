#include "cUranusModel.h"

using namespace std;

void cUranusModel::RungeKuttaUran(){
    cout << endl << "      ATURAN: RungeKuttaUran" << endl;
    double 
        kt1, ku1, kv1, kw1, kc1, kcloud1, kice1, 
        kh2s1, kh2s_cloud1, kh2s_ice1, 
        knh31, knh3_cloud1, knh3_ice1, 
        kch41, kch4_cloud1, kch4_ice1, 
//        knh4sh1, knh4sh_cloud1,
        knh4sh1,
        
        kt2, ku2, kv2, kw2, kc2, kcloud2, kice2, 
        kh2s2, kh2s_cloud2, kh2s_ice2, 
        knh32, knh3_cloud2, knh3_ice2, 
        kch42, kch4_cloud2, kch4_ice2, 
//        knh4sh2, knh4sh_cloud2,
        knh4sh2,

        kt3, ku3, kv3, kw3, kc3, kcloud3, kice3, 
        kh2s3, kh2s_cloud3, kh2s_ice3, 
        knh33, knh3_cloud3, knh3_ice3, 
        kch43, kch4_cloud3, kch4_ice3, 
//        knh4sh3, knh4sh_cloud3,
        knh4sh3,

        kt4, ku4, kv4, kw4, kc4, kcloud4, kice4, 
        kh2s4, kh2s_cloud4, kh2s_ice4,
        knh34, knh3_cloud4, knh3_ice4,
        kch44, kch4_cloud4, kch4_ice4,
//        knh4sh4, knh4sh_cloud4,
        knh4sh4;


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

    // Physical temperature bounds in normalised units (T / t_ref).
    // t_min prevents negative absolute temperatures that blow up the buoyancy
    // denominator; t_max caps runaway heating to ~2× the surface temperature.
    const double t_min = 0.1;                        // 7.6 K physical
    const double t_max = 10.0;                       // 760 K physical

    #pragma omp parallel for private(kt1, ku1, kv1, kw1, kc1, kcloud1, kice1, kch41, kch4_cloud1, kch4_ice1, kh2s1, kh2s_cloud1, kh2s_ice1, knh31, knh3_cloud1, knh3_ice1, knh4sh1, kt2, ku2, kv2, kw2, kc2, kcloud2, kice2, kch42, kch4_cloud2, kch4_ice2, kh2s2, kh2s_cloud2, kh2s_ice2, knh32, knh3_cloud2, knh3_ice2, knh4sh2, kt3, ku3, kv3, kw3, kc3, kcloud3, kice3, kch43, kch4_cloud3, kch4_ice3, kh2s3, kh2s_cloud3, kh2s_ice3, knh33, knh3_cloud3, knh3_ice3, knh4sh3, kt4, ku4, kv4, kw4, kc4, kcloud4, kice4, kch44, kch4_cloud4, kch4_ice4, kh2s4, kh2s_cloud4, kh2s_ice4, knh34, knh3_cloud4, knh3_ice4, knh4sh4)

    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){

            // Build geometry struct once per (i,j) — reused for all k
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

            for(int k = 1; k < km-1; k++){

                cUranusModel::RHSUran(i, j, k, geo);

                kt1 = rhs_t.x[i][j][k];

                ku1 = rhs_u.x[i][j][k];
                kv1 = rhs_v.x[i][j][k];
                kw1 = rhs_w.x[i][j][k];

                kc1 = rhs_h2o.x[i][j][k];
                kcloud1 = rhs_h2o_cloud.x[i][j][k];
                kice1 = rhs_h2o_ice.x[i][j][k];

                kh2s1 = rhs_h2s.x[i][j][k];
                kh2s_cloud1 = rhs_h2s_cloud.x[i][j][k];
                kh2s_ice1 = rhs_h2s_ice.x[i][j][k];

                knh31 = rhs_nh3.x[i][j][k];
                knh3_cloud1 = rhs_nh3_cloud.x[i][j][k];
                knh3_ice1 = rhs_nh3_ice.x[i][j][k];

                kch41 = rhs_ch4.x[i][j][k];
                kch4_cloud1 = rhs_ch4_cloud.x[i][j][k];
                kch4_ice1 = rhs_ch4_ice.x[i][j][k];

                knh4sh1 = rhs_nh4sh.x[i][j][k];
//                knh4sh_cloud1 = rhs_nh4sh_cloud.x[i][j][k];

                t.x[i][j][k] = AtomUtils::clamp(tn.x[i][j][k] + kt1 * 0.5 * dt, t_min, t_max);

                u.x[i][j][k] = un.x[i][j][k] + ku1 * 0.5 * dt;
                v.x[i][j][k] = vn.x[i][j][k] + kv1 * 0.5 * dt;
                w.x[i][j][k] = wn.x[i][j][k] + kw1 * 0.5 * dt;

                h2o.x[i][j][k] = h2on.x[i][j][k] + kc1 * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + kcloud1 * 0.5 * dt;
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + kice1 * 0.5 * dt;

                h2s.x[i][j][k] = h2sn.x[i][j][k] + kh2s1 * 0.5 * dt;
                h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + kcloud1 * 0.5 * dt;
                h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + kice1 * 0.5 * dt;

                nh3.x[i][j][k] = nh3n.x[i][j][k] + knh31 * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + knh3_cloud1 * 0.5 * dt;
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + knh3_ice1 * 0.5 * dt;

                ch4.x[i][j][k] = ch4n.x[i][j][k] + kch41 * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + kch4_cloud1 * 0.5 * dt;
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + kch4_ice1 * 0.5 * dt;

                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + knh4sh1 * 0.5 * dt;
//                nh4sh_cloud.x[i][j][k] = nh4sh_cloudn.x[i][j][k] + knh4sh_cloud1 * 0.5 * dt;

                cUranusModel::RHSUran(i, j, k, geo);

                kt2 = rhs_t.x[i][j][k];

                ku2 = rhs_u.x[i][j][k];
                kv2 = rhs_v.x[i][j][k];
                kw2 = rhs_w.x[i][j][k];

                kc2 = rhs_h2o.x[i][j][k];
                kcloud2 = rhs_h2o_cloud.x[i][j][k];
                kice2 = rhs_h2o_ice.x[i][j][k];

                kh2s2 = rhs_h2s.x[i][j][k];
                kh2s_cloud2 = rhs_h2s_cloud.x[i][j][k];
                kh2s_ice2 = rhs_h2s_ice.x[i][j][k];

                knh32 = rhs_nh3.x[i][j][k];
                knh3_cloud2 = rhs_nh3_cloud.x[i][j][k];
                knh3_ice2 = rhs_nh3_ice.x[i][j][k];

                kch42 = rhs_ch4.x[i][j][k];
                kch4_cloud2 = rhs_ch4_cloud.x[i][j][k];
                kch4_ice2 = rhs_ch4_ice.x[i][j][k];

                knh4sh2 = rhs_nh4sh.x[i][j][k];
//                knh4sh_cloud2 = rhs_nh4sh_cloud.x[i][j][k];

                t.x[i][j][k] = AtomUtils::clamp(tn.x[i][j][k] + kt2 * 0.5 * dt, t_min, t_max);

                u.x[i][j][k] = un.x[i][j][k] + ku2 * 0.5 * dt;
                v.x[i][j][k] = vn.x[i][j][k] + kv2 * 0.5 * dt;
                w.x[i][j][k] = wn.x[i][j][k] + kw2 * 0.5 * dt;

                h2o.x[i][j][k] = h2on.x[i][j][k] + kc2 * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + kcloud2 * 0.5 * dt;
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + kice2 * 0.5 * dt;

                h2s.x[i][j][k] = h2sn.x[i][j][k] + kh2s2 * 0.5 * dt;
                h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + kcloud2 * 0.5 * dt;
                h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + kice2 * 0.5 * dt;

                nh3.x[i][j][k] = nh3n.x[i][j][k] + knh32 * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + knh3_cloud2 * 0.5 * dt;
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + knh3_ice2 * 0.5 * dt;

                ch4.x[i][j][k] = ch4n.x[i][j][k] + kch42 * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + kch4_cloud2 * 0.5 * dt;
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + kch4_ice2 * 0.5 * dt;

                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + knh4sh2 * 0.5 * dt;

                cUranusModel::RHSUran(i, j, k, geo);

                kt3 = rhs_t.x[i][j][k];

                ku3 = rhs_u.x[i][j][k];
                kv3 = rhs_v.x[i][j][k];
                kw3 = rhs_w.x[i][j][k];

                kc3 = rhs_h2o.x[i][j][k];
                kcloud3 = rhs_h2o_cloud.x[i][j][k];
                kice3 = rhs_h2o_ice.x[i][j][k];

                kh2s3 = rhs_h2s.x[i][j][k];
                kh2s_cloud3 = rhs_h2s_cloud.x[i][j][k];
                kh2s_ice3 = rhs_h2s_ice.x[i][j][k];

                knh33 = rhs_nh3.x[i][j][k];
                knh3_cloud3 = rhs_nh3_cloud.x[i][j][k];
                knh3_ice3 = rhs_nh3_ice.x[i][j][k];

                kch43 = rhs_ch4.x[i][j][k];
                kch4_cloud3 = rhs_ch4_cloud.x[i][j][k];
                kch4_ice3 = rhs_ch4_ice.x[i][j][k];

                knh4sh3 = rhs_nh4sh.x[i][j][k];

                t.x[i][j][k] = AtomUtils::clamp(tn.x[i][j][k] + kt3 * dt, t_min, t_max);

                u.x[i][j][k] = un.x[i][j][k] + ku3 * dt;
                v.x[i][j][k] = vn.x[i][j][k] + kv3 * dt;
                w.x[i][j][k] = wn.x[i][j][k] + kw3 * dt;

                h2o.x[i][j][k] = h2on.x[i][j][k] + kc3 * dt;
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + kcloud3 * dt;
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + kice3 * dt;

                h2s.x[i][j][k] = h2sn.x[i][j][k] + kh2s3 * dt;
                h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] + kcloud3 * dt;
                h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] + kice3 * dt;

                nh3.x[i][j][k] = nh3n.x[i][j][k] + knh33 * dt;
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + knh3_cloud3 * dt;
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + knh3_ice3 * dt;

                ch4.x[i][j][k] = ch4n.x[i][j][k] + kch43 * dt;
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + kch4_cloud3 * dt;
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + kch4_ice3 * dt;

                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + knh4sh3 * dt;

                cUranusModel::RHSUran(i, j, k, geo);

                kt4 = rhs_t.x[i][j][k];

                ku4 = rhs_u.x[i][j][k];
                kv4 = rhs_v.x[i][j][k];
                kw4 = rhs_w.x[i][j][k];


                kc4 = rhs_h2o.x[i][j][k];
                kcloud4 = rhs_h2o_cloud.x[i][j][k];
                kice4 = rhs_h2o_ice.x[i][j][k];

                kh2s4 = rhs_h2s.x[i][j][k];
                kh2s_cloud4 = rhs_h2s_cloud.x[i][j][k];
                kh2s_ice4 = rhs_h2s_ice.x[i][j][k];

                knh34 = rhs_nh3.x[i][j][k];
                knh3_cloud4 = rhs_nh3_cloud.x[i][j][k];
                knh3_ice4 = rhs_nh3_ice.x[i][j][k];

                kch44 = rhs_ch4.x[i][j][k];
                kch4_cloud4 = rhs_ch4_cloud.x[i][j][k];
                kch4_ice4 = rhs_ch4_ice.x[i][j][k];

                knh4sh4 = rhs_nh4sh.x[i][j][k];

                t.x[i][j][k] = AtomUtils::clamp(tn.x[i][j][k] + dt * (kt1 + 2.0 * kt2
                    + 2.0 * kt3 + kt4) / 6.0, t_min, t_max);
                u.x[i][j][k] = un.x[i][j][k] + dt * (ku1 + 2.0 * ku2 
                    + 2.0 * ku3 + ku4)/6.0;
                v.x[i][j][k] = vn.x[i][j][k] + dt * (kv1 + 2.0 * kv2 
                    + 2.0 * kv3 + kv4)/6.0;
                w.x[i][j][k] = wn.x[i][j][k] + dt * (kw1 + 2.0 * kw2 
                    + 2.0 * kw3 + kw4)/6.0;

                h2o.x[i][j][k] = h2on.x[i][j][k] 
                    + dt * (kc1 + 2.0 * kc2 
                    + 2.0 * kc3 + kc4)/6.0;
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] 
                    + dt * (kcloud1 + 2.0 * kcloud2 
                    + 2.0 * kcloud3 + kcloud4)/6.0;
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] 
                    + dt * (kice1 + 2.0 * kice2 
                    + 2.0 * kice3 + kice4)/6.0;

                h2s.x[i][j][k] = h2sn.x[i][j][k] 
                    + dt * (kh2s1 + 2.0 * kh2s2 
                    + 2.0 * kh2s3 + kh2s4)/6.0;

                h2s_cloud.x[i][j][k] = h2s_cloudn.x[i][j][k] 
                    + dt * (kh2s_cloud1 + 2.0 * kh2s_cloud2 
                    + 2.0 * kh2s_cloud3 + kh2s_cloud4)/6.0;
                h2s_ice.x[i][j][k] = h2s_icen.x[i][j][k] 
                    + dt * (kh2s_ice1 + 2.0 * kh2s_ice2 
                    + 2.0 * kh2s_ice3 + kh2s_ice4)/6.0;

                nh3.x[i][j][k] = nh3n.x[i][j][k] 
                    + dt * (knh31 + 2.0 * knh32 
                    + 2.0 * knh33 + knh34)/6.0;
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] 
                    + dt * (knh3_cloud1 + 2.0 * knh3_cloud2 
                    + 2.0 * knh3_cloud3 + knh3_cloud4)/6.0;
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] 
                    + dt * (knh3_ice1 + 2.0 * knh3_ice2 
                    + 2.0 * knh3_ice3 + knh3_ice4)/6.0;

                ch4.x[i][j][k] = ch4n.x[i][j][k] 
                    + dt * (kch41 + 2.0 * kch42 
                    + 2.0 * kch43 + kch44)/6.0;
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] 
                    + dt * (kch4_cloud1 + 2.0 * kch4_cloud2 
                    + 2.0 * kch4_cloud3 + kch4_cloud4)/6.0;
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] 
                    + dt * (kch4_ice1 + 2.0 * kch4_ice2 
                    + 2.0 * kch4_ice3 + kch4_ice4)/6.0;

                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] 
                    + dt * (knh4sh1 + 2.0 * knh4sh2 
                    + 2.0 * knh4sh3 + knh4sh4)/6.0;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: RungeKuttaUran ended" << endl;
    return;
}
/*
*
*/

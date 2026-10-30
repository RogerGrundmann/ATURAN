/*
 * 2026-10-30, row ICE-FROZEN-FLOW (user: "test a small longitude-varying temperature seed"). Both default off.
 *
 * ATURAN_T_SEED=<A> [K]: once, before the time loop, T += A cos(m phi) sin(colat) in every interior layer
 *   (not the floor i=0, not the lid im-1); m = ATURAN_T_SEED_M (default 4). Deterministic, so two runs agree.
 *   Breaks the zonal symmetry of the start; whether the flow then grows or damps the seed is the question.
 * ATURAN_ZDEV=<n>: every n iterations, print-only, the largest departure from the zonal mean of T [K], u, v, w
 *   [m/s] over the whole shell, with where it sits (layer i, latitude). Zonal mean over the km-1 distinct
 *   meridians.
 */
#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>

inline void cUranusModel::tSeed(){
    static const double A = [](){ const char* e = getenv("ATURAN_T_SEED"); return e ? atof(e) : 0.0; }();
    static const int m = [](){ const char* e = getenv("ATURAN_T_SEED_M"); return e ? atoi(e) : 4; }();
    if(A == 0.0) return;
    for(int i = 1; i < im - 1; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                t.x[i][j][k] += A / t_ref * std::cos(m * phi.z[k]) * std::sin(the.z[j]);
    printf("      ATURAN: T_SEED %.3f K, m = %d, cos(m phi) sin(colat) in layers 1..%d\n", A, m, im - 2);
}

inline void cUranusModel::zdevDiag(){
    static const int every = [](){ const char* e = getenv("ATURAN_ZDEV"); return e ? atoi(e) : 0; }();
    if(every <= 0 || iter_n % every != 0) return;
    struct F { const char* name; Array* a; double s; } f[] = { { "T", &t, t_ref }, { "u", &u, u_0 }, { "v", &v, u_0 }, { "w", &w, u_0 } };
    printf("      ATURAN: ZDEV it %d max |X - zonal mean|:", iter_n);
    for(const F& q : f){
        double best = 0.0; int bi = 0, bj = 0;
        for(int i = 0; i < im; i++) for(int j = 0; j < jm; j++){
            double s = 0.0;
            for(int k = 0; k < km - 1; k++) s += q.a->x[i][j][k];
            const double mean = s / (km - 1);
            for(int k = 0; k < km - 1; k++){
                const double d = std::fabs(q.a->x[i][j][k] - mean) * q.s;
                if(d > best){ best = d; bi = i; bj = j; }
            }
        }
        printf("  %s %.4e (i %d, lat %d)", q.name, best, bi, 90 - bj);
    }
    printf("\n");
}

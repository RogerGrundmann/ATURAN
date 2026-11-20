/*
 * SHARED HEADER (planet/SHARED.md5): byte-identical in ATJUP, ATSAT, ATNEPT and ATURAN.
 *
 * <TAG>_PROC_DIAG=<n> -- PROCESS fields, OUTPUT-ONLY, default 0 = off. Written for ATSAT on 2026-11-20
 * (ATSAT c5a6f09, P2 of the visualisation proposal; user: "start with P1, then P2"), shared the same day
 * (user: "port ProcSat to the other giants").
 *
 * WHY. Every VTK file holds the STATE (vapour, cloud, ice, NH4SH, w, ...). Where NH4SH forms or decomposes,
 * where vapour condenses or evaporates, where precipitation takes condensate out and where the moist convective
 * adjustment fires could only be guessed from two snapshots. This instrument charges every change of t and
 * of the species to the stage of the time loop that made it -- as the TB / CWB diagnostics do for level means
 * and columns -- but keeps the field, zonally averaged:
 *
 *   stage          fields (zonal mean, time mean over the n iterations before the write)
 *   SatAdj         dT_satadj; cond_<g> (> 0) and evap_<g> (< 0) = cell-wise sign split of d(cloud + ice)
 *   NH4SH equil.   dT_chem; nh4sh_form (> 0), nh4sh_decomp (< 0)   (Chemistry*::EquilibriumNH4SH, ATJUP / ATSAT)
 *   Precipitation  the precipitation does not write the fields: it leaves tendencies that RHS applies over the
 *                  next iterations (S_precip_*, nondimensional = kg/m3/s x L_atm[m] / u_0; Q_precip W/m3), all
 *                  times <TAG>_PRECIP_COUPLING. Sampled after each call and averaged over the calls:
 *                  conv_<g> = (S_c + S_i) u_0 / L (< 0: condensate taken out); revap_<g> = S_v u_0 / L (> 0:
 *                  re-evaporation of falling precipitation); dT_precip = Q_precip / (rho cp_mix), with the
 *                  density the model's rhs_t uses (passed in by the model; without rhs_t's clamp);
 *                  Pflux_<g> = rain + snow + graupel flux, Pflux_NH4SH = the NH4SH settling flux
 *   MCA            dT_mca; dq_<g>_mca = d(vapour + cloud + ice); dNH4SH_mca; mca_freq = share of the MCA calls
 *                  that changed t in the cell (0..1)
 *
 * <g> runs over the condensable gases the model passes in (H2O NH3 CH4, plus H2S on the ice giants). Fields of a
 * stage the model never called in the interval are left out of the file (the ice giants have no NH4SH
 * equilibrium stage). Rates are per PHYSICAL second: the summed change over the interval / (n dt L_atm 1e3 / u_0).
 * Units: K/s; g/m3/s (species x 1e3, the scale of the slices and of printMinMax); P_* kg/m2/s as the arrays
 * hold. Means are over k = 0..km-2 (k = km-1 is the seam copy).
 *
 * Output: <output_path>/<Planet>_proczm_0_<iter>.vtk, legacy ASCII structured grid with the geometry of the
 * zonal slice (DIMENSIONS jm im 1, x = 0.1 i, y = 0.05 j), plus 'height' (km) and 'PressureStat' (zonal mean at
 * the write), so ParaView and tools/render.py read it as a zonal slice. ~2 MB per file.
 *
 * Cost when off: one branch per call. When on: (2 + number of gases) full-field snapshots and one pass per
 * stage. It only READS the model, so the prognostic fields are byte-identical with and without it. The
 * accumulators are not in the restart: a resumed run starts a fresh interval.
 */

#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include "Array.h"

template<class Model>
class ProcessDiag {
public:
    // One condensable gas. Precipitation pointers may be null where the model has no such array.
    struct Gas {
        const char* name;
        Array *vap, *cloud, *ice;
        Array *P_r, *P_s, *P_g;          // precipitation fluxes [kg/m2/s]
        Array *S_v, *S_c, *S_i;          // precipitation tendencies [nondim]
    };
    typedef std::function<double(int, int, int)> RhoFn;

    static int stride(){
        static const int v = [](){
            const std::string e_ = std::string(Model::planet_tag()) + "_PROC_DIAG";
            const char* e = getenv(e_.c_str()); return e ? atoi(e) : 0; }();
        return v;
    }

    ProcessDiag(Model& m, const std::vector<Gas>& gases, RhoFn rho_heat)
        : m(m), gas(gases), rho_heat(rho_heat), on(stride() > 0), n_it(0), n_precip(0), n_mca(0),
          seen_satadj(false), seen_chem(false), seen_precip(false), seen_mca(false) {
        if(!on) return;
        const int G = (int)gas.size();
        // field layout: the order the file is written in
        F_DT_SATADJ = add("dT_satadj");
        for(int g = 0; g < G; g++){ F_COND.push_back(add("cond_" + std::string(gas[g].name))); F_EVAP.push_back(add("evap_" + std::string(gas[g].name))); }
        F_DT_CHEM = add("dT_chem"); F_NH4SH_FORM = add("nh4sh_form"); F_NH4SH_DECOMP = add("nh4sh_decomp");
        F_DT_PRECIP = add("dT_precip");
        for(int g = 0; g < G; g++) F_CONV.push_back(add("conv_" + std::string(gas[g].name)));
        for(int g = 0; g < G; g++) F_REVAP.push_back(add("revap_" + std::string(gas[g].name)));
        for(int g = 0; g < G; g++) F_PFLUX.push_back(add("Pflux_" + std::string(gas[g].name)));
        F_PFLUX_NH4SH = add("Pflux_NH4SH");
        F_DT_MCA = add("dT_mca");
        for(int g = 0; g < G; g++) F_DQ_MCA.push_back(add("dq_" + std::string(gas[g].name) + "_mca"));
        F_DNH4SH_MCA = add("dNH4SH_mca"); F_MCA_FREQ = add("mca_freq");
        n3 = (size_t)m.im * m.jm * m.km;
        acc.assign(names.size() * m.im * m.jm, 0.0);
        snap.assign(2 + G, std::vector<double>(n3, 0.0));
        printf("\n      %s: [PROC] %s_PROC_DIAG=%d: process fields, zonal means every %d iterations -> %s/%s_proczm_0_<iter>.vtk\n",
               Model::planet_tag(), Model::planet_tag(), stride(), stride(), m.output_path.c_str(), Model::planet_name());
    }
    bool active() const { return on; }

    // ---- SatAdj: before the first SaturationAdjustment call / after the last (before the NH4SH equilibrium)
    void satadj_begin(){ if(!on) return; take_T(); for(size_t g = 0; g < gas.size(); g++) take(1 + g, gas[g].cloud, gas[g].ice, nullptr); }
    void satadj_end(){
        if(!on) return; seen_satadj = true;
        diff_T(F_DT_SATADJ, -1);
        for(size_t g = 0; g < gas.size(); g++) diff(1 + g, gas[g].cloud, gas[g].ice, nullptr, F_COND[g], F_EVAP[g]);
    }
    // ---- NH4SH equilibrium (ATJUP / ATSAT)
    void chem_begin(){ if(!on) return; take_T(); take(1 + gas.size(), &m.nh4sh, nullptr, nullptr); }
    void chem_end(){ if(!on) return; seen_chem = true; diff_T(F_DT_CHEM, -1); diff(1 + gas.size(), &m.nh4sh, nullptr, nullptr, F_NH4SH_FORM, F_NH4SH_DECOMP); }
    // ---- Precipitation: sampled after the call
    void precip_begin(){}
    void precip_end(){
        if(!on) return; seen_precip = true;
        const double c = coupling();
        const double cp = m.cp_mix;
        sample(F_DT_PRECIP, [&](int i, int j, int k) -> double {
            const double r = rho_heat(i, j, k); return (r > 0.0 && cp > 0.0) ? c * m.Q_precip.x[i][j][k] / (r * cp) : 0.0; });
        for(size_t g = 0; g < gas.size(); g++){
            const Gas& s = gas[g];
            sample(F_CONV[g],  [&](int i, int j, int k) -> double { return 1e3 * c * (at(s.S_c, i, j, k) + at(s.S_i, i, j, k)) * m.u_0 / (m.L_atm * 1e3); });
            sample(F_REVAP[g], [&](int i, int j, int k) -> double { return 1e3 * c * at(s.S_v, i, j, k) * m.u_0 / (m.L_atm * 1e3); });
            sample(F_PFLUX[g], [&](int i, int j, int k) -> double { return at(s.P_r, i, j, k) + at(s.P_s, i, j, k) + at(s.P_g, i, j, k); });
        }
        sample(F_PFLUX_NH4SH, [&](int i, int j, int k) -> double { return m.P_nh4sh.x[i][j][k]; });
        n_precip++;
    }
    // ---- MCA
    void mca_begin(){
        if(!on) return;
        take_T(); for(size_t g = 0; g < gas.size(); g++) take(1 + g, gas[g].vap, gas[g].cloud, gas[g].ice);
        take(1 + gas.size(), &m.nh4sh, nullptr, nullptr);
    }
    void mca_end(){
        if(!on) return; seen_mca = true;
        diff_T(F_DT_MCA, F_MCA_FREQ);
        for(size_t g = 0; g < gas.size(); g++) diff(1 + g, gas[g].vap, gas[g].cloud, gas[g].ice, F_DQ_MCA[g], F_DQ_MCA[g]);
        diff(1 + gas.size(), &m.nh4sh, nullptr, nullptr, F_DNH4SH_MCA, F_DNH4SH_MCA);
        n_mca++;
    }
    // ---- end of every iteration: count it; write + reset at the stride
    void iteration_done(int iter_n){
        if(!on) return;
        n_it++;
        if(iter_n % stride() == 0 && n_it > 0){ write(iter_n); reset(); }
    }

private:
    Model& m;
    std::vector<Gas> gas;
    RhoFn rho_heat;
    bool on;
    size_t n3;
    int n_it, n_precip, n_mca;
    bool seen_satadj, seen_chem, seen_precip, seen_mca;
    std::vector<std::string> names;
    int F_DT_SATADJ, F_DT_CHEM, F_NH4SH_FORM, F_NH4SH_DECOMP, F_DT_PRECIP, F_PFLUX_NH4SH, F_DT_MCA, F_DNH4SH_MCA, F_MCA_FREQ;
    std::vector<int> F_COND, F_EVAP, F_CONV, F_REVAP, F_PFLUX, F_DQ_MCA;
    std::vector<double> acc;                         // [field][i][j], sums over k and over the calls
    std::vector<std::vector<double> > snap;          // 0 = T, 1..G = per gas, G+1 = NH4SH

    int add(const std::string& n){ names.push_back(n); return (int)names.size() - 1; }
    static double coupling(){
        static const double v = [](){
            const std::string e_ = std::string(Model::planet_tag()) + "_PRECIP_COUPLING";
            const char* e = getenv(e_.c_str()); return e ? atof(e) : 1.0; }();
        return v;
    }
    static double at(const Array* a, int i, int j, int k){ return a ? a->x[i][j][k] : 0.0; }
    size_t id(int i, int j, int k) const { return ((size_t)i * m.jm + j) * m.km + k; }
    double& a(int f, int i, int j){ return acc[((size_t)f * m.im + i) * m.jm + j]; }
    // species in g/m3: 1e3 x (a + b + c)
    static double sp(const Array* a, const Array* b, const Array* c, int i, int j, int k){
        return 1e3 * (at(a, i, j, k) + at(b, i, j, k) + at(c, i, j, k));
    }

    void take_T(){
        std::vector<double>& b = snap[0];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++) for(int k = 0; k < m.km; k++) b[id(i, j, k)] = m.t.x[i][j][k] * m.t_ref;
    }
    void take(size_t s, const Array* a1, const Array* a2, const Array* a3){
        std::vector<double>& b = snap[s];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++) for(int k = 0; k < m.km; k++) b[id(i, j, k)] = sp(a1, a2, a3, i, j, k);
    }
    // each (i, j) is owned by one thread and sums its k in order: deterministic at any thread count
    void diff_T(int f, int f_freq){
        const std::vector<double>& b = snap[0];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double d = 0.0, c = 0.0;
            for(int k = 0; k < m.km - 1; k++){ const double x = m.t.x[i][j][k] * m.t_ref - b[id(i, j, k)]; d += x; if(x != 0.0) c += 1.0; }
            a(f, i, j) += d; if(f_freq >= 0) a(f_freq, i, j) += c;
        }
    }
    // f_pos == f_neg: net change; otherwise the cell-wise sign split
    void diff(size_t s, const Array* a1, const Array* a2, const Array* a3, int f_pos, int f_neg){
        const std::vector<double>& b = snap[s];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double p = 0.0, n = 0.0;
            for(int k = 0; k < m.km - 1; k++){ const double x = sp(a1, a2, a3, i, j, k) - b[id(i, j, k)]; if(f_pos == f_neg) p += x; else if(x > 0.0) p += x; else n += x; }
            a(f_pos, i, j) += p; if(f_pos != f_neg) a(f_neg, i, j) += n;
        }
    }
    void sample(int f, const std::function<double(int, int, int)>& q){
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double d = 0.0;
            for(int k = 0; k < m.km - 1; k++) d += q(i, j, k);
            a(f, i, j) += d;
        }
    }
    void reset(){
        std::fill(acc.begin(), acc.end(), 0.0); n_it = n_precip = n_mca = 0;
        seen_satadj = seen_chem = seen_precip = seen_mca = false;
    }

    // which stage a field belongs to: 0 satadj, 1 chem, 2 precip, 3 mca
    int stage_of(int fi) const {
        if(fi < F_DT_CHEM) return 0;
        if(fi < F_DT_PRECIP) return 1;
        if(fi < F_DT_MCA) return 2;
        return 3;
    }
    void write(int iter_n){
        const double s_it = m.dt * m.L_atm * 1.0e3 / m.u_0;     // physical seconds per iteration
        const double t_int = n_it * s_it, nk = m.km - 1;
        const std::string name = m.output_path + "/" + Model::planet_name() + "_proczm_0_" + std::to_string(iter_n) + ".vtk";
        std::ofstream f(name);
        if(!f.is_open()){ fprintf(stderr, "%s: [PROC] cannot open %s\n", Model::planet_tag(), name.c_str()); return; }
        char buf[64];
        f << "# vtk DataFile Version 3.0\n";
        f << "Process_zonal_mean_" << Model::planet_short() << " iterations " << (iter_n - n_it + 1) << ".." << iter_n << " (" << t_int << " s)\n";
        f << "ASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS " << m.jm << " " << m.im << " 1\nPOINTS " << m.im * m.jm << " float\n";
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){ snprintf(buf, sizeof buf, "%.4f %.4f 0.0\n", 0.1 * i, 0.05 * j); f << buf; }
        f << "POINT_DATA " << m.im * m.jm << "\n";
        auto field = [&](const char* nm, const std::function<double(int, int)>& val){
            f << "SCALARS " << nm << " float 1\nLOOKUP_TABLE default\n";
            for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){ snprintf(buf, sizeof buf, "%.6e\n", val(i, j)); f << buf; }
        };
        field("height", [&](int i, int) -> double { return m.get_layer_height(i); });
        field("PressureStat", [&](int i, int j) -> double { double s = 0.0; for(int k = 0; k < m.km - 1; k++) s += m.p_stat.x[i][j][k]; return s / nk; });
        const bool seen[4] = { seen_satadj, seen_chem, seen_precip, seen_mca };
        std::string left;
        for(int fi = 0; fi < (int)names.size(); fi++){
            const int st = stage_of(fi);
            if(!seen[st]){ if(left.find(names[fi]) == std::string::npos) left += " " + names[fi]; continue; }
            double norm;
            if(st == 2)                norm = n_precip > 0 ? n_precip * nk : 0.0;   // mean over the calls
            else if(fi == F_MCA_FREQ)  norm = n_mca > 0 ? n_mca * nk : 0.0;         // share of MCA calls
            else                       norm = t_int * nk;                           // rate per second
            field(names[fi].c_str(), [&](int i, int j) -> double { return norm > 0.0 ? acc[((size_t)fi * m.im + i) * m.jm + j] / norm : 0.0; });
        }
        f.close();
        printf("      %s: [PROC] wrote %s (%d iterations = %.1f s, %d precipitation calls, %d MCA calls)%s%s\n",
               Model::planet_tag(), name.c_str(), n_it, t_int, n_precip, n_mca, left.empty() ? "" : "; stages not called, left out:", left.c_str());
    }
};

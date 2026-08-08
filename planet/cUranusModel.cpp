/*
 * Uranus Atmosphere Circulation Model (JACM) applied to laminar flow
 * program for the computation of uranus-atmospherical circulating flows in a spherical shell
 * modeling of the atmosphere with gases, their cloud and ice formation: H2O, H2S, NH3 and NH4SH
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/

#include "cUranusModel.h"

#include "ChemistryUran.h"
#include "PressureSolverUran.h"
#include "SaturationAdjustmentUran.h"
#include "BC_Uran.h"
#include "VelocityInitializerUran.h"
#include "ConvectiveAdjustmentUran.h"
#include "TurbulenceUran.h"
#include "RadiationUran.h"
#include "PrecipitationUran.h"
#include "PressureSolver.h"

using namespace std;
using namespace tinyxml2;

// Dry convective adjustment (ConvectiveAdjustmentUran), the SHARED ConvectiveAdjustment<Planet>
// that ATSAT, ATJUP and ATNEPT already run. DEFAULT OFF, so every existing ATURAN run stays
// bit-identical; ATURAN_CONV_ADJ=1 switches it on. It restores any superadiabatic column to the
// dry adiabat while conserving the column's mass-weighted enthalpy, and nothing else in this
// model does that. Whether Uranus develops such columns at all is unmeasured here — switching it
// on and reading the per-iteration report (columns touched, layers mixed, max dT, enthalpy drift)
// is how to find out. Same convention as every other ported module: gated, off, measured later.
//
// Note the adiabat it compares against is g/cp_mix, so it reads the cp_mix corrected in 1992397 —
// 1.0008 K/km rather than the 1.3620 the old mass weighting implied.
// Turbulence closure (TurbulenceUran), the SHARED Turbulence<Planet> that ATSAT, ATJUP and ATNEPT
// run. DEFAULT OFF; ATURAN_TURB=1 switches it on, and turb_model = "none" also disables it, so
// BOTH have to allow it — on ATSAT those were once independent switches and only one of them
// decided anything, which is worth not repeating.
//
// turb_model is a real configuration entry here (uranus/config_aturan.xml, default k_omega_SST,
// with k_omega, k_epsilon and none kept as commented alternatives in param.py), which is how
// ATSAT, ATJUP and ATOM do it; ATNEPT hard-codes the string instead. ATURAN_TURB_MODEL overrides
// the config value at runtime without regenerating anything.
//
// STAGE ONE OF THREE. This fills nue* and the closure diagnostics from the velocity field. It does
// NOT integrate k* and dis* — RungeKuttaUran predates the closure and has no stages for them — and
// nue* reaches no momentum or scalar equation. With the knob ON it moves the diagnostic arrays and
// nothing else; with it OFF it does not run.
//
// ATSAT's warning applies here and has NOT been checked on Uranus: ATSAT's `re` is 1000, which made
// its eddy viscosity ~77x SMALLER than the molecular background — the opposite of ATJUP's situation
// and the opposite of what ATJUP's comment claimed. ATURAN's `re` is 1000 too. The first thing to
// do after switching this on is compare nue* against 1/re.
// Grey multi-layer radiation (RadiationUran), the SHARED Radiation<Planet> that the other three
// already run. DEFAULT OFF, so every existing ATURAN run stays byte-identical; ATURAN_RADIATION=1
// switches it on. It fills the DIAGNOSTIC arrays radiation / epsilon / Q_rad and touches neither t
// nor any rhs — wiring the heating into the temperature equation is a separate step, as it was on
// the other three models.
//
// Uranus is the hard case for this scheme in the opposite direction to Neptune. It absorbs
// S*(1-A) = 3.696*0.700 = 2.59 W/m2 of sunlight against an intrinsic flux of 0.042 W/m2 — 1.6 % of
// the budget, and consistent with zero within its own error bar. Neptune's internal flux is nearly
// a third of its budget and Jupiter's, which this scheme was calibrated on, is about half at forty
// times the absolute flux. A grey scheme tuned there has no claim on a planet in near radiative
// equilibrium, and ATNEPT's own port found it emitting 23x Neptune's entire budget. Expect trouble
// and measure it rather than assuming the port carried.
//
// MEASURED: 2.22x the budget at iteration 1, 2.51x at iteration 2, and NOT because the column is
// clamped. 75e7403's message guessed at the t_max ceiling; that guess is wrong at this run length.
// Peak temperature is 411.8 K against a ceiling of 764 K (t_max 10.0 nondim x t_ref 76.4), so the
// clamp never engages. The 1477 K in the earlier as-is reference run came from the PRE-29693b6
// integrator, which was not RK4; the separated one peaks 1065 K lower. Whatever drives the
// over-emission is in the scheme or in the photosphere temperature at 0.17 bar, not in a clamp.
// Precipitation microphysics (PrecipitationUran), the SHARED Precipitation<Planet>. DEFAULT OFF.
// ATSAT runs its precipitation ON by an explicit decision taken after measuring it; ATURAN's has
// not been measured, and a module that will eventually feed back into the condensate fields is not
// something to switch on by inheritance.
//
// AS PORTED IT IS DIAGNOSTIC ONLY, and that is worth stating because ATNEPT's equivalent commit
// says otherwise. Precipitation.h writes NO condensate array: it fills the P_* fluxes, Q_precip
// and the S_precip_* source terms, and in ATURAN nothing reads S_precip_* at all. ATSAT and ATJUP
// consume them in their RHS; ATNEPT only allocates them, as here. So the coupling is a further
// step, exactly as it was for the turbulence closure.
//
// The call still sits AFTER the SaturationAdjustment block. That is precautionary rather than
// currently load-bearing — once the S_precip_* terms reach the RHS the ordering will matter, and
// putting it in the right place now costs nothing.
static int precip_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_PRECIP"); return e ? atoi(e) : 0; }();
    return v;
}

static int radiation_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_RADIATION"); return e ? atoi(e) : 0; }();
    return v;
}

static int turb_env_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_TURB"); return e ? atoi(e) : 0; }();
    return v;
}

static int conv_adj_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_CONV_ADJ"); return e ? atoi(e) : 0; }();
    return v;
}

// Which pressure solver runs. DEFAULT 0 = ATURAN's own PressureSolverUran, so every existing run
// stays byte-identical; ATURAN_PRESS_SOLVER=1 selects the SHARED PressureSolver<Planet> that
// ATSAT and ATJUP run.
//
// THESE ARE NOT THE SAME ALGORITHM, which is why this is a knob rather than a replacement. The
// shared version carries the red-black ordering that makes a sweep mean one thing on any thread
// count, the obstacle and rigid-lid handling, and the metric-radius and coordinate-stretching
// hooks. What ATURAN's own solver does that the shared one may not is not established here — that
// comparison is a physics question, not a refactor, and the same one ATNEPT left open.
static int press_solver_shared(){
    static const int v = [](){ const char* e = getenv("ATURAN_PRESS_SOLVER"); return e ? atoi(e) : 0; }();
    return v;
}

// README item 1's instrument. ATURAN_TBUDGET=1 turns it on; ATURAN_TBUDGET_J and _K choose the
// column, defaulting to the equator at mid-longitude, which is the column item 1's own table was
// taken from. See the members in cUranusModel.h for why it is one column and one RK stage.
bool cUranusModel::tbudget_enabled(){
    static const bool v = [](){ const char* e = getenv("ATURAN_TBUDGET"); return e && atoi(e) != 0; }();
    return v;
}
int cUranusModel::tbudget_column_j(){
    static const int v = [](){ const char* e = getenv("ATURAN_TBUDGET_J"); return e ? atoi(e) : 90; }();
    return v;
}
int cUranusModel::tbudget_column_k(){
    static const int v = [](){ const char* e = getenv("ATURAN_TBUDGET_K"); return e ? atoi(e) : 180; }();
    return v;
}

bool cUranusModel::tattrib_enabled(){
    static const bool v = [](){ const char* e = getenv("ATURAN_TATTRIB"); return e && atoi(e) != 0; }();
    return v;
}

// Snapshot the column at the top of an iteration. Everything after this is attributed to whichever
// stage ran between two consecutive tattrib() calls.
void cUranusModel::tattribBegin(){
    if(!tattrib_enabled()) return;
    const int jb = tbudget_column_j(), kb = tbudget_column_k();
    tatt_prev.assign(im, 0.0);
    for(int i = 0; i < im; i++) tatt_prev[i] = t.x[i][jb][kb] * t_ref;
}

// Difference the column against the previous snapshot and charge the difference to `label`.
// Keyed by name, not by call order: the iter_n % 2 gate means odd and even iterations execute
// different stage sequences, and an index would silently mis-attribute every odd iteration.
void cUranusModel::tattrib(const char* label){
    if(!tattrib_enabled() || tatt_prev.empty()) return;
    const int jb = tbudget_column_j(), kb = tbudget_column_k();

    int idx = -1;
    for(size_t s = 0; s < tatt_names.size(); s++)
        if(tatt_names[s] == label){ idx = (int)s; break; }
    if(idx < 0){
        tatt_names.push_back(label);
        tatt_dsum.push_back(std::vector<double>(im, 0.0));
        idx = (int)tatt_names.size() - 1;
    }

    for(int i = 0; i < im; i++){
        const double cur = t.x[i][jb][kb] * t_ref;
        tatt_dsum[idx][i] += cur - tatt_prev[i];
        tatt_prev[i]       = cur;
    }
}

// Cumulative kelvin charged to each stage since the run began, at the layers item 1 quotes. The
// column mean is the last column: a stage that only redistributes shows large opposite entries at
// the deep and the top and ~0 there, and a stage that adds or removes heat does not.
void cUranusModel::printTemperatureAttribution(){
    if(!tattrib_enabled() || tatt_names.empty()) return;

    const int probe[] = {0, 5, 20, 30, 38, 39, 40};
    const int np = (int)(sizeof(probe)/sizeof(probe[0]));

    std::cout << std::endl
        << "      ATURAN: who moved T, column j = " << tbudget_column_j()
        << " k = " << tbudget_column_k() << ", cumulative K since iteration 1" << std::endl
        << "      stage                ";
    for(int q = 0; q < np; q++) printf("     i=%-2d ", probe[q]);
    printf("  col.mean\n");

    std::vector<double> tot(im, 0.0);
    for(size_t s = 0; s < tatt_names.size(); s++){
        printf("      %-20s", tatt_names[s].c_str());
        for(int q = 0; q < np; q++) printf(" %9.3f", tatt_dsum[s][probe[q]]);
        double mean = 0.0;
        for(int i = 0; i < im; i++){ mean += tatt_dsum[s][i]; tot[i] += tatt_dsum[s][i]; }
        printf("  %9.4f\n", mean / im);
    }
    printf("      %-20s", "ALL STAGES");
    for(int q = 0; q < np; q++) printf(" %9.3f", tot[probe[q]]);
    double tmean = 0.0;
    for(int i = 0; i < im; i++) tmean += tot[i];
    printf("  %9.4f\n", tmean / im);
}

bool cUranusModel::vbudget_enabled(){
    static const bool v = [](){ const char* e = getenv("ATURAN_VBUDGET"); return e && atoi(e) != 0; }();
    return v;
}
int cUranusModel::vbudget_level_i(){
    static const int v = [](){ const char* e = getenv("ATURAN_VBUDGET_I"); return e ? atoi(e) : 20; }();
    return v;
}
int cUranusModel::vbudget_meridian_k(){
    static const int v = [](){ const char* e = getenv("ATURAN_VBUDGET_K"); return e ? atoi(e) : 180; }();
    return v;
}

// The six terms of rhs_v against LATITUDE, at one level and one meridian. Signs are as they enter
// the sum, so the six add to the total. dphdthe is the only baroclinic term — the horizontal
// gradient of p_hydro — and is identically zero unless ATURAN_HYDRO_SPLIT=1.
void cUranusModel::printMeridionalBudget(){
    if(!vbudget_enabled() || vbud_tot.empty()) return;
    const int ib = vbudget_level_i(), kb = vbudget_meridian_k();

    std::cout << std::endl
        << "      ATURAN: meridional budget, level i = " << ib << " meridian k = " << kb
        << ", RK stage 0, nondimensional dv/dt" << std::endl
        << "        lat        v        dp_dyn      dp_hydro     transport"
        << "     diffusion      Coriolis   centrifugal        dv/dt" << std::endl;

    double a_dp=0, a_dph=0, a_tr=0, a_di=0, a_co=0, a_ce=0, a_to=0;
    for(int j = 1; j < jm-1; j += 10){
        printf("      %5.0f %9.5f %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",
               90.0 - (double)j, v.x[ib][j][kb] * u_0,
               vbud_dp[j], vbud_dph[j], vbud_trans[j],
               vbud_diff[j], vbud_cor[j], vbud_cent[j], vbud_tot[j]);
    }
    for(int j = 1; j < jm-1; j++){
        a_dp += fabs(vbud_dp[j]);   a_dph += fabs(vbud_dph[j]);
        a_tr += fabs(vbud_trans[j]); a_di += fabs(vbud_diff[j]);
        a_co += fabs(vbud_cor[j]);  a_ce  += fabs(vbud_cent[j]);
        a_to += fabs(vbud_tot[j]);
    }
    printf("      %5s %9s %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",
           "SUM|.|", "", a_dp, a_dph, a_tr, a_di, a_co, a_ce, a_to);
    std::cout << "      the baroclinic share is dp_hydro / sum of all six = "
              << (a_dp+a_dph+a_tr+a_di+a_co+a_ce > 0.0
                  ? a_dph / (a_dp+a_dph+a_tr+a_di+a_co+a_ce) : 0.0) << std::endl;
}

// One line per layer: the five terms of rhs_t, converted to the kelvin they contribute to THIS
// iteration. The integrator forms y_{n+1} = y_n + dt/6 * (k1 + 2k2 + 2k3 + k4); these are k1, so
// a term's contribution is dt * term * t_ref kelvin if the four stages agree, which is the right
// scale to read even where they do not. The five columns add to dT/dt by construction — the point
// is not the sum but WHICH term carries it, and with what sign, at each height.
void cUranusModel::printTemperatureBudget(){
    if(!tbudget_enabled() || tbud_tot.empty()) return;

    const int    jb = tbudget_column_j(), kb = tbudget_column_k();
    const double f  = dt * t_ref;                 // nondimensional tendency -> K per iteration

    std::cout << std::endl
        << "      ATURAN: temperature budget, column j = " << jb << " k = " << kb
        << ", RK stage 0, K per iteration" << std::endl
        << "         i     T[K]      pressure     transport     diffusion   thermalmass"
        << "     radiation        dT/dt" << std::endl;

    double s_pres = 0.0, s_trans = 0.0, s_diff = 0.0, s_tmf = 0.0, s_rad = 0.0, s_tot = 0.0;
    for(int i = 0; i < im; i++){
        printf("      %4d %8.2f  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
               i, t.x[i][jb][kb] * t_ref,
               f * tbud_pres[i], f * tbud_trans[i], f * tbud_diff[i],
               f * tbud_tmf[i],  f * tbud_rad[i],   f * tbud_tot[i]);
        s_pres += f * tbud_pres[i];   s_trans += f * tbud_trans[i];
        s_diff += f * tbud_diff[i];   s_tmf   += f * tbud_tmf[i];
        s_rad  += f * tbud_rad[i];    s_tot   += f * tbud_tot[i];
    }
    printf("      %4s %8s  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
           "sum", "", s_pres, s_trans, s_diff, s_tmf, s_rad, s_tot);

    // The column sum is the closure item 1 actually cares about. If the column neither gains nor
    // loses heat while its top warms and its deep cools, this line is ~0 while the per-layer dT/dt
    // above is not — that is redistribution, and it says so directly rather than by inference from
    // two checkpoints 27 iterations apart.
    std::cout << "      column sum of dT/dt = " << s_tot
              << " K/iteration over " << im << " layers; mean per layer = "
              << (im > 0 ? s_tot / im : 0.0) << " K/iteration" << std::endl;
}

using namespace AtomUtils;

cUranusModel* cUranusModel::m_model = NULL;

const double cUranusModel::pi180 = 180.0/M_PI;                          // pi180 = 57.3

const double cUranusModel::the_degree = 1.0;                            // compares to 1° step size laterally
const double cUranusModel::phi_degree = 1.0;                            // compares to 1° step size longitudinally

const double cUranusModel::the0 = 0.0;                                  // North Pole
const double cUranusModel::phi0 = 0.0;                                  // zero meridian in Greenwich

const double cUranusModel::r0 = 1.0;                                    // Uranus radius 72000km

const double cUranusModel::dr = 0.025;                                  // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
const double cUranusModel::dthe = the_degree/pi180; 
const double cUranusModel::dphi = phi_degree/pi180;

cUranusModel::cUranusModel():
    j_ellipse(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    has_welcome_msg_printed(false){
//    if(PythonStream::is_enable()){
//        backup = std::cout.rdbuf();
//        std::cout.rdbuf(&ps);
//    }
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
    m_model = this;
    rad.initArray_1D(im, 0);                                            // radial coordinate direction
    the.initArray_1D(jm, 0);                                            // lateral coordinate direction
    phi.initArray_1D(km, 0);                                            // longitudinal coordinate direction
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
    init_layer_heights();
}

cUranusModel::~cUranusModel(){
    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }
    m_model = NULL;
}

#include "cUranusDefaults.cpp.inc"
/*
*
*/
void cUranusModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    if(err){
        doc.PrintError();
        throw std::invalid_argument("   couldn't load config file inside cUranusModel");
    }
    XMLElement *aturan = doc.FirstChildElement("aturan");
    if(!aturan){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("aturan")->FirstChildElement("common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_uranus = doc.FirstChildElement("aturan")->FirstChildElement("uranus");
    if(!elem_uranus){
        return;
    }
#include "UranusLoadConfig.cpp.inc"
}
/*
*
*/
void cUranusModel::Run(){

    #ifdef _OPENMP
        printf("\n\n   number of processors: %d\n\n", omp_get_num_procs());

    #pragma omp parallel
        {
            printf("   thread %d of %d in \"Desktop-Dell-XPS 8960\"\n", 
                omp_get_thread_num(), omp_get_num_threads());
        }
    #else
        printf("   OpenMP is not supported\n");
    #endif
        printf("   ended\n\n");

    mkdir(output_path.c_str(), 0777);

    m_model = this;

    cout.precision(6);
    cout.setf(ios::fixed);

    if(!has_welcome_msg_printed)
        print_welcome_msg();

    initMsg();

    panorama_cnt = 0;

    resetArrays();

    dt = 0.001;                                                         //  no dimension

    // THE turbulence gate, resolved once. ATURAN_TURB and turb_model must BOTH allow it: on ATSAT
    // these were two independent switches and only one of them decided anything, so a run with
    // turb_model = "none" still ran the closure. Both are consulted here and the answer is a bool
    // the rest of the run reads. ATURAN_TURB_MODEL overrides the configured string first.
    if(const char* tm = getenv("ATURAN_TURB_MODEL")) turb_model = tm;
    turb_active = (turb_env_enabled() != 0) && (turb_model != "none");

    init_layer_heights();
    TropopauseLocation();
    init_tropopause_layers();
    VelocityInitializerUran(*this).compute();                           // construction of zonal initial velocities from measurements

    AtomUtils::damp_wiggles(u, nullptr, true, true, true);
    AtomUtils::damp_wiggles(v, nullptr, true, true, true);
    AtomUtils::damp_wiggles(w, nullptr, true, true, true);

//    goto Printout;

    ChemistryUran(*this).ThermalPropertiesUran();
    gam = g * 1.0e3 / cp_mix;                                           // dry adiabatic lapse rate [K/km]: g[m/s²]*1000 / cp[J/(kg·K)]

    init_temperature();
    init_PressureStatic();

//    goto Printout;

    init_vapour_cloud_ice("CH4", ch4_tropopause, coeff_ch4_A, coeff_ch4_B,
        coeff_ch4_A_i, coeff_ch4_B_i, t_0_ch4, t_00_ch4,
        ep_ch4, r_ch4, m_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4,
        ch4, ch4_cloud, ch4_ice, cloudiness_ch4);

    AtomUtils::damp_wiggles(ch4,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("H2O", h2o_tropopause, coeff_h2o_A, coeff_h2o_B,
        coeff_h2o_A_i, coeff_h2o_B_i, t_0_h2o, t_00_h2o,
        ep_h2o, r_h2o, m_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o,
        h2o, h2o_cloud, h2o_ice, cloudiness_h2o);

    AtomUtils::damp_wiggles(h2o,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("H2S", h2s_tropopause, coeff_h2s_A, coeff_h2s_B,
        coeff_h2s_A_i, coeff_h2s_B_i, t_0_h2s, t_00_h2s,
        ep_h2s, r_h2s, m_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s,
        h2s, h2s_cloud, h2s_ice, cloudiness_h2s);

    AtomUtils::damp_wiggles(h2s,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("NH3", nh3_tropopause, coeff_nh3_A, coeff_nh3_B,
        coeff_nh3_A_i, coeff_nh3_B_i, t_0_nh3, t_00_nh3,
        ep_nh3, r_nh3, m_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3,
        nh3, nh3_cloud, nh3_ice, cloudiness_nh3);

    AtomUtils::damp_wiggles(nh3,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_ice,   nullptr, true, true, true);

//    goto Printout;

    SaturationAdjustmentUran(*this).run("H2O",
        coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
        t_0_h2o, t_00_h2o,
        ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
        h2o, h2o_cloud, h2o_ice);

    AtomUtils::damp_wiggles(h2o,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_ice,   nullptr, true, true, true);

    SaturationAdjustmentUran(*this).run("NH3",
        coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
        nh3, nh3_cloud, nh3_ice);

    AtomUtils::damp_wiggles(nh3,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_ice,   nullptr, true, true, true);

    SaturationAdjustmentUran(*this).run("H2S",
        coeff_h2s_A, coeff_h2s_B, coeff_h2s_A_i, coeff_h2s_B_i,
        t_0_h2s, t_00_h2s,
        ep_h2s, lv_h2s, ls_h2s, cp_h2s, r_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, m_h2s,
        h2s, h2s_cloud, h2s_ice);

    AtomUtils::damp_wiggles(h2s,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2s_ice,   nullptr, true, true, true);

    SaturationAdjustmentUran(*this).run("CH4",
        coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
        ch4, ch4_cloud, ch4_ice);

    AtomUtils::damp_wiggles(ch4,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_ice,   nullptr, true, true, true);


    ChemistryUran(*this).ChemMassRateUran();
    ChemistryUran(*this).FluxLimiterNH4SH();

    AtomUtils::damp_wiggles(massflux_h2s,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(massflux_nh3,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(massflux_nh4sh, nullptr, true, true, true);
    AtomUtils::damp_wiggles(fluxlim_nh4sh,  nullptr, true, true, true);

    ChemistryUran(*this).DiffMassFluxUran();

    AtomUtils::damp_wiggles(difflux_h2s,     nullptr, true, true, true);
    AtomUtils::damp_wiggles(difflux_nh3,     nullptr, true, true, true);
    AtomUtils::damp_wiggles(difflux_nh4sh,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(thermalmassflux, nullptr, true, true, true);

    Forces();
    Latent_Heat();

    init_PressureDynamic();

//    goto Printout;

    BC_Uran(*this).bcRadius();                                          // extrapolation in i-direction alomg grid boundaries
    BC_Uran(*this).bcTheta();                                           // extrapolation in j-direction alomg grid boundaries
    BC_Uran(*this).bcPhi();                                             // extrapolation in k-direction alomg grid boundaries

    restoreVar(1.0);

    // The closure SEEDS k* and dis* algebraically (Turbulence<Planet>::init_fields) before the
    // transport equations start carrying them. It belongs exactly here — after Forces,
    // Latent_Heat, the three boundary passes and the first restoreVar — because init_fields reads
    // the velocity field and the layer geometry, and none of that exists at the top of Run().
    // Putting the call there segfaults; this is ATSAT's slot (cSaturnModel.cpp:482, right after
    // its own pre-loop restoreVar) and ATNEPT's after it, for the same reason.
    if(turb_active) TurbulenceUran(*this).init();

//    goto Printout;



    for(iter_n = 1; iter_n <= nm; iter_n++){

        tattribBegin();   // no-op unless ATURAN_TATTRIB is set

        auto begin = std::chrono::high_resolution_clock::now();

        cout << endl << endl;
        cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        cout << " 3D Uranus iterational process" << endl;
        cout << " present state of the computation " << endl << endl
             << " ======================== iteration number n = " << iter_n << " ========================= " << endl << endl
             << "    max total iteration number nm = " << nm << endl
             << "    checkpoint when to write 3D-panorama = " << checkpoint << endl
             << "    panorama_print = " << panorama_print << endl << endl;

        if(iter_n % 2 == 0){

            if(press_solver_shared()) PressureSolver<cUranusModel>(*this).run();
            else                      PressureSolverUran(*this).run();
            AtomUtils::damp_wiggles(p_dyn, nullptr, true, true, true);
            tattrib("pressure+damp");

            SaturationAdjustmentUran(*this).run("H2O",
                coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
                t_0_h2o, t_00_h2o,
                ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
                C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
                h2o, h2o_cloud, h2o_ice);

            SaturationAdjustmentUran(*this).run("H2S",
                coeff_h2s_A, coeff_h2s_B, coeff_h2s_A_i, coeff_h2s_B_i,
                t_0_h2s, t_00_h2s,
                ep_h2s, lv_h2s, ls_h2s, cp_h2s, r_h2s,
                C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, m_h2s,
                h2s, h2s_cloud, h2s_ice);

            SaturationAdjustmentUran(*this).run("NH3",
                coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
                t_0_nh3, t_00_nh3,
                ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
                C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
                nh3, nh3_cloud, nh3_ice);

            SaturationAdjustmentUran(*this).run("CH4",
                coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
                t_0_ch4, t_00_ch4,
                ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
                C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
                ch4, ch4_cloud, ch4_ice);

            tattrib("SaturationAdjust");

            // Must follow the saturation adjustments: it removes condensate in place.
            if(precip_enabled()) PrecipitationUran(*this).run();
            tattrib("Precipitation");

            ChemistryUran(*this).DiffMassFluxUran();

            AtomUtils::damp_wiggles(difflux_h2s,     nullptr, true, true, true);
            AtomUtils::damp_wiggles(difflux_nh3,     nullptr, true, true, true);
            AtomUtils::damp_wiggles(difflux_nh4sh,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(thermalmassflux, nullptr, true, true, true);

            ChemistryUran(*this).ChemMassRateUran();
            ChemistryUran(*this).FluxLimiterNH4SH();

            AtomUtils::damp_wiggles(massflux_h2s,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(massflux_nh3,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(massflux_nh4sh, nullptr, true, true, true);
            AtomUtils::damp_wiggles(fluxlim_nh4sh,  nullptr, true, true, true);

            tattrib("Chemistry");

            Forces();
            tattrib("Forces");
            Latent_Heat();
            tattrib("Latent_Heat");

        }  // if loop

        RungeKuttaUran();
        tattrib("RungeKutta");

        // ===== THE VERTICAL PASS OF THE t FILTER HAS ITS OWN STRENGTH =====
        //
        // aa6f739 measured this call supplying 98 % of the photosphere's warming with a column
        // mean of -0.012 K: a Shapiro 1-2-1 filter for grid-scale noise acting as the model's
        // dominant vertical heat transport. ATURAN_DAMP_T_VERT scales the VERTICAL pass only; the
        // horizontal passes are untouched.
        //
        // DEFAULT 1.0, AND THAT IS BIT-IDENTICAL. damp_wiggles applies its axes in the order
        // k, j, i, so splitting the single (true,true,true) call into horizontal-then-vertical
        // preserves the order exactly. Verified rather than argued: 7 of 7 output files identical
        // with the variable unset, and the sweep's own 1.0 row reproduces README item 8's
        // corrected-metric numbers, OLR/in 29.58 and T(tau=1) 135.98.
        //
        // 0.0 disables the vertical pass entirely and is the diagnostic extreme, not a proposal.
        // README item 1 carries the sweep and what it costs.
        static const double damp_t_vert = [](){
            const char* e = getenv("ATURAN_DAMP_T_VERT"); return e ? atof(e) : 1.0; }();
        AtomUtils::damp_wiggles(t, nullptr, false, true, true);              // horizontal, unchanged
        if(damp_t_vert != 0.0)
            AtomUtils::damp_wiggles(t, nullptr, true, false, false, damp_t_vert);
        tattrib("damp_wiggles(t)");
        AtomUtils::damp_wiggles(u, nullptr, true, true, true);
        AtomUtils::damp_wiggles(v, nullptr, true, true, true);
        AtomUtils::damp_wiggles(w, nullptr, true, true, true);

        BC_Uran(*this).bcRadius();                                      // extrapolation in i-direction along grid boundaries
        BC_Uran(*this).bcTheta();                                       // extrapolation in j-direction along grid boundaries
        BC_Uran(*this).bcPhi();                                         // extrapolation in k-direction along grid boundaries
        tattrib("BoundaryConds");

        // How far the run is from a steady state, and WHERE. MUST run BEFORE restoreVar: it
        // differences each field against the n-copy restoreVar is about to overwrite, so after
        // that call every number it prints is identically zero. That is the reason it was never
        // called from anywhere useful, and the reason it is called here.
        //
        // Same cadence as printMinMax. ATURAN_STEADY=0 switches it off. This ADDS LINES TO THE
        // LOG that ATURAN has never printed — the routine has never run in this model — but it
        // writes nothing and changes no output file.
        static const int steady_on = [](){
            const char* e = getenv("ATURAN_STEADY"); return e ? atoi(e) : 1; }();
        if(steady_on && iter_n % checkpoint == 0) steadyQuery();

        restoreVar(1.0);
        tattrib("restoreVar");

        // After the state has been advanced and the boundaries applied: put any
        // superadiabatic column back on the dry adiabat. Off by default (ATURAN_CONV_ADJ).
        if(radiation_enabled()) RadiationUran(*this).run();
        if(turb_active) TurbulenceUran(*this).run();
        if(conv_adj_enabled()) ConvectiveAdjustmentUran(*this).run();

        tattrib("Rad/Turb/ConvAdj");

        panorama_cnt++;

        if(iter_n % checkpoint == 0){
            printMinMax();
            printMeridionalBudget();         // no-op unless ATURAN_VBUDGET is set
            printTemperatureBudget();        // no-op unless ATURAN_TBUDGET is set
            printTemperatureAttribution();   // no-op unless ATURAN_TATTRIB is set
            writeData();
        }

        if(panorama_cnt == panorama_print) panorama_cnt = 0;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for one time step\n", elapsed.count() * 1e-9);

    }  // end for iter_n

    cout << endl << "      Uranus: run_3D_loop atm ended ..........................." << endl;


/*
    Printout:
        printMinMax();
        writeData();
*/


    print_final_msg();

    return;
}
/*
*
*/
void cUranusModel::resetArrays(){
    cout << endl << "      ATURAN: resetArrays" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    Topography.initArray_2D(jm, km, 0.0); // topography
    LatentHeat.initArray_2D(jm, km, 0.0);            // areas of higher latent heat
    Precipitation.initArray_2D(jm, km, 0.0);         // areas of higher precipitation
    precipitable_water.initArray_2D(jm, km, 0.0);    // areas of precipitable water in the air
    nh3_total.initArray_2D(jm, km, 0.0);             // areas of higher nh3 concentration
    nh3_cloud_total.initArray_2D(jm, km, 0.0);       // areas of higher nh3_cloud concentration
    nh3_ice_total.initArray_2D(jm, km, 0.0);         // areas of higher nh3_ice concentration
    aux_2D_v.initArray_2D(jm, km, 0.0);              // auxilliar field v
    aux_2D_w.initArray_2D(jm, km, 0.0);              // auxilliar field w
    tropopause_height.initArray_2D(jm, km, 0.0); // local height of the tropopause

    t.initArray(im, jm, km, ta);                    // temperature
    u.initArray(im, jm, km, ua);                    // u-component velocity component in r-direction
    v.initArray(im, jm, km, va);                    // v-component velocity component in theta-direction
    w.initArray(im, jm, km, wa);                    // w-component velocity component in phi-direction

    h2o.initArray(im, jm, km, 0.0);                    // water vapour
    h2o_cloud.initArray(im, jm, km, 0.0);                // cloud water
    h2o_ice.initArray(im, jm, km, 0.0);                  // cloud ice
    ch4.initArray(im, jm, km, 0.0);                    // methane vapour
    ch4_cloud.initArray(im, jm, km, 0.0);                // methane cloud
    ch4_ice.initArray(im, jm, km, 0.0);                  // methane ice
    h2s.initArray(im, jm, km, 0.0);                    // water vapour
    h2s_cloud.initArray(im, jm, km, 0.0);                // cloud water
    h2s_ice.initArray(im, jm, km, 0.0);                  // cloud ice
    nh3.initArray(im, jm, km, 0.0);                 // nh3-vapour
    nh3_cloud.initArray(im, jm, km, 0.0);           // nh3-cloud
    nh3_ice.initArray(im, jm, km, 0.0);             // nh3-ice
    nh4sh.initArray(im, jm, km, 0.0);                 // nh4sh-vapour

    tn.initArray(im, jm, km, ta);                    // temperature new
    un.initArray(im, jm, km, ua);                    // u-velocity component in r-direction new
    vn.initArray(im, jm, km, va);                    // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, wa);                    // w-velocity component in phi-direction new

    ch4n.initArray(im, jm, km, 0.0);                    // methane vapour new
    ch4_cloudn.initArray(im, jm, km, 0.0);                // methane cloud new
    ch4_icen.initArray(im, jm, km, 0.0);                    // methane ice new
    h2on.initArray(im, jm, km, 0.0);                    // water vapour new
    h2o_cloudn.initArray(im, jm, km, 0.0);                // cloud water new
    h2o_icen.initArray(im, jm, km, 0.0);                    // cloud ice new
    h2sn.initArray(im, jm, km, 0.0);                    // water vapour new
    h2s_cloudn.initArray(im, jm, km, 0.0);                // cloud water new
    h2s_icen.initArray(im, jm, km, 0.0);                    // cloud ice new
    nh3n.initArray(im, jm, km, 0.0);                // nh3 new
    nh3_cloudn.initArray(im, jm, km, 0.0);            // nh3_cloud new
    nh3_icen.initArray(im, jm, km, 0.0);            // nh3_ice new
    nh4shn.initArray(im, jm, km, 0.0);                // nh4sh new

    massflux_h2s.initArray(im, jm, km, 0.0);    // mass flux h2s
    massflux_nh3.initArray(im, jm, km, 0.0);    // mass flux nh3
    massflux_nh4sh.initArray(im, jm, km, 0.0);  // mass flux nh4sh
    fluxlim_nh4sh.initArray(im, jm, km, 0.0);   // TVD flux-limiter correction

    // Item 1's budget instrument. Allocated only when it is switched on, so the default build
    // carries six empty vectors and the recording branch in RHSUran is a size check that fails.
    if(vbudget_enabled()){
        vbud_dp.assign(jm, 0.0);   vbud_dph.assign(jm, 0.0);  vbud_trans.assign(jm, 0.0);
        vbud_diff.assign(jm, 0.0); vbud_cor.assign(jm, 0.0);  vbud_cent.assign(jm, 0.0);
        vbud_tot.assign(jm, 0.0);
    }
    if(tbudget_enabled()){
        tbud_pres.assign(im, 0.0);  tbud_trans.assign(im, 0.0);
        tbud_diff.assign(im, 0.0);  tbud_tmf.assign(im, 0.0);
        tbud_rad.assign(im, 0.0);   tbud_tot.assign(im, 0.0);
    }

    difflux_h2s.initArray(im, jm, km, 0.0);   // diffusive flux h2s
    difflux_nh3.initArray(im, jm, km, 0.0);   // diffusive flux nh3
    difflux_nh4sh.initArray(im, jm, km, 0.0); // diffusive flux nh4sh

    thermalmassflux.initArray(im, jm, km, 0.0);   // thermal massflux_h2s

    p_dyn.initArray(im, jm, km, pa);                // dynamic pressure
    p_hydro.initArray(im, jm, km, 0.0);             // hydrostatic pressure perturbation
    p_dynn.initArray(im, jm, km, pa);               // dynamic pressure, previous iteration

    // Turbulence closure fields (stage one). Allocated unconditionally so the arrays exist for
    // printMinMax and for the shared BC lists whether or not the closure runs; they stay zero
    // while turb_active is false.
    // RK4 stage accumulators for the separated integrator.
    acc_t.initArray(im, jm, km, 0.0);
    acc_u.initArray(im, jm, km, 0.0);
    acc_v.initArray(im, jm, km, 0.0);
    acc_w.initArray(im, jm, km, 0.0);
    acc_ch4.initArray(im, jm, km, 0.0);
    acc_ch4_cloud.initArray(im, jm, km, 0.0);
    acc_ch4_ice.initArray(im, jm, km, 0.0);
    acc_h2o.initArray(im, jm, km, 0.0);
    acc_h2o_cloud.initArray(im, jm, km, 0.0);
    acc_h2o_ice.initArray(im, jm, km, 0.0);
    acc_h2s.initArray(im, jm, km, 0.0);
    acc_h2s_cloud.initArray(im, jm, km, 0.0);
    acc_h2s_ice.initArray(im, jm, km, 0.0);
    acc_nh3.initArray(im, jm, km, 0.0);
    acc_nh3_cloud.initArray(im, jm, km, 0.0);
    acc_nh3_ice.initArray(im, jm, km, 0.0);
    acc_nh4sh.initArray(im, jm, km, 0.0);
    acc_tke.initArray(im, jm, km, 0.0);
    acc_dis.initArray(im, jm, km, 0.0);

    P_rain.initArray(im, jm, km, 0.0);
    P_snow.initArray(im, jm, km, 0.0);
    P_graupel.initArray(im, jm, km, 0.0);
    P_nh3_rain.initArray(im, jm, km, 0.0);
    P_nh3_snow.initArray(im, jm, km, 0.0);
    P_nh3_graupel.initArray(im, jm, km, 0.0);
    P_ch4_rain.initArray(im, jm, km, 0.0);
    P_ch4_snow.initArray(im, jm, km, 0.0);
    P_ch4_graupel.initArray(im, jm, km, 0.0);
    P_nh4sh.initArray(im, jm, km, 0.0);
    Q_precip.initArray(im, jm, km, 0.0);
    S_precip_h2o.initArray(im, jm, km, 0.0);
    S_precip_h2o_cloud.initArray(im, jm, km, 0.0);
    S_precip_h2o_ice.initArray(im, jm, km, 0.0);
    S_precip_nh3.initArray(im, jm, km, 0.0);
    S_precip_nh3_cloud.initArray(im, jm, km, 0.0);
    S_precip_nh3_ice.initArray(im, jm, km, 0.0);
    S_precip_ch4.initArray(im, jm, km, 0.0);
    S_precip_ch4_cloud.initArray(im, jm, km, 0.0);
    S_precip_ch4_ice.initArray(im, jm, km, 0.0);
    precip_srf_h2o.initArray_2D(jm, km, 0.0);
    precip_srf_nh3.initArray_2D(jm, km, 0.0);
    precip_srf_ch4.initArray_2D(jm, km, 0.0);
    precip_srf_nh4sh.initArray_2D(jm, km, 0.0);
    precip_srf_total.initArray_2D(jm, km, 0.0);

    radiation.initArray(im, jm, km, 0.0);            // net thermal radiative flux [W/m2]
    epsilon.initArray(im, jm, km, 0.0);              // layer emissivity
    Q_rad.initArray(im, jm, km, 0.0);                // radiative heating rate [W/m3]

    tke.initArray(im, jm, km, 0.0);
    dis.initArray(im, jm, km, 0.0);
    tken.initArray(im, jm, km, 0.0);
    disn.initArray(im, jm, km, 0.0);
    nue.initArray(im, jm, km, 0.0);
    nue_t.initArray(im, jm, km, 0.0);
    prod.initArray(im, jm, km, 0.0);
    tke_source.initArray(im, jm, km, 0.0);
    dis_source.initArray(im, jm, km, 0.0);
    vel_star.initArray_2D(jm, km, 0.0);
    p_stat.initArray(im, jm, km, pa);                // static pressure
    rho_mix.initArray(im, jm, km, 0.0);             // local mixture density

    rhs_t.initArray(im, jm, km, 0.0);                // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.0);                // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.0);                // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.0);                // auxilliar field RHS w-velocity component

    rhs_ch4.initArray(im, jm, km, 0.0);                // auxilliar field RHS methane vapour
    rhs_ch4_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS methane cloud
    rhs_ch4_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS methane ice
    rhs_h2o.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_h2o_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS cloud water
    rhs_h2o_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS cloud ice
    rhs_h2s.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_h2s_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS cloud water
    rhs_h2s_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS cloud ice
    rhs_nh3.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh3
    rhs_nh3_cloud.initArray(im, jm, km, 0.0);        // auxilliar field RHS nh3_cloud
    rhs_nh3_ice.initArray(im, jm, km, 0.0);            // auxilliar field RHS nh3_ice
    rhs_nh4sh.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh4sh
    rhs_tke.initArray(im, jm, km, 0.0);
    rhs_dis.initArray(im, jm, km, 0.0);

    aux.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_u.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.0);                // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.0);                // auxilliar field w-velocity component

    Q_Latent.initArray(im, jm, km, 0.0);                // latent heat
    Q_Sensible.initArray(im, jm, km, 0.0);            // sensible heat
    CoriolisForce.initArray(im, jm, km, 0.0);        // coriolis force
    BuoyancyForce.initArray(im, jm, km, 0.0);        // buoyancy force, Boussinesque approximation
    PresGradForce.initArray(im, jm, km, 0.0);// pressure gradient force
    CentrifugalForce.initArray(im, jm, km, 0.0);             // centrifugal force
    SeaMount.initArray(im, jm, km, 0.0);             // sea mount contour

    cloudiness_ch4.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_h2o.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_h2s.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_nh3.initArray(im, jm, km, 0.0); // cloudiness, N in literature

    w_nh3.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh3
    w_h2s.initArray(im, jm, km, 0.0);                // chemical reaction rate of h2s
    w_nh4sh.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh4sh

    j_nh3.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of nh3
    j_h2s.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of h2s
    j_nh4sh.initArray(im, jm, km, 0.0);              // ordinary-diffusion mass flux of nh4sh

    jT_nh3.initArray(im, jm, km, 0.0);               // thermo-diffusion mass flux of nh3
    jT_h2s.initArray(im, jm, km, 0.0);               // thermo-diffusion mass flux of h2s
    jT_nh4sh.initArray(im, jm, km, 0.0);             // thermo-diffusion mass flux of nh4sh

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for resetArrays\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: resetArrays ended" << endl;
    return;
}
/*
*
*/
void cUranusModel::restoreVar(double coeff){
//    cout << endl << "      ATURAN: restoreVar" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                // k* and dis* take the same start-of-step copies the other prognostic fields
                // get; without them the RK4 stages would integrate from a moving base. With the
                // closure off both fields are identically zero, so these are zero too.
                tken.x[i][j][k] = coeff * tke.x[i][j][k];
                disn.x[i][j][k] = coeff * dis.x[i][j][k];
                // p_dynn is the previous-iteration dynamic pressure. Nothing maintained it
                // and nothing allocated it, while steadyQuery's pressure case sat commented
                // out because of that. It belongs with the other n-copies.
                p_dynn.x[i][j][k] = coeff * p_dyn.x[i][j][k];
                un.x[i][j][k] = coeff * u.x[i][j][k];
                vn.x[i][j][k] = coeff * v.x[i][j][k];
                wn.x[i][j][k] = coeff * w.x[i][j][k];
                ch4n.x[i][j][k] = coeff * ch4.x[i][j][k];
                ch4_cloudn.x[i][j][k] = coeff * ch4_cloud.x[i][j][k];
                ch4_icen.x[i][j][k] = coeff * ch4_ice.x[i][j][k];
                h2on.x[i][j][k] = coeff * h2o.x[i][j][k];
                h2o_cloudn.x[i][j][k] = coeff * h2o_cloud.x[i][j][k];
                h2o_icen.x[i][j][k] = coeff * h2o_ice.x[i][j][k];
                h2sn.x[i][j][k] = coeff * h2s.x[i][j][k];
                h2s_cloudn.x[i][j][k] = coeff * h2s_cloud.x[i][j][k];
                h2s_icen.x[i][j][k] = coeff * h2s_ice.x[i][j][k];
                nh3n.x[i][j][k] = coeff * nh3.x[i][j][k];
                nh3_cloudn.x[i][j][k] = coeff * nh3_cloud.x[i][j][k];
                nh3_icen.x[i][j][k] = coeff * nh3_ice.x[i][j][k];
                nh4shn.x[i][j][k] = coeff * nh4sh.x[i][j][k];
//                nh4sh_cloudn.x[i][j][k] = coeff * nh4sh_cloud.x[i][j][k];
            }
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for restoreVar\n", elapsed.count() * 1e-9);

//    cout << "      ATURAN: restoreVar ended" << endl;

    return;
}
/*
*
*/
/*
    cout << endl << "      ATURAN: thermodynamic fft_gaussian_filter_3d in Run begin ......................." << endl;

    fft_gaussian_filter_3d(nh3,1);
    fft_gaussian_filter_3d(nh3_cloud,1);
    fft_gaussian_filter_3d(nh3_ice,1);

    fft_gaussian_filter_3d(nh4sh,1);
//    fft_gaussian_filter_3d(nh4sh_cloud,1);
//    fft_gaussian_filter_3d(nh4sh_ice,1);

    fft_gaussian_filter_3d(h2s,1);
    fft_gaussian_filter_3d(h2s_cloud,1);
    fft_gaussian_filter_3d(h2s_ice,1);

    fft_gaussian_filter_3d(h2o,1);
    fft_gaussian_filter_3d(h2o_cloud,1);
    fft_gaussian_filter_3d(h2o_ice,1);

    cout << endl << "      ATURAN: thermodynamic fft_gaussian_filter_3d in Run end ......................." << endl;
*/



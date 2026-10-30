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
#include "TbUran.h"
#include "CwbUran.h"        // ATURAN_CWB_DIAG, print-only
#include "UbUran.h"          // ATURAN_UB_DIAG, print-only

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
#include "ChromUran.h"
#include "ZseedUran.h"          // ATURAN_T_SEED / ATURAN_ZDEV (row ICE-FROZEN-FLOW)         // ATURAN_CHROM, the chromophore tracer

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
// already run. DEFAULT ON since 2026-10-11 on the user's word; ATURAN_RADIATION=0 restores the run
// without it (not byte-checked). It fills the DIAGNOSTIC arrays radiation / epsilon / Q_rad and touches neither t
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
// DEFAULT 1 SINCE 2026-10-12 on the user's word ("on Uranus exists no rain of any species"; the
// same on Neptune): the scheme runs, with the rate scale calibrated to the emitted flux
// (Planet::precip_rate_scale()), the corrected t_00 values of NH3, CH4 and H2S (cUranusModel.h)
// and the coupling to the species and temperature equations (RHS_Uran_Turb.cpp,
// ATURAN_PRECIP_COUPLING). The notes above describe the scheme as it stood before that.
// ATURAN_PRECIP=0 with ATURAN_T00_LEGACY=1 and ATURAN_LV_H2S=3.534e6 is the model before (byte-identical
// except for the all-zero Precip_* block in the horizontal VTK).
static int precip_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_PRECIP"); return e ? atoi(e) : 1; }();
    return v;
}

static int radiation_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_RADIATION"); return e ? atoi(e) : 1; }();
    return v;
}

// ATURAN_TURB: -1 = not set. turb_model IS the switch since 2026-10-11 (the user's word), as on ATJUP
// since 2026-10-09: the closure runs unless turb_model (or ATURAN_TURB_MODEL, which overrides it) is
// "none" / "laminar" / empty, and param.py's default is k_omega_SST, so a stock run has it ON.
// Until then the closure needed ATURAN_TURB=1 as well and the configured model alone did nothing.
// ATURAN_TURB=0 still switches it off; any other value leaves the decision to turb_model.
static int turb_env_enabled(){
    static const int v = [](){ const char* e = getenv("ATURAN_TURB"); return e ? atoi(e) : -1; }();
    return v;
}

// ATURAN_MCA=<N> (2026-10-26; DEFAULT 1 since 2026-10-27, user: "port the MCA default to the giants"; 0 = off): the moist convective adjustment with the Ledoux criterion
// (MoistConvectiveAdjustmentUran, shared ConvectiveAdjustment.h, where its note and knobs are) every N
// iterations, right after the dry one. Ported from ATJUP. It mass-weights with rho_at() and reads mass
// fractions as c / rho_at(), so it means what it says only with ATURAN_LOCAL_RHO=1 (with the constant r_mix
// every layer weighs the same). No tracer carries a formation enthalpy here (the NH4SH chemistry
// releases none).
static int mca_stride(){
    static const int v = [](){ const char* e = getenv("ATURAN_MCA"); return e ? std::max(0, atoi(e)) : 1; }();
    return v;
}

// ATURAN_SVP_FIT: assign the shared sets (ATPhys.h). The H2O liquid, CH4 and H2S sets are not touched.
void cUranusModel::setSvpFit(){
    const ATPhys::SvpSet hi = ATPhys::svp_fit_h2o_ice(), ni = ATPhys::svp_fit_nh3_ice(), nl = ATPhys::svp_fit_nh3_liq();
    C_h2o_ice = hi.C; L0_h2o_ice = hi.L0; del_alf_h2o_ice = hi.del_alf; del_bet_h2o_ice = hi.del_bet;
    C_nh3_ice = ni.C; L0_nh3_ice = ni.L0; del_alf_nh3_ice = ni.del_alf; del_bet_nh3_ice = ni.del_bet;
    C_nh3     = nl.C; L0_nh3     = nl.L0; del_alf_nh3     = nl.del_alf; del_bet_nh3     = nl.del_bet;
    printf("      ATURAN: ATURAN_SVP_FIT: NH3 liquid / NH3 ice / H2O ice vapour pressures from the shared fits (ATPhys.h)\n");
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
    // every 10th row, plus the five next to each pole (2026-10-16: the residuum sits at 89 S)
    for(int j = 1; j < jm-1; j++){
        if(!(j % 10 == 1 || j <= 5 || j >= jm - 6)) continue;
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
// ATURAN_INIT_SPECIES_LOCAL_RHO: see the note at rho_init() in cUranusModel.h.
void cUranusModel::fill_rho_mix_init(){
    if(!init_species_local_rho()) return;
    double rmin = 1.0e30, rmax = 0.0;
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++){
                const double r = rho_init(i, j, k);
                rho_mix.x[i][j][k] = r;
                rmin = std::min(rmin, r); rmax = std::max(rmax, r);
            }
    printf("      %s: INIT_SPECIES_LOCAL_RHO: R_mix = %.6g, r_mix = %.6g, set-up density %.6g .. %.6g kg/m3\n",
           "ATURAN", R_mix, r_mix, rmin, rmax);
}

void cUranusModel::damp_species_init(Array& f){
    if(!init_species_local_rho()){ AtomUtils::damp_wiggles(f, nullptr, true, true, true); return; }
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++) f.x[i][j][k] /= rho_init(i, j, k);
    AtomUtils::damp_wiggles(f, nullptr, true, true, true);
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++) f.x[i][j][k] *= rho_init(i, j, k);
}

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

    // ATURAN_BC_PHI_PERIODIC (2026-10-13; DEFAULT 1 on the user's word, 0 restores byte-identically;
    // ATPhys.h has the description): the seam
    // meridian k = 0 = km-1 is integrated instead of being reconstructed from its neighbours.
    // The knob is read by every loop that differences in phi; the Shapiro filter lives in lib/
    // and is told here.
    AtomUtils::set_phi_seam_periodic(ATPhys::phi_periodic<cUranusModel>());
    if(ATPhys::phi_periodic<cUranusModel>())
        printf("      ATURAN: the phi seam is periodic: k = 0 is integrated, k = km-1 is its copy (ATURAN_BC_PHI_PERIODIC=0 restores the reconstructed seam)\n");

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
    turb_active = (turb_env_enabled() != 0) && !(turb_model.empty() || turb_model == "none" || turb_model == "laminar");

    init_layer_heights();
    TropopauseLocation();
    init_tropopause_layers();
    VelocityInitializerUran(*this).compute();                           // construction of zonal initial velocities from measurements

    AtomUtils::damp_wiggles(u, nullptr, true, true, true);
    AtomUtils::damp_wiggles(v, nullptr, true, true, true);
    AtomUtils::damp_wiggles(w, nullptr, true, true, true);

//    goto Printout;

    if(svp_fit()) setSvpFit();                       // ATURAN_SVP_FIT, before any SVP is evaluated
    ChemistryUran(*this).ThermalPropertiesUran();
    gam = g * 1.0e3 / cp_mix;                                           // dry adiabatic lapse rate [K/km]: g[m/s²]*1000 / cp[J/(kg·K)]

    // ATURAN_T_INIT_ANCHOR=1 (cUranusModel.h has the note): the reference densities follow the
    // anchored bottom state, one factor on all of them, so that r_mix * R_mix * T_bottom is the
    // bottom pressure. cp_mix, R_mix, mue_mix, k_mix and gam are ratios and stay.
    if(t_init_anchor() != 0){
        const double T_bottom = t_bottom_anchor_K();
        const double f = p_bottom_bar() * 1.0e5 / (R_mix * T_bottom) / r_mix;
        r_h2 *= f; r_he *= f; r_ch4 *= f; r_nh3 *= f; r_h2s *= f; r_h2o *= f; r_nh4sh *= f;
        c_h2 *= f; c_he *= f; c_ch4 *= f; c_nh3 *= f; c_h2s *= f; c_h2o *= f; c_nh4sh *= f;
        r_mix *= f; c_mix *= f;
        const double n = g / (gam * R_ref);
        printf("      ATURAN: T_INIT_ANCHOR: adiabat through %.1f K at 1 bar, n = %.4f; bottom %.2f K at %.2f bar;"
               " reference densities x %.4f (r_mix %.4f kg/m3); L_atm %.1f km\n",
               t_1bar_K(), n, T_bottom, p_bottom_bar(), f, r_mix, L_atm);
        for(int e = 0; e < 2; e++){
            const double T_tr = e == 0 ? t_equator : t_pole;
            const double h_tr = (T_bottom - T_tr) / gam;
            const double p_tr = p_bottom_bar() * std::pow(T_tr / T_bottom, n);
            const double p_lid = L_atm > h_tr ? p_tr * std::exp(-(L_atm - h_tr) * g / (R_ref * T_tr))
                                              : p_bottom_bar() * std::pow((T_bottom - gam * L_atm) / T_bottom, n);
            printf("      ATURAN: T_INIT_ANCHOR: %s: isothermal at %.1f K above %.1f km (%.4f bar); lid at %.4e bar\n",
                   e == 0 ? "equator" : "pole   ", T_tr, h_tr, p_tr, p_lid);
        }
    }

    init_temperature();
    balanceThermalWind();                                               // ATURAN_TW_BALANCE (on with ATURAN_GEO_BALANCE); needs t and the prescribed w
    init_PressureStatic();
    fill_rho_mix_init();

//    goto Printout;

    init_vapour_cloud_ice("CH4", ch4_tropopause, coeff_ch4_A, coeff_ch4_B,
        coeff_ch4_A_i, coeff_ch4_B_i, t_0_ch4, t_00_ch4,
        ep_ch4, r_ch4, m_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4,
        ch4, ch4_cloud, ch4_ice, cloudiness_ch4);

    damp_species_init(ch4);
    damp_species_init(ch4_cloud);
    damp_species_init(ch4_ice);

    init_vapour_cloud_ice("H2O", h2o_tropopause, coeff_h2o_A, coeff_h2o_B,
        coeff_h2o_A_i, coeff_h2o_B_i, t_0_h2o, t_00_h2o,
        ep_h2o, r_h2o, m_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o,
        h2o, h2o_cloud, h2o_ice, cloudiness_h2o);

    damp_species_init(h2o);
    damp_species_init(h2o_cloud);
    damp_species_init(h2o_ice);

    init_vapour_cloud_ice("H2S", h2s_tropopause, coeff_h2s_A, coeff_h2s_B,
        coeff_h2s_A_i, coeff_h2s_B_i, t_0_h2s, t_00_h2s,
        ep_h2s, r_h2s, m_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s,
        h2s, h2s_cloud, h2s_ice, cloudiness_h2s);

    damp_species_init(h2s);
    damp_species_init(h2s_cloud);
    damp_species_init(h2s_ice);

    init_vapour_cloud_ice("NH3", nh3_tropopause, coeff_nh3_A, coeff_nh3_B,
        coeff_nh3_A_i, coeff_nh3_B_i, t_0_nh3, t_00_nh3,
        ep_nh3, r_nh3, m_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3,
        nh3, nh3_cloud, nh3_ice, cloudiness_nh3);

    damp_species_init(nh3);
    damp_species_init(nh3_cloud);
    damp_species_init(nh3_ice);

//    goto Printout;

    SaturationAdjustmentUran(*this).run("H2O",
        coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
        t_0_h2o, t_00_h2o,
        ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
        h2o, h2o_cloud, h2o_ice);

    damp_species_init(h2o);
    damp_species_init(h2o_cloud);
    damp_species_init(h2o_ice);

    SaturationAdjustmentUran(*this).run("NH3",
        coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
        nh3, nh3_cloud, nh3_ice);

    damp_species_init(nh3);
    damp_species_init(nh3_cloud);
    damp_species_init(nh3_ice);

    SaturationAdjustmentUran(*this).run("H2S",
        coeff_h2s_A, coeff_h2s_B, coeff_h2s_A_i, coeff_h2s_B_i,
        t_0_h2s, t_00_h2s,
        ep_h2s, lv_h2s, ls_h2s, cp_h2s, r_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, m_h2s,
        h2s, h2s_cloud, h2s_ice);

    damp_species_init(h2s);
    damp_species_init(h2s_cloud);
    damp_species_init(h2s_ice);

    SaturationAdjustmentUran(*this).run("CH4",
        coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
        ch4, ch4_cloud, ch4_ice);

    damp_species_init(ch4);
    damp_species_init(ch4_cloud);
    damp_species_init(ch4_ice);


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

    // ATURAN_INIT_SATADJ_LAST=1 (2026-10-16, default 0 = as before): the four saturation adjustments
    // once more, in the loop's order and with no species filter after them, so that the state
    // init_PressureDynamic freezes p_baro from is the adjusted one. The adjustments above are each
    // followed by a Shapiro pass on vapour / cloud / ice, which leaves cells supersaturated, and the
    // first adjustment of the loop (iteration 2) then released latent heat: with ATURAN_HYDRO_SPLIT=1
    // ATURAN_BH_DIAG charged it -1.3e-5 (NH3) -4.5e-6 (H2S) -6.9e-6 (CH4) of the buoyancy pressure
    // force at 36 N, level 37 (NSA2), +0.08 K at level 20 -- measured on ATNEPT, ported 2026-10-16.
    static const bool init_satadj_last = [](){ const char* e = getenv("ATURAN_INIT_SATADJ_LAST"); return e && atoi(e) != 0; }();
    if(init_satadj_last){
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
    }

    // ATURAN_PDYN_GEO_BC=1 (2026-10-15, DEFAULT with ATURAN_GEO_BALANCE, 0 = the order before): the three boundary passes
    // also run BEFORE init_PressureDynamic, so that the balanced start (ATURAN_PDYN_GEO=2) is built
    // from the fields the first RK stage sees. Without it p_baro is frozen against the raw pole row,
    // bcTheta then resets t[0] = 4/3 t[1] - 1/3 t[2] and v = w = 0 there, the p_hydro the RK loop
    // rebuilds moves at j = 0, and the row next to the pole starts unbalanced: at level 38, 89 N
    // dp_hydro 2.969e-2 at init against 2.633e-2 at stage 0 of iteration 1 (dv/dt +3.4e-3, UDBG1).
    static const bool pdyn_geo_bc = [](){ const char* e = getenv("ATURAN_PDYN_GEO_BC"); return e ? atoi(e) != 0 : geo_balance(); }();
    if(pdyn_geo_bc){
        BC_Uran(*this).bcRadius();
        BC_Uran(*this).bcTheta();
        BC_Uran(*this).bcPhi();
    }

    init_PressureDynamic();

//    goto Printout;

    BC_Uran(*this).bcRadius();                                          // extrapolation in i-direction alomg grid boundaries
    BC_Uran(*this).bcTheta();                                           // extrapolation in j-direction alomg grid boundaries
    BC_Uran(*this).bcPhi();                                             // extrapolation in k-direction alomg grid boundaries

    tSeed();                                                            // ATURAN_T_SEED: off unless set
    restoreVar(1.0);

    // The closure SEEDS k* and dis* algebraically (Turbulence<Planet>::init_fields) before the
    // transport equations start carrying them. It belongs exactly here — after Forces,
    // Latent_Heat, the three boundary passes and the first restoreVar — because init_fields reads
    // the velocity field and the layer geometry, and none of that exists at the top of Run().
    // Putting the call there segfaults; this is ATSAT's slot (cSaturnModel.cpp:482, right after
    // its own pre-loop restoreVar) and ATNEPT's after it, for the same reason.
    if(turb_active) TurbulenceUran(*this).init();

//    goto Printout;



    // ATURAN_TB_DIAG=<n> (print-only, see TbUran.h): the level-mean temperature after every stage of
    // the loop that can write t, charged to that stage, in K; printed every n iterations.
    static const int tb_every = [](){ const char* e = getenv("ATURAN_TB_DIAG"); return e ? atoi(e) : 0; }();
    const bool tb_on = tb_every > 0 && im <= UranTb::MAXLEV;
    UranTb::st().on = tb_on;
    enum { TB_PRE = 0, TB_SATADJ, TB_PRECIP, TB_PHYS, TB_RK_RESET, TB_RK_REST, TB_BC, TB_CLAMP, TB_AFTER, TB_N };
    static const char* tb_name[TB_N] = { "pressure", "sat.adjust", "precip", "chem/latent", "RK reset", "RK rest",
        "BC", "clamp", "turb/conv" };
    static const char* tb_term[UranTb::NTERM] = { "RK press.work", "RK transp r", "RK transp the", "RK transp phi",
        "RK diffusion", "RK th.massflux", "RK radiation", "RK precip" };
    std::vector<std::vector<double>> tb_sum(TB_N, std::vector<double>(im, 0.0)), tb_tsum(UranTb::NTERM, std::vector<double>(im, 0.0));
    std::vector<double> tb_last(im, 0.0), tb_first(im, 0.0);
    double tb_norm = 0.0;
    auto tb_lev = [&](Array& a) -> std::vector<double> {
        std::vector<double> m_(im, 0.0);
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < im; i++){
            double s_ = 0.0;
            for(int j = 0; j < jm; j++){
                const double sj = std::sin(the.z[j]);
                for(int k = 0; k < km; k++) s_ += sj * a.x[i][j][k];
            }
            m_[i] = s_ * t_ref / tb_norm;
        }
        return m_; };
    auto tb_mark = [&](int stage){
        if(!tb_on) return;
        const std::vector<double> w_ = tb_lev(t);
        for(int i = 0; i < im; i++){ tb_sum[stage][i] += w_[i] - tb_last[i]; tb_last[i] = w_[i]; } };
    if(tb_on){
        for(int j = 0; j < jm; j++) tb_norm += std::sin(the.z[j]) * km;
        tb_last = tb_first = tb_lev(t);
    }

    // ATURAN_UB_DIAG=<n> (print-only, see UbUran.h): the level mean of the radial wind after every
    // stage of the loop that can write u, charged to that stage, in m/s; the Runge-Kutta step split
    // into the terms of rhs_u; and the radial pressure force with what the pressure solver and the
    // filter on p_dyn did to it. For the whole level and for the equatorial band (5 S .. 5 N by
    // default, ATURAN_UB_BAND=<degrees>). Printed every n iterations.
    static const int ub_every = [](){ const char* e = getenv("ATURAN_UB_DIAG"); return e ? atoi(e) : 0; }();
    static const int ub_band  = [](){ const char* e = getenv("ATURAN_UB_BAND"); return e ? atoi(e) : 5; }();
    const bool ub_on = ub_every > 0 && im <= UranUb::MAXLEV;
    UranUb::st().on = ub_on;
    // ATURAN_UB_LAT=<degrees north> (2026-10-13, default 0): the centre of the band, so that it can
    // be put where the residuum sits when that is not the equator.
    static const int ub_lat = [](){ const char* e = getenv("ATURAN_UB_LAT"); return e ? atoi(e) : 0; }();
    UranUb::st().j_lo = std::max(1,      (jm - 1) / 2 - ub_lat - ub_band);
    UranUb::st().j_hi = std::min(jm - 2, (jm - 1) / 2 - ub_lat + ub_band);
    // ATURAN_BH_DIAG=<n> (2026-10-16, print-only): the BUOYANCY part of rhs_v's -dp_hydro/dthe / r
    // at the ATURAN_VBUDGET level and meridian, recomputed after every stage of the loop that can
    // write t, the change charged to that stage, cumulative, printed every n iterations. Written for
    // ATURAN_HYDRO_SPLIT=1: the buoyancy integral moves in the first step after p_baro is frozen
    // (NRH_b2); ported 2026-10-16. The column integral, its factor and the level mean are those of
    // computeHydrostaticPressure(); nothing is written to the model.
    static const int bh_every = [](){ const char* e = getenv("ATURAN_BH_DIAG"); return e ? atoi(e) : 0; }();
    enum { BH_SOLVER = 0, BH_SAT_H2O, BH_SAT_H2S, BH_SAT_NH3, BH_SATADJ, BH_PRECIP, BH_PHYS, BH_RK, BH_TFILT, BH_BC, BH_CLAMP, BH_AFTER, BH_N };
    static const char* bh_name[BH_N] = { "solver/p filt", "sat.adj H2O", "sat.adj H2S", "sat.adj NH3", "sat.adj CH4", "precip", "chem/latent", "RK step",
        "t filter", "BC", "clamp", "turb/rad/conv" };
    const int bh_i = std::min(im - 2, std::max(1, vbudget_level_i()));
    const int bh_k = std::min(km - 2, std::max(1, vbudget_meridian_k()));
    std::vector<std::vector<double>> bh_sum(BH_N, std::vector<double>(jm, 0.0));
    std::vector<double> bh_last(jm, 0.0), bh_first(jm, 0.0);
    // the buoyancy column integral of meridian bh_k up to level bh_i, every row; minus its
    // meridional centred difference / r = the buoyancy's share of the T_DPH column of VBUDGET
    auto bh_term = [&]() -> std::vector<double> {
        const bool rho_level_ = [](){ const char* e = getenv("ATURAN_BUOY_RHO_LEVEL"); return e ? atoi(e) != 0 : geo_balance(); }();
        const bool nd_km_ = [](){ const char* e = getenv("ATURAN_HYDRO_ND_KM"); return e && atoi(e) != 0; }();
        const double bscale_ = [](){ const char* e = getenv("ATURAN_BUOY_SCALE"); return e ? atof(e) : 1.0; }();
        const double nd_ = 1.0e5 * (L_atm * (nd_km_ ? 1.0 : 1.0e3)) / (u_0 * u_0);
        const double nd_rel_ = g * (L_atm * 1.0e3) / (u_0 * u_0);
        std::vector<double> ref_(bh_i + 1, 0.0);
        #pragma omp parallel for schedule(static)
        for(int i = 0; i <= bh_i; i++){
            double s_ = 0.0, w_ = 0.0;
            for(int j = 0; j < jm; j++){
                const double wgt = sin(the.z[j]);
                for(int k = 0; k < km; k++){
                    if(!(t.x[i][j][k] > 0.0)) continue;
                    const double b = g * p_stat.x[i][j][k] / (r_mix * R_mix * t.x[i][j][k] * t_ref);
                    if(!std::isfinite(b)) continue;
                    s_ += wgt * b; w_ += wgt;
                }
            }
            ref_[i] = w_ > 0.0 ? s_ / w_ : 0.0;
        }
        std::vector<double> B_(jm, 0.0), out_(jm, 0.0);
        for(int j = 0; j < jm; j++){
            auto bu_ = [&](int i)->double{
                const double T = t.x[i][j][bh_k];
                if(!(T > 0.0)) return 0.0;
                const double b = g * p_stat.x[i][j][bh_k] / (r_mix * R_mix * T * t_ref) - ref_[i];
                const double f = rho_level_ ? -nd_rel_ * bscale_ * buoyancy * b / ref_[i]
                                            : -nd_ * bscale_ * buoyancy * b;
                return std::isfinite(f) ? f : 0.0; };
            double P_ = 0.0;
            for(int i = 1; i <= bh_i; i++) P_ += 0.5 * (bu_(i-1) + bu_(i)) * dr;
            B_[j] = P_;
        }
        const double rm_ = metricRadius(rad.z[bh_i]);
        for(int j = 1; j < jm - 1; j++) out_[j] = - (B_[j+1] - B_[j-1]) / (2.0 * dthe) / rm_;
        return out_; };
    auto bh_mark = [&](int stage){
        if(bh_every <= 0) return;
        const std::vector<double> w_ = bh_term();
        for(int j = 0; j < jm; j++){ bh_sum[stage][j] += w_[j] - bh_last[j]; bh_last[j] = w_[j]; } };
    if(bh_every > 0) bh_last = bh_first = bh_term();
    auto bh_print = [&](){
        static const int rows_[] = { 1, 2, 4, 6, 10, 20, 30, 40, 54, 70, 90, 110, 126, 140, 150, 160, 170, 174, 176, 178, 179 };
        printf("\n      ATURAN_BH_DIAG: buoyancy part of -dp_hydro/dthe / r at level i = %d, meridian k = %d, nondimensional,"
               " cumulative change over iterations 1..%d by stage\n        lat       start", bh_i, bh_k, iter_n);
        for(int s_ = 0; s_ < BH_N; s_++) printf(" %13s", bh_name[s_]);
        printf("         now\n");
        for(int r_ : rows_){
            if(r_ < 1 || r_ > jm - 2) continue;
            printf("      %5.0f %11.3e", 90.0 - (double)r_, bh_first[r_]);
            for(int s_ = 0; s_ < BH_N; s_++) printf(" %13.3e", bh_sum[s_][r_]);
            printf(" %11.3e\n", bh_last[r_]);
        }
    };

    enum { UB_RK_RESET = 0, UB_RK_REST, UB_FILTER, UB_BC, UB_N };
    static const char* ub_name[UB_N] = { "RK reset", "RK rest", "filter on u", "BC" };
    static const char* ub_term[UranUb::NTERM] = { "RK -dp/dr", "RK +(v2+w2)/r", "RK -advection", "RK diffusion",
        "RK Coriolis", "RK centrifugal", "RK buoyancy", "RK sponge" };
    static const char* ub_reg[UranUb::NREG] = { "whole level", "latitude band (ATURAN_UB_LAT, default the equator)" };
    typedef std::vector<std::vector<double>> ub_tab;                       // [region][level]
    std::vector<ub_tab> ub_sum(UB_N, ub_tab(UranUb::NREG, std::vector<double>(im, 0.0)));
    std::vector<ub_tab> ub_tsum(UranUb::NTERM, ub_tab(UranUb::NREG, std::vector<double>(im, 0.0)));
    ub_tab ub_last(UranUb::NREG, std::vector<double>(im, 0.0)), ub_first = ub_last;
    ub_tab ub_f_first = ub_last, ub_m_first = ub_last, ub_f_solver = ub_last, ub_f_filter = ub_last;
    double ub_norm[UranUb::NREG] = { 0.0, 0.0 };
    // the cells the Runge-Kutta step integrates in phi: k = 0..km-2 with the periodic seam
    const int ub_k0 = ATPhys::phi_first(ATPhys::phi_periodic<cUranusModel>());
    auto ub_in = [&](int r_, int j) -> bool {
        return r_ == UranUb::R_ALL || (j >= UranUb::st().j_lo && j <= UranUb::st().j_hi); };
    // level mean of a field [m/s]
    auto ub_lev = [&](Array& a) -> ub_tab {
        ub_tab m_(UranUb::NREG, std::vector<double>(im, 0.0));
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < im; i++){
            double s_[UranUb::NREG] = { 0.0, 0.0 };
            for(int j = 1; j < jm-1; j++){
                const double sj = std::sin(the.z[j]);
                double row_ = 0.0;
                for(int k = ub_k0; k < km-1; k++) row_ += a.x[i][j][k];
                for(int r_ = 0; r_ < UranUb::NREG; r_++) if(ub_in(r_, j)) s_[r_] += sj * row_;
            }
            for(int r_ = 0; r_ < UranUb::NREG; r_++) m_[r_][i] = s_[r_] * u_0 / ub_norm[r_];
        }
        return m_; };
    // level mean of the radial pressure force -dp_dyn/dr and of the metric term (v2 + w2)/r as
    // rhs_u forms them, in m/s per iteration
    auto ub_force = [&](bool metric) -> ub_tab {
        ub_tab m_(UranUb::NREG, std::vector<double>(im, 0.0));
        const double scale_p_ = 1.0e5 / (r_mix * u_0 * u_0);
        #pragma omp parallel for schedule(static)
        for(int i = 1; i < im-1; i++){
            double s_[UranUb::NREG] = { 0.0, 0.0 };
            const double inv_rm_ = 1.0 / metricRadius(rad.z[i]);
            for(int j = 1; j < jm-1; j++){
                const double sj = std::sin(the.z[j]);
                double row_ = 0.0;
                for(int k = ub_k0; k < km-1; k++)
                    row_ += metric ? (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k]) * inv_rm_
                                   : -(p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) / (2.0 * dr) * scale_p_;
                for(int r_ = 0; r_ < UranUb::NREG; r_++) if(ub_in(r_, j)) s_[r_] += sj * row_;
            }
            for(int r_ = 0; r_ < UranUb::NREG; r_++) m_[r_][i] = s_[r_] * dt * u_0 / ub_norm[r_];
        }
        return m_; };
    auto ub_mark = [&](int stage){
        if(!ub_on) return;
        const ub_tab w_ = ub_lev(u);
        for(int r_ = 0; r_ < UranUb::NREG; r_++) for(int i = 0; i < im; i++){
            ub_sum[stage][r_][i] += w_[r_][i] - ub_last[r_][i]; ub_last[r_][i] = w_[r_][i]; } };
    if(ub_on){
        for(int j = 1; j < jm-1; j++) for(int r_ = 0; r_ < UranUb::NREG; r_++)
            if(ub_in(r_, j)) ub_norm[r_] += std::sin(the.z[j]) * (km - 1 - ub_k0);
        ub_last = ub_first = ub_lev(u);
        ub_f_first = ub_force(false);
        ub_m_first = ub_force(true);
    }

    // ATURAN_CWB_DIAG=<n> (print-only, see CwbUran.h): the H2O, NH3, CH4 and H2S columns after
    // every stage of the loop that can write them, charged to that stage, in g/m2; printed every
    // n iterations. The stages are the ones ATURAN_TB_DIAG marks.
    static const int cwb_every = [](){ const char* e = getenv("ATURAN_CWB_DIAG"); return e ? atoi(e) : 0; }();
    const bool cwb_on = cwb_every > 0;
    UranCwb::st().on = cwb_on;
    static const char* cw_gas[UranCwb::NGAS] = { "H2O", "NH3", "CH4", "H2S" };
    static const char* cw_name[TB_N] = { "before the adjustment (pressure)", "SaturationAdjustment (all calls)",
        "Precipitation (in place)", "scalar sync, chemistry, forces, radiation", "RungeKutta: reset to the n-level copies",
        "RungeKutta: floor, limiters and remainder", "boundary conditions", "clampNegativeSpecies", "turbulence, convective adjustment, filters" };
    static const char* cw_term[UranCwb::NTERM] = { "RungeKutta: radial transport   -u dq/dr", "RungeKutta: meridional transport",
        "RungeKutta: zonal transport", "RungeKutta: diffusion", "RungeKutta: chemistry source", "RungeKutta: precipitation source",
        "RungeKutta: flux-form term  -s q div(u)", "(reference, not a term) q*div(u)" };
    Array* cw_f[UranCwb::NGAS][3]  = { { &h2o,  &h2o_cloud,  &h2o_ice  }, { &nh3,  &nh3_cloud,  &nh3_ice  }, { &ch4,  &ch4_cloud,  &ch4_ice  }, { &h2s,  &h2s_cloud,  &h2s_ice  } };
    Array* cw_fn[UranCwb::NGAS][3] = { { &h2on, &h2o_cloudn, &h2o_icen }, { &nh3n, &nh3_cloudn, &nh3_icen }, { &ch4n, &ch4_cloudn, &ch4_icen }, { &h2sn, &h2s_cloudn, &h2s_icen } };
    double cw_sum[UranCwb::NGAS][TB_N] = {{0.0}}, cw_tsum[UranCwb::NGAS][UranCwb::NTERM] = {{0.0}};
    double cw_norm = 0.0, cw_last[UranCwb::NGAS] = {0.0}, cw_first[UranCwb::NGAS] = {0.0}, cw_first0[UranCwb::NGAS] = {0.0};
    // i_hi = 0: every level with a thickness; i_hi = 1: the plane i = 0 alone.
    auto cwb_col = [&](Array* f[3], int i_hi) -> double {
        double s_ = 0.0;
        const int ie = i_hi > 0 ? i_hi : im;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:s_)
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                const double sj = std::sin(the.z[j]);
                for(int i = 0; i < ie; i++){
                    const double dz = layer_thickness_m(i);
                    if(dz > 0.0) s_ += sj * dz * (f[0]->x[i][j][k] + f[1]->x[i][j][k] + f[2]->x[i][j][k]);
                }
            }
        }
        return s_ * 1.0e3 / cw_norm; };
    auto cwb_mark = [&](int stage){
        if(!cwb_on) return;
        for(int g_ = 0; g_ < UranCwb::NGAS; g_++){
            const double w_ = cwb_col(cw_f[g_], 0);
            cw_sum[g_][stage] += w_ - cw_last[g_]; cw_last[g_] = w_;
        } };
    if(cwb_on){
        for(int j = 0; j < jm; j++) cw_norm += std::sin(the.z[j]) * km;
        for(int g_ = 0; g_ < UranCwb::NGAS; g_++){
            cw_last[g_] = cw_first[g_] = cwb_col(cw_f[g_], 0);
            cw_first0[g_] = cwb_col(cw_f[g_], 1);
        }
    }

    chromInfo();                                                        // ATURAN_CHROM: one line when the tracer is on
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

            ub_tab ub_f0, ub_f1;
            if(ub_on) ub_f0 = ub_force(false);
            if(press_solver_shared()) PressureSolver<cUranusModel>(*this).run();
            else                      PressureSolverUran(*this).run();
            if(ub_on) ub_f1 = ub_force(false);
            // ATURAN_DAMP_P_HORIZ=<s> (2026-10-14, ported from ATNEPT; default 0 with
            // ATURAN_GEO_BALANCE, 1.0 = the filter before): strength of the two horizontal passes of
            // the Shapiro filter on p_dyn; the vertical pass is unchanged. On ATNEPT the horizontal
            // passes halved the balancing meridional gradient in 224 iterations. damp_wiggles
            // applies its axes in the order k, j, i, so the split keeps the order.
            static const double damp_p_horiz = [](){
                const char* e = getenv("ATURAN_DAMP_P_HORIZ"); return e ? atof(e) : (geo_balance() ? 0.0 : 1.0); }();
            if(damp_p_horiz == 1.0)
                AtomUtils::damp_wiggles(p_dyn, nullptr, true, true, true);
            else {
                if(damp_p_horiz != 0.0) AtomUtils::damp_wiggles(p_dyn, nullptr, false, true, true, damp_p_horiz);
                AtomUtils::damp_wiggles(p_dyn, nullptr, true, false, false);
            }
            if(ub_on){
                const ub_tab ub_f2 = ub_force(false);
                for(int r_ = 0; r_ < UranUb::NREG; r_++) for(int i = 0; i < im; i++){
                    ub_f_solver[r_][i] += ub_f1[r_][i] - ub_f0[r_][i];
                    ub_f_filter[r_][i] += ub_f2[r_][i] - ub_f1[r_][i];
                }
            }
            tb_mark(TB_PRE); cwb_mark(TB_PRE); bh_mark(BH_SOLVER);
            tattrib("pressure+damp");

            SaturationAdjustmentUran(*this).run("H2O",
                coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
                t_0_h2o, t_00_h2o,
                ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
                C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
                h2o, h2o_cloud, h2o_ice);
            bh_mark(BH_SAT_H2O);

            SaturationAdjustmentUran(*this).run("H2S",
                coeff_h2s_A, coeff_h2s_B, coeff_h2s_A_i, coeff_h2s_B_i,
                t_0_h2s, t_00_h2s,
                ep_h2s, lv_h2s, ls_h2s, cp_h2s, r_h2s,
                C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, m_h2s,
                h2s, h2s_cloud, h2s_ice);
            bh_mark(BH_SAT_H2S);

            SaturationAdjustmentUran(*this).run("NH3",
                coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
                t_0_nh3, t_00_nh3,
                ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
                C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
                nh3, nh3_cloud, nh3_ice);
            bh_mark(BH_SAT_NH3);

            SaturationAdjustmentUran(*this).run("CH4",
                coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
                t_0_ch4, t_00_ch4,
                ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
                C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
                ch4, ch4_cloud, ch4_ice);

            tattrib("SaturationAdjust");

            // Must follow the saturation adjustments: it removes condensate in place.
            tb_mark(TB_SATADJ); cwb_mark(TB_SATADJ); bh_mark(BH_SATADJ);
            if(precip_enabled()) PrecipitationUran(*this).run();
            tb_mark(TB_PRECIP); cwb_mark(TB_PRECIP); bh_mark(BH_PRECIP);

            // ATURAN_RK_SCALAR_SYNC=<0|1|2> (2026-10-11; DEFAULT 2 on the user's word, 0 restores):
            // ATJUP's repair of 2026-10-08, as ported to ATSAT the same day. RungeKuttaUran
            // integrates every scalar from its n-level copy (h2on, ..., tn), and restoreVar
            // refreshed those at the END of the previous iteration -- before the saturation
            // adjustment above ran. So everything the adjustment writes directly (the condensate
            // it forms, the latent heat it releases) was overwritten by the step that follows.
            // 1 copies the twelve condensable fields into their n-level copies, so the step starts
            // from the adjusted state; 2 copies the temperature too, so the latent heat goes with
            // the condensate. The in-place precipitation above is kept by it as well.
            static const int rk_scalar_sync = [](){ const char* e = getenv("ATURAN_RK_SCALAR_SYNC"); return e ? atoi(e) : 2; }();
            if(rk_scalar_sync >= 1){
                #pragma omp parallel for collapse(2) schedule(static)
                for(int i = 0; i < im; i++){
                    for(int j = 0; j < jm; j++){
                        for(int k = 0; k < km; k++){
                            h2on.x[i][j][k]       = h2o.x[i][j][k];
                            h2o_cloudn.x[i][j][k] = h2o_cloud.x[i][j][k];
                            h2o_icen.x[i][j][k]   = h2o_ice.x[i][j][k];
                            h2sn.x[i][j][k]       = h2s.x[i][j][k];
                            h2s_cloudn.x[i][j][k] = h2s_cloud.x[i][j][k];
                            h2s_icen.x[i][j][k]   = h2s_ice.x[i][j][k];
                            nh3n.x[i][j][k]       = nh3.x[i][j][k];
                            nh3_cloudn.x[i][j][k] = nh3_cloud.x[i][j][k];
                            nh3_icen.x[i][j][k]   = nh3_ice.x[i][j][k];
                            ch4n.x[i][j][k]       = ch4.x[i][j][k];
                            ch4_cloudn.x[i][j][k] = ch4_cloud.x[i][j][k];
                            ch4_icen.x[i][j][k]   = ch4_ice.x[i][j][k];
                            if(rk_scalar_sync >= 2) tn.x[i][j][k] = t.x[i][j][k];
                        }
                    }
                }
            }
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
            tb_mark(TB_PHYS); cwb_mark(TB_PHYS); bh_mark(BH_PHYS);

        }  // if loop

        std::vector<double> tb_n;
        if(tb_on){
            tb_n = tb_lev(tn);                                          // what the step starts from
            for(int t_ = 0; t_ < UranTb::MAXTHR; t_++) for(int q_ = 0; q_ < UranTb::NTERM; q_++)
                for(int i = 0; i < im; i++) UranTb::st().acc[t_][q_][i] = 0.0;
        }

        double cw_n[UranCwb::NGAS] = {0.0};
        if(cwb_on){
            for(int g_ = 0; g_ < UranCwb::NGAS; g_++) cw_n[g_] = cwb_col(cw_fn[g_], 0);    // what the step starts from
            for(int t_ = 0; t_ < UranCwb::MAXTHR; t_++) for(int g_ = 0; g_ < UranCwb::NGAS; g_++)
                for(int q_ = 0; q_ < 8; q_++) UranCwb::st().acc[t_][g_][q_] = 0.0;
        }

        ub_tab ub_n;
        if(ub_on){
            ub_n = ub_lev(un);                                          // what the step starts from
            for(int t_ = 0; t_ < UranUb::MAXTHR; t_++) for(int r_ = 0; r_ < UranUb::NREG; r_++)
                for(int q_ = 0; q_ < UranUb::NTERM; q_++) for(int i = 0; i < im; i++) UranUb::st().acc[t_][r_][q_][i] = 0.0;
        }

        chromUV();                                        // ATURAN_CHROM: the UV field the source reads
        RungeKuttaUran();
        zdevDiag();                                       // ATURAN_ZDEV: print-only
        bh_mark(BH_RK);

        if(ub_on){
            const ub_tab w_ = ub_lev(u);
            for(int r_ = 0; r_ < UranUb::NREG; r_++) for(int i = 0; i < im; i++){
                double step_terms = 0.0;
                for(int q_ = 0; q_ < UranUb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < UranUb::MAXTHR; t_++) a_ += UranUb::st().acc[t_][r_][q_][i];
                    a_ *= dt / 6.0 * u_0 / ub_norm[r_];
                    ub_tsum[q_][r_][i] += a_; step_terms += a_;
                }
                ub_sum[UB_RK_RESET][r_][i] += ub_n[r_][i] - ub_last[r_][i];
                ub_sum[UB_RK_REST][r_][i]  += w_[r_][i] - ub_n[r_][i] - step_terms;
                ub_last[r_][i] = w_[r_][i];
            }
        }


        if(cwb_on){
            for(int g_ = 0; g_ < UranCwb::NGAS; g_++){
                const double w_ = cwb_col(cw_f[g_], 0);
                double step_terms = 0.0;
                for(int q_ = 0; q_ < UranCwb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < UranCwb::MAXTHR; t_++) a_ += UranCwb::st().acc[t_][g_][q_];
                    a_ *= dt / 6.0 * 1.0e3 / cw_norm;
                    cw_tsum[g_][q_] += a_;
                    if(q_ != UranCwb::T_QDIV) step_terms += a_;
                }
                cw_sum[g_][TB_RK_RESET] += cw_n[g_] - cw_last[g_];
                cw_sum[g_][TB_RK_REST]  += w_ - cw_n[g_] - step_terms;
                cw_last[g_] = w_;
            }
        }

        if(tb_on){
            const std::vector<double> w_ = tb_lev(t);
            for(int i = 0; i < im; i++){
                double step_terms = 0.0;
                for(int q_ = 0; q_ < UranTb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < UranTb::MAXTHR; t_++) a_ += UranTb::st().acc[t_][q_][i];
                    a_ *= dt / 6.0 * t_ref / tb_norm;
                    tb_tsum[q_][i] += a_; step_terms += a_;
                }
                tb_sum[TB_RK_RESET][i] += tb_n[i] - tb_last[i];
                tb_sum[TB_RK_REST][i]  += w_[i] - tb_n[i] - step_terms;
                tb_last[i] = w_[i];
            }
        }
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
        // DEFAULT 0.0 SINCE 2026-10-11 on the user's word (1.0 restores the filter of every run
        // before that). MEASURED with ATURAN_TB_DIAG, 224 iterations, level-mean change in K at levels
        // 1 / 20 / 35 / 39 (satchk/giants/URANtbB against URANtbE):
        //     vertical pass 1.0, thermal mass flux 1     -65.8 / -7.8 / +37.5 / +69.5
        //     vertical pass 0.0, thermal mass flux 0     -2.9 / -0.5 / +1.1 / +2.3
        // mean OLR 21.1 -> 1.51 W/m2 against 0.689 in, photosphere 137.0 -> 69.1 K; exit 0, no non-finite
        // value. With the pass off and the thermal mass flux still on, that term digs a cold hole
        // (-17 K at level 28), which is why ATURAN_THERMAL_MASSFLUX went to 0 in the same change.
        // What is left is the diffusion term at the two radial ends (the 2-point wall value does
        // not continue a straight lapse) and a rougher profile: the README's sweep measured the
        // mean second difference 3.1 times larger at 0.0 than at 1.0, with nothing oscillating.
        static const double damp_t_vert = [](){
            const char* e = getenv("ATURAN_DAMP_T_VERT"); return e ? atof(e) : 0.0; }();
        // ATURAN_DAMP_T_HORIZ=<strength> (2026-10-15, default 1.0 = the pass of every run before):
        // the horizontal pass on t. The latitude pass clamps at the poles, j = 0 takes
        // 0.25 (T[1] - T[0]) every iteration, so the cold pole row under the lid is pulled toward
        // its neighbour -- the suspect for the polar warming of the GEO_BALANCE stack (UST1).
        // ATURAN_DAMP_T_J=0 keeps the longitude pass and drops the latitude one. DEFAULT 0 WITH
        // ATURAN_GEO_BALANCE since 2026-10-15 (1 otherwise). MEASURED (UST1 / UDTJ0_4, 4 it): the
        // pole row at the lid warmed 0.045 K per iteration with the pass, 3e-5 without, and dv/dt
        // at 89 N stopped growing (3.4 -> 6.7e-3 in 4 it with it, flat without).
        static const double damp_t_horiz = [](){
            const char* e = getenv("ATURAN_DAMP_T_HORIZ"); return e ? atof(e) : 1.0; }();
        static const bool damp_t_j = [](){
            const char* e = getenv("ATURAN_DAMP_T_J"); return e ? atoi(e) != 0 : !geo_balance(); }();
        if(damp_t_horiz != 0.0)
            AtomUtils::damp_wiggles(t, nullptr, false, damp_t_j, true, damp_t_horiz);
        if(damp_t_vert != 0.0)
            AtomUtils::damp_wiggles(t, nullptr, true, false, false, damp_t_vert);
        bh_mark(BH_TFILT);
        tattrib("damp_wiggles(t)");
        // ATURAN_VEL_SHAPIRO=<strength> (2026-10-14, ported from ATNEPT; default 0 with
        // ATURAN_GEO_BALANCE, 1.0 = the pass of every run before): the per-iteration Shapiro pass on
        // u, v, w; ATURAN_VEL_SHAPIRO_ORDER=4 takes ATJUP's 4th-order pass. On ATNEPT the 1-2-1 pass
        // took two thirds of the high-latitude jets in 224 iterations and the 4th-order one still
        // moved the corners of the piecewise-linear jets. ATJUP runs no pass in its loop by default.
        static const double vel_shapiro = [](){
            const char* e = getenv("ATURAN_VEL_SHAPIRO"); return e ? atof(e) : (geo_balance() ? 0.0 : 1.0); }();
        static const int vel_order = [](){
            const char* e = getenv("ATURAN_VEL_SHAPIRO_ORDER"); return e ? atoi(e) : 2; }();
        auto vel_filter = [&](Array& f){
            if(vel_order == 4) AtomUtils::damp_wiggles_ho(f, nullptr, true, true, true, vel_shapiro);
            else               AtomUtils::damp_wiggles   (f, nullptr, true, true, true, vel_shapiro); };
        if(vel_shapiro != 0.0) vel_filter(u);
        ub_mark(UB_FILTER);
        if(vel_shapiro != 0.0){
            vel_filter(v);
            vel_filter(w);
        }

        BC_Uran(*this).bcRadius();                                      // extrapolation in i-direction along grid boundaries
        BC_Uran(*this).bcTheta();                                       // extrapolation in j-direction along grid boundaries
        BC_Uran(*this).bcPhi();                                         // extrapolation in k-direction along grid boundaries
        ub_mark(UB_BC);
        tb_mark(TB_BC); cwb_mark(TB_BC); bh_mark(BH_BC);
        clampNegativeSpecies();                                         // ATURAN_SPECIES_CLAMP, default on (BC_Uran.h)
        tb_mark(TB_CLAMP); cwb_mark(TB_CLAMP); bh_mark(BH_CLAMP);
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
        if(mca_stride() > 0 && iter_n % mca_stride() == 0){
            typedef MoistConvectiveAdjustmentUran::Gas    Gas_;
            typedef MoistConvectiveAdjustmentUran::Tracer Tr_;
            // H2S has no ice coefficient set here; its liquid set stands for both phases, as in the model's own adjustment.
            std::vector<Gas_> gases = {
                { "H2O", &h2o, &h2o_cloud, &h2o_ice, ep_h2o, lv_h2o, ls_h2o, t_0_h2o,
                  C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, C_h2o_ice, L0_h2o_ice, del_alf_h2o_ice, del_bet_h2o_ice },
                { "NH3", &nh3, &nh3_cloud, &nh3_ice, ep_nh3, lv_nh3, ls_nh3, t_0_nh3,
                  C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, C_nh3_ice, L0_nh3_ice, del_alf_nh3_ice, del_bet_nh3_ice },
                { "H2S", &h2s, &h2s_cloud, &h2s_ice, ep_h2s, lv_h2s, ls_h2s, t_0_h2s,
                  C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, C_h2s, L0_h2s, del_alf_h2s, del_bet_h2s },
                { "CH4", &ch4, &ch4_cloud, &ch4_ice, ep_ch4, lv_ch4, ls_ch4, t_0_ch4,
                  C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, C_ch4_ice, L0_ch4_ice, del_alf_ch4_ice, del_bet_ch4_ice } };
            std::vector<Tr_> tracers = { { "NH4SH", &nh4sh, 0.0, true, 0.0 } };
            if(chrom_on()) tracers.push_back({ "chrom", &chrom, 0.0, false, 0.0 });   // ATURAN_CHROM: mixed as ATJUP mixes it
            MoistConvectiveAdjustmentUran(*this, gases, tracers).run();
        }
        tb_mark(TB_AFTER); cwb_mark(TB_AFTER); bh_mark(BH_AFTER);
        if(bh_every > 0 && iter_n % bh_every == 0) bh_print();
        if(tb_on && iter_n % tb_every == 0){
            static const int lv_[] = { 0, 1, 2, 3, 5, 8, 12, 16, 20, 24, 28, 32, 35, 37, 38, 39, 40 };
            const int nl_ = (int)(sizeof(lv_) / sizeof(lv_[0]));
            printf("\n      ATURAN: [TB] level-mean temperature, cumulative over iterations %d..%d, K\n      ATURAN: [TB] %-15s", 1, iter_n, "level");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %8d", lv_[l]);
            printf("\n      ATURAN: [TB] %-15s", "T now");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %8.2f", tb_last[lv_[l]]);
            printf("\n      ATURAN: [TB] %-15s", "change");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_last[lv_[l]] - tb_first[lv_[l]]);
            for(int q_ = 0; q_ < TB_N; q_++){
                printf("\n      ATURAN: [TB] %-15s", tb_name[q_]);
                for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_sum[q_][lv_[l]]);
                if(q_ == TB_RK_RESET) for(int r_ = 0; r_ < UranTb::NTERM; r_++){
                    printf("\n      ATURAN: [TB] %-15s", tb_term[r_]);
                    for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_tsum[r_][lv_[l]]);
                }
            }
            printf("\n");
        }

        tattrib("Rad/Turb/ConvAdj");

        if(ub_on && iter_n % ub_every == 0){
            static const int lv_[] = { 0, 1, 2, 3, 5, 8, 12, 16, 20, 24, 28, 32, 35, 37, 38, 39, 40 };
            const int nl_ = (int)(sizeof(lv_) / sizeof(lv_[0]));
            const ub_tab f_now = ub_force(false), m_now = ub_force(true);
            for(int r_ = 0; r_ < UranUb::NREG; r_++){
                auto row_ = [&](const char* name, const std::vector<double>& x_, const char* fmt){
                    printf("\n      ATURAN: [UB] %-16s", name);
                    for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(fmt, x_[lv_[l]]); };
                printf("\n      ATURAN: [UB] radial wind, %s (rows %d..%d), cumulative over iterations %d..%d, m/s\n      ATURAN: [UB] %-16s",
                       ub_reg[r_], r_ == UranUb::R_ALL ? 1 : UranUb::st().j_lo, r_ == UranUb::R_ALL ? jm-2 : UranUb::st().j_hi, 1, iter_n, "level");
                for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %9d", lv_[l]);
                row_("u now", ub_last[r_], " %+9.5f");
                std::vector<double> d_(im, 0.0);
                for(int i = 0; i < im; i++) d_[i] = ub_last[r_][i] - ub_first[r_][i];
                row_("change", d_, " %+9.5f");
                for(int q_ = 0; q_ < UB_N; q_++){
                    row_(ub_name[q_], ub_sum[q_][r_], " %+9.5f");
                    if(q_ == UB_RK_RESET) for(int t_ = 0; t_ < UranUb::NTERM; t_++) row_(ub_term[t_], ub_tsum[t_][r_], " %+9.5f");
                }
                printf("\n      ATURAN: [UB] the radial force in the same region, m/s per iteration (x 1e3)");
                std::vector<double> x_(im, 0.0);
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * f_now[r_][i];                 row_("-dp/dr now", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * m_now[r_][i];                 row_("(v2+w2)/r now", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * (f_now[r_][i] + m_now[r_][i]); row_("sum of the two", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_m_first[r_][i];            row_("(v2+w2)/r start", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_f_first[r_][i];            row_("-dp/dr at start", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_f_solver[r_][i];           row_("  solver", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_f_filter[r_][i];           row_("  filter on p", x_, " %+9.4f");
                printf("\n");
            }
        }

        if(cwb_on && iter_n % cwb_every == 0){
            for(int g_ = 0; g_ < UranCwb::NGAS; g_++){
                double sum_ = 0.0;
                printf("\n      ATURAN: [CWB] %s column (vapour+cloud+ice), cumulative over iterations 1..%d, g/m2:  start %.9e  now %.9e  change %+.6e (%+.4f %%)\n",
                       cw_gas[g_], iter_n, cw_first[g_], cw_last[g_], cw_last[g_] - cw_first[g_], 100.0 * (cw_last[g_] / cw_first[g_] - 1.0));
                const double p0_ = cwb_col(cw_f[g_], 1) - cw_first0[g_], in_ = (cw_last[g_] - cw_first[g_]) - p0_;
                printf("      ATURAN: [CWB] %s   of the change: plane i=0 (set by bcRadius) %+.6e, integrated cells i>=1 %+.6e (%+.4f %% of their start)\n",
                       cw_gas[g_], p0_, in_, 100.0 * in_ / (cw_first[g_] - cw_first0[g_]));
                for(int q_ = 0; q_ < TB_N; q_++){
                    printf("      ATURAN: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], cw_name[q_], cw_sum[g_][q_]); sum_ += cw_sum[g_][q_];
                    if(q_ == TB_RK_RESET)
                        for(int r_ = 0; r_ < UranCwb::NTERM; r_++){
                            printf("      ATURAN: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], cw_term[r_], cw_tsum[g_][r_]);
                            if(r_ != UranCwb::T_QDIV) sum_ += cw_tsum[g_][r_];
                        }
                }
                printf("      ATURAN: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], "unattributed", cw_last[g_] - cw_first[g_] - sum_);
            }
        }

        panorama_cnt++;

        // ATURAN_VTK_STRIDE=<n> (2026-10-12, ATJUP's and ATSAT's knob): write the VTK files every n
        // iterations instead of at every checkpoint, so the printed extrema can stay dense while
        // the files do not. Default 0 = follow `checkpoint`, as before. This model writes no
        // periodic restart files, so there is no restart stride to go with it.
        static const int vtk_stride = [](){
            const char* e = getenv("ATURAN_VTK_STRIDE"); return e ? atoi(e) : 0; }();
        if(iter_n % checkpoint == 0){
            printMinMax();
            printMeridionalBudget();         // no-op unless ATURAN_VBUDGET is set
            printTemperatureBudget();        // no-op unless ATURAN_TBUDGET is set
            printTemperatureAttribution();   // no-op unless ATURAN_TATTRIB is set
        }
        if(vtk_stride > 0 ? (iter_n % vtk_stride == 0) : (iter_n % checkpoint == 0)) writeData();
        // ATURAN_CLOUDVIEW_STRIDE=<n> > 0 (cUranusModel.h): the cloudview on its own cadence; unset rides with the panorama in writeData.
        if(cloudview_stride() > 0 && iter_n % cloudview_stride() == 0) paraview_cloudview(iter_n);

        if(panorama_cnt == panorama_print) panorama_cnt = 0;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for one time step\n", elapsed.count() * 1e-9);

    }  // end for iter_n

    clampNegativeReport();

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
    chrom.initArray(im, jm, km, 0.0);                 // chromophore (ATURAN_CHROM)
    chrom_uv.initArray(im, jm, km, 0.0);              // its UV fraction

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
    chromn.initArray(im, jm, km, 0.0);                // chromophore new

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
    acc_chrom.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
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
    P_h2s_rain.initArray(im, jm, km, 0.0);
    P_h2s_snow.initArray(im, jm, km, 0.0);
    P_h2s_graupel.initArray(im, jm, km, 0.0);
    S_precip_h2s.initArray(im, jm, km, 0.0);
    S_precip_h2s_cloud.initArray(im, jm, km, 0.0);
    S_precip_h2s_ice.initArray(im, jm, km, 0.0);
    precip_srf_h2s.initArray_2D(jm, km, 0.0);
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
    rhs_chrom.initArray(im, jm, km, 0.0);                // auxilliar field RHS chromophore
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
                chromn.x[i][j][k] = coeff * chrom.x[i][j][k];
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



/*
 * SHARED DIAGNOSTICS — the min/max report machinery and the steady-state query. MUST BE
 * BYTE-IDENTICAL IN EVERY MODEL THAT USES IT; `make check-shared` verifies that against
 * planet/SHARED.md5.
 *
 * This is the first shared header that is not physics. It exists because of a measurement: of
 * the ATSAT code that was not yet shared, 58 % still matched ATJUP line-for-line after planet
 * names were normalised away, and the largest blocks of that were not the physics — they were
 * the reporting scaffolding. PrintMsg_Sat.cpp and PrintMsg_Jup.cpp overlapped 83.8 %, the
 * highest of any pair in either model. That is the code a third and fourth planet copy again,
 * so it is the code whose duplication compounds fastest.
 *
 * ===== WHAT IS HERE AND WHAT IS NOT =====
 *
 * Here: searchMinMax_3D, searchMinMax_2D and steadyQuery — pure machinery, no planet in them.
 *
 * NOT here, deliberately: printMinMax() itself. That routine is a LIST of which fields to print,
 * in what unit, with what coefficient, and the two models genuinely disagree about it — ATJUP
 * prints the species with coefficient 1.0 as kg/m3, ATSAT multiplies by 1e3*r_mix and labels
 * g/m3. ATJUP's note records that its own factor was wrong by 1.2844x and was corrected; ATSAT
 * has not had that correction, and its temperature line scales t by 273.15 where t_ref (134) is
 * meant. Sharing printMinMax would have to pick one reading of the species units for both models
 * silently. It is a real question and it gets its own commit; see the note left in printMinMax.
 *
 * Also not here: print_welcome_msg, initMsg, print_final_msg — a few lines each, almost entirely
 * planet-specific text.
 *
 * ===== THE THREE DIVERGENCES THIS ABSORBS =====
 *
 * 1. THE METRIC RADIUS. The continuity residual used
 *        ATSAT   metricRadius(rad.z[i])      ATJUP   rad.z[i]
 *    and metricRadius() is the established hook for exactly this — FluxLimiter and PressureSolver
 *    ask the model the same question. ATJUP's is the identity by design (it shifts rad.z itself at
 *    initialisation), so the shared call is bit-identical in both. See FluxLimiter.h, which
 *    explains why the two models answer this question from opposite ends.
 *
 * 2. THE POLAR COSINE. ATJUP flips the sign of cos(theta) in the southern hemisphere
 *    (`if(costhe_abs() && j > 90) costhe = -costhe`); ATSAT does not. That is already a named
 *    static on cJupiterModel, so ATSAT gains the same accessor returning false and the branch is
 *    dead there — the divergence becomes a declared answer rather than a missing line.
 *
 * 3. THE COLUMN WIDTHS. ATJUP prints the unit column at setw(12) and separates the max and min
 *    halves with three spaces; ATSAT uses setw(6) and ten spaces. ATJUP widened its column
 *    because its unit strings are longer (" kg/(m3s)"). This is cosmetic and the two are
 *    preserved exactly, through minmax_unit_width() and minmax_separator(), so that this commit
 *    changes NO output in either model. Unifying them would be an improvement and is a separate,
 *    visible decision — not something to slip in under a refactor.
 *
 * ===== WHAT CHANGED IN THE CODE ITSELF =====
 *
 * steadyQuery was 270 lines of which about 200 were the same twelve-line block written twelve
 * times: scan a field against its n-copy, keep the running maximum and its cell, then a
 * thirteen-case switch to print them back in order. The field list is IDENTICAL in both models —
 * same twelve fields, same thirteen labels, same order — so it is now written once as a table and
 * walked twice. This is the same treatment BC_Sat.cpp had in 34205f2 ("say the boundary
 * conditions once each, as field lists, instead of 110 times longhand") and for the same reason:
 * a list that is data can be checked against the model's field set, and a list that is code
 * cannot.
 *
 * The scan order is unchanged and so is the >= comparison, which matters: >= means the LAST cell
 * holding the extreme value wins, so the reported location depends on the traversal order. The
 * table is walked inside the same i,j,k loop nest the twelve copies shared, in the same order.
 */

#pragma once

#include "ATPhys.h"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

template<class Planet>
class Reporting{

    Planet &m;

    // Latitude/longitude in the report's own convention: j runs from the north pole, k from the
    // prime meridian eastward, and both halves of every line are labelled with a hemisphere.
    struct HemisphereCoords{
        double lat, lon;
        std::string east_or_west, north_or_south;
    };

    static HemisphereCoords convert_coords(double lon, double lat){
        HemisphereCoords ret;
        if(lat > 90){
            ret.lat = lat - 90;
            ret.north_or_south = "°S";
        }else{
            ret.lat = 90 - lat;
            ret.north_or_south = "°N";
        }
        if(lon > 180){
            ret.lon = 360 - lon;
            ret.east_or_west = "°W";
        }else{
            ret.lon = lon;
            ret.east_or_west = "°E";
        }
        return ret;
    }

    // The one line both searchMinMax routines print. Pulled out because it was written twice
    // verbatim and the two copies had already drifted in the 2D version's level handling.
    void print_pair(const std::string &name_maxValue, const std::string &name_minValue,
        const std::string &name_unitValue, double maxValue, double minValue,
        int imax_level, int jmax, int kmax, int imin_level, int jmin, int kmin) const
    {
        const int uw = Planet::minmax_unit_width();
        const char *sep = Planet::minmax_separator();
        const std::string level = "km";

        HemisphereCoords cmax = convert_coords(kmax, jmax);
        HemisphereCoords cmin = convert_coords(kmin, jmin);

        std::cout << std::setiosflags(std::ios::left) << std::setw(26) << std::setfill('.')
            << name_maxValue << " = "
            << std::resetiosflags(std::ios::left) << std::setw(12) << std::fixed
            << std::setfill(' ') << maxValue << std::setw(uw)
            << name_unitValue << std::setw(5) << (int)cmax.lat << std::setw(3)
            << cmax.north_or_south << std::setw(4) << (int)cmax.lon
            << std::setw(3) << cmax.east_or_west << std::setw(6) << imax_level
            << std::setw(2) << level << sep
            << std::setiosflags(std::ios::left) << std::setw(26) << std::setfill('.')
            << name_minValue << " = "
            << std::resetiosflags(std::ios::left) << std::setw(12) << std::fixed
            << std::setfill(' ') << minValue << std::setw(uw)
            << name_unitValue << std::setw(5) << (int)cmin.lat << std::setw(3)
            << cmin.north_or_south << std::setw(4) << (int)cmin.lon
            << std::setw(3) << cmin.east_or_west << std::setw(6) << imin_level
            << std::setw(2) << level << std::endl;
    }

public:

    explicit Reporting(Planet &m_) : m(m_) {}


    /*
     * The min/max report's COMMON LIST — the 39 fields every model carries, in one place.
     *
     * This could not be shared until the units were settled, and they disagreed on 35 of these 39
     * rows: a species density counted twice, W/m2 against W/m3, four force prefixes, three
     * readings of thermalmassflux, and a temperature scaled by 273.15 instead of t_ref. Each was
     * settled on its own terms in its own commit; what is left is a list the three models agree
     * on, so it lives here.
     *
     * NOT here: each planet's own fields. ATSAT and ATJUP print radiation, turbulence and
     * precipitation arrays ATNEPT has none of; ATNEPT prints rho_mix and an H2S condensate pair
     * the others lack. Each model calls this and then prints its own, which is why the report now
     * has ONE ORDER everywhere — shared sections first, the planet's own after. That reordering is
     * the only visible change to any log.
     */
    void print_minmax_common(){
        std::cout << std::endl << std::endl << " Courant time step   dt = " << m.dt
                  << std::endl << std::endl;

        std::cout << std::endl << std::endl << " Temperatures " << std::endl;
        searchMinMax_3D(" max 3D temperature ", " min 3D temperature ", " degC",
            m.t, m.t_ref, [](double i)->double{return i - 273.15;}, true);
        searchMinMax_3D(" max 3D thermalflux ", " min 3D thermalflux ", " W/m3", m.thermalmassflux, 1.0);
        std::cout << std::endl;

        std::cout << std::endl << " Velocities " << std::endl;
        searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", " m/s", m.u, m.u_0);
        searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", " m/s", m.v, m.u_0);
        searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", " m/s", m.w, m.u_0);
        std::cout << std::endl;

        std::cout << std::endl << " Pressures " << std::endl;
        searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", " bar", m.p_dyn, m.p_dyn_to_bar());
        searchMinMax_3D(" max 3D pressure static ", " min 3D pressure static ", " bar", m.p_stat, 1.0);
        std::cout << std::endl;

        // ---- DISPLAY UNITS, settled across all four models 2026-08-04 ----
        //
        // One unit per quantity, used by this report AND by every model's ParaView field lists,
        // so an array name means the same thing wherever it is read. Before this the four writers
        // disagreed with each other and with this table: ATJUP and ATNEPT passed 1.0 (kg/m3) for
        // the species while ATSAT and ATURAN passed 1e3 (g/m3), and NH4SH was 1e3 in ATSAT and
        // 1e6 in the other three.
        //
        //     species, clouds, ices and their fluxes   1e3   g/m3, g/(m3s), g/m4
        //     the NH4SH family                         1e9   ug/m3, ug/(m3s), ug/m4
        //     thermalmassflux                          1.0   W/m3   (a heating rate, not a species)
        //
        // THE UNITS ARE CHOSEN SO EVERY MODEL IS LEGIBLE UNDER THE WRITER'S precision(4) FIXED
        // FORMAT, which is the real constraint. Measured maxima at the time of the settlement:
        //
        //     h2o   [g/m3]   ATNEPT 18.76   ATURAN 80.00   ATSAT 90.00   ATJUP 90.86
        //     nh4sh [ug/m3]  ATNEPT  0.034  ATURAN  1.512  ATJUP 6730.5  ATSAT 4764984.6
        //
        // NH4SH spans EIGHT ORDERS OF MAGNITUDE across the four, which is why it needs a unit of
        // its own: at the species' 1e3 it would print 0.0000 on Neptune and Uranus both, and at
        // the 1e6 this table used before, Neptune's still did. ug/m3 is the first unit that keeps
        // all four non-zero. If a fifth planet ever breaks that, the answer is the writer's
        // precision, not another unit.
        std::cout << std::endl << " Water " << std::endl;
        searchMinMax_3D(" max 3D h2o ",  " min 3D h2o ", " g/m3", m.h2o, 1e3);
        searchMinMax_3D(" max 3D h2o_cloud ", " min 3D h2o_cloud ", " g/m3", m.h2o_cloud, 1e3);
        searchMinMax_3D(" max 3D h2o_ice ", " min 3D h2o_ice ", " g/m3", m.h2o_ice, 1e3);
        std::cout << std::endl;

        std::cout << std::endl << " Methane " << std::endl;
        searchMinMax_3D(" max 3D ch4 ",  " min 3D ch4 ", " g/m3", m.ch4, 1e3);
        searchMinMax_3D(" max 3D ch4_cloud ", " min 3D ch4_cloud ", " g/m3", m.ch4_cloud, 1e3);
        searchMinMax_3D(" max 3D ch4_ice ", " min 3D ch4_ice ", " g/m3", m.ch4_ice, 1e3);
        std::cout << std::endl;

        std::cout << std::endl << " Hydrogen Sulfide " << std::endl;
        searchMinMax_3D(" max 3D h2s ",  " min 3D h2s ", " g/m3", m.h2s, 1e3);
        searchMinMax_3D(" max 3D w_h2s ", " min 3D w_h2s ", " g/(m3s)", m.w_h2s, 1e3);
        searchMinMax_3D(" max 3D j_h2s ", " min 3D j_h2s ", " g/m4", m.j_h2s, 1e3);
        searchMinMax_3D(" max 3D jT_h2s ", " min 3D jT_h2s ", " g/m4", m.jT_h2s, 1e3);
        searchMinMax_3D(" max 3D massflux_h2s ", " min 3D massflux_h2s ", " g/(m3s)", m.massflux_h2s, 1e3);
        searchMinMax_3D(" max 3D diff_h2s ", " min 3D diff_h2s ", " g/(m3s)", m.difflux_h2s, 1e3);
        std::cout << std::endl;

        std::cout << std::endl << " Ammonia " << std::endl;
        searchMinMax_3D(" max 3D nh3 ",  " min 3D nh3 ", " g/m3", m.nh3, 1e3);
        searchMinMax_3D(" max 3D nh3_cloud ", " min 3D nh3_cloud ", " g/m3", m.nh3_cloud, 1e3);
        searchMinMax_3D(" max 3D nh3_ice ", " min 3D nh3_ice ", " g/m3", m.nh3_ice, 1e3);
        searchMinMax_3D(" max 3D w_nh3 ", " min 3D w_nh3 ", " g/(m3s)", m.w_nh3, 1e3);
        searchMinMax_3D(" max 3D j_nh3 ", " min 3D j_nh3 ", " g/m4", m.j_nh3, 1e3);
        searchMinMax_3D(" max 3D jT_nh3 ", " min 3D jT_nh3 ", " g/m4", m.jT_nh3, 1e3);
        searchMinMax_3D(" max 3D massflux_nh3 ", " min 3D massflux_nh3 ", " g/(m3s)", m.massflux_nh3, 1e3);
        searchMinMax_3D(" max 3D diff_nh3 ", " min 3D diff_nh3 ", " g/(m3s)", m.difflux_nh3, 1e3);
        std::cout << std::endl;

        std::cout << std::endl << " Ammonia Hydrosulfide " << std::endl;
        searchMinMax_3D(" max 3D nh4sh ",  " min 3D nh4sh ", " ug/m3", m.nh4sh, 1e9);
        searchMinMax_3D(" max 3D w_nh4sh ", " min 3D w_nh4sh ", " ug/(m3s)", m.w_nh4sh, 1e9);
        searchMinMax_3D(" max 3D j_nh4sh ", " min 3D j_nh4sh ", " ug/m4", m.j_nh4sh, 1e9);
        searchMinMax_3D(" max 3D jT_nh4sh ", " min 3D jT_nh4sh ", " ug/m4", m.jT_nh4sh, 1e9);
        searchMinMax_3D(" max 3D massflux_nh4sh ", " min 3D massflux_nh4sh ", " ug/(m3s)", m.massflux_nh4sh, 1e9);
        searchMinMax_3D(" max 3D diff_nh4sh ", " min 3D diff_nh4sh ", " ug/(m3s)", m.difflux_nh4sh, 1e9);
        std::cout << std::endl;

        std::cout << std::endl << " Forces " << std::endl;
        searchMinMax_3D(" max 3D Coriolis force ", " min 3D Coriolis force ", " mN/m3", m.CoriolisForce, 1e3);
        searchMinMax_3D(" max 3D centrifugal force ", " min 3D centrifugal force ", " mN/m3", m.CentrifugalForce, 1e3);
        searchMinMax_3D(" max 3D buoyancy force ", " min 3D buoyancy force ", " N/m3", m.BuoyancyForce, 1.0);
        searchMinMax_3D(" max 3D presgrad force ", " min 3D presgrad force ", " N/m3", m.PresGradForce, 1.0);
        std::cout << std::endl;

        std::cout << std::endl << " Energies " << std::endl;
        searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", " W/m3", m.Q_Sensible, 1.0);
        searchMinMax_3D(" max 3D latent heat ", " min 3D latent heat ", " W/m3", m.Q_Latent, 1.0);
        std::cout << std::endl;
    }

    // The lambda and print_heading default exactly as the models' own member versions do, so a
    // row that needs neither reads as five arguments here too.
    void searchMinMax_3D(std::string name_maxValue, std::string name_minValue,
        std::string name_unitValue, Array &value_D, double coeff,
        std::function< double(double) > lambda = [](double i)->double{ return i; },
        bool print_heading = false)
    {
        static const std::string heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
        static const std::string heading_2 = " results based on three dimensional considerations of the problem";

        double maxValue = value_D.x[0][0][0];
        double minValue = value_D.x[0][0][0];
        int imax = 0, jmax = 0, kmax = 0;
        int imin = 0, jmin = 0, kmin = 0;
        for(int j = 0; j < m.jm; j++){
            for(int k = 0; k < m.km; k++){
                for(int i = 0; i < m.im; i++){
                    if(value_D.x[i][j][k] > maxValue){
                        maxValue = value_D.x[i][j][k];
                        imax = i; jmax = j; kmax = k;
                    }else if(value_D.x[i][j][k] < minValue){
                        minValue = value_D.x[i][j][k];
                        imin = i; jmin = j; kmin = k;
                    }
                }
            }
        }
        const int imax_level = imax * (int)m.L_atm/(m.im-1);
        const int imin_level = imin * (int)m.L_atm/(m.im-1);
        std::cout.precision(6);
        if(print_heading){
            std::cout << std::endl << heading_1 << std::endl << heading_2
                      << std::endl << std::endl;
        }
        maxValue = lambda(maxValue * coeff);
        minValue = lambda(minValue * coeff);
        print_pair(name_maxValue, name_minValue, name_unitValue, maxValue, minValue,
                   imax_level, jmax, kmax, imin_level, jmin, kmin);
    }

    void searchMinMax_2D(std::string name_maxValue, std::string name_minValue,
        std::string name_unitValue, Array_2D &value, double coeff)
    {
        double minValue = value.y[0][0];
        double maxValue = value.y[0][0];
        int jmax = 0, kmax = 0, jmin = 0, kmin = 0;
        for(int j = 1; j < m.jm-1; j++){
            for(int k = 1; k < m.km-1; k++){
                if(value.y[j][k] > maxValue){
                    maxValue = value.y[j][k];
                    jmax = j; kmax = k;
                }else if(value.y[j][k] < minValue){
                    minValue = value.y[j][k];
                    jmin = j; kmin = k;
                }
            }
        }
        std::cout.precision(6);
        maxValue = maxValue * coeff;
        minValue = minValue * coeff;
        // A 2D field has no level, and both models print the column as a hard 0.
        print_pair(name_maxValue, name_minValue, name_unitValue, maxValue, minValue,
                   0, jmax, kmax, 0, jmin, kmin);
    }

    /*
     * How far is the run from a steady state, and WHERE.
     *
     * For each prognostic field it reports max|f - f_n| over the grid with the cell it occurs in,
     * f_n being the copy restoreVar() made at the end of the previous iteration — so each number
     * is the largest change that field underwent during ONE iteration. Alongside them, the
     * largest residual of the continuity equation, which is a statement about the pressure
     * projection rather than about convergence in time.
     *
     * IT MUST RUN BEFORE restoreVar(). After it, every f_n equals its f and all thirteen numbers
     * are identically zero.
     *
     * The routine existed in ATJUP, ATSAT, ATNEPT and ATURAN and was called by none of them, so
     * none of this had ever been printed — and in ATSAT it could not have been, because p_dynn
     * was declared and never allocated (a NULL dereference waiting for the first caller). Three
     * further defects, identical in both models:
     *   - the pressure line assigned p_dynn = p_dyn and then differenced against it, so dp was
     *     zero by construction; p_dynn is now allocated and maintained by restoreVar with the
     *     other n-copies, and is back in the restart list;
     *   - the continuity residual compared a magnitude against a signed accumulator;
     *   - min_nh4sh/max_nh4sh were never initialised and case 13 printed stack garbage.
     */
    void steadyQuery(){
        // The field table, in the order the report prints it. Identical in both models, which is
        // why it can live here rather than being asked of the planet.
        struct Entry{
            const char *name;
            Array *f;
            Array *fn;
            double val;
            int i, j, k;
        };
        Entry e[] = {
            {" dp: pressure Poisson equation ",    &m.p_dyn,      &m.p_dynn,      0.0, 0, 0, 0},
            {" du: Navier Stokes equation ",       &m.u,          &m.un,          0.0, 0, 0, 0},
            {" dv: Navier Stokes equation ",       &m.v,          &m.vn,          0.0, 0, 0, 0},
            {" dw: Navier Stokes equation ",       &m.w,          &m.wn,          0.0, 0, 0, 0},
            {" dt: energy transport equation ",    &m.t,          &m.tn,          0.0, 0, 0, 0},
            {" dh2o: h2o vap transport equation ", &m.h2o,        &m.h2on,        0.0, 0, 0, 0},
            {" dh2oc: h2o cl transport equation ", &m.h2o_cloud,  &m.h2o_cloudn,  0.0, 0, 0, 0},
            {" dh2oi: h2o ic transport equation ", &m.h2o_ice,    &m.h2o_icen,    0.0, 0, 0, 0},
            {" dnh3: nh3 transport equation ",     &m.nh3,        &m.nh3n,        0.0, 0, 0, 0},
            {" dnh3: nh3 cl transport equation ",  &m.nh3_cloud,  &m.nh3_cloudn,  0.0, 0, 0, 0},
            {" dnh3i: nh3 i transport equation ",  &m.nh3_ice,    &m.nh3_icen,    0.0, 0, 0, 0},
            {" dnh4sh: nh4sh transport equation ", &m.nh4sh,      &m.nh4shn,      0.0, 0, 0, 0},
        };
        const int ne = (int)(sizeof(e)/sizeof(e[0]));

        // Continuity residual, over the interior only.
        double minimum = 0.0;
        for(int i = 1; i < m.im-1; i++){
            for(int j = 1; j < m.jm-1; j++){
                const double sinthe = sin(m.the.z[j]);
                double costhe = cos(m.the.z[j]);
                if(Planet::costhe_abs() && j > 90) costhe = - costhe;
                const double rm_met = m.metricRadius(m.rad.z[i]);
                const double rmsinthe = rm_met * sinthe;
                for(int k = 1; k < m.km-1; k++){
                    const double dudr = (m.u.x[i+1][j][k] - m.u.x[i-1][j][k])/(2. * m.dr);
                    const double dvdthe = (m.v.x[i][j+1][k] - m.v.x[i][j-1][k])/(2. * m.dthe);
                    const double dwdphi = (m.w.x[i][j][k+1] - m.w.x[i][j][k-1])/(2. * m.dphi);
                    const double residuum = dudr + 2. * m.u.x[i][j][k]/rm_met + dvdthe/rm_met
                        + costhe/rmsinthe * m.v.x[i][j][k] + dwdphi/rmsinthe;
                    // fabs on BOTH sides. It used to compare fabs(residuum) against a `minimum`
                    // holding the SIGNED value, so one negative residual made every later cell
                    // compare true and the reported location became the last cell scanned.
                    if(fabs(residuum) >= minimum){
                        minimum = fabs(residuum);
                        m.i_res = i; m.j_res = j; m.k_res = k;
                    }
                }
            }
        }

        // Per-field change since the previous iteration. >= keeps the LAST cell holding the
        // extreme value, so the traversal order is part of the answer and is left alone.
        for(int i = 0; i < m.im; i++){
            for(int j = 0; j < m.jm; j++){
                for(int k = 0; k < m.km; k++){
                    for(int a = 0; a < ne; a++){
                        const double d = fabs(e[a].f->x[i][j][k] - e[a].fn->x[i][j][k]);
                        if(d >= e[a].val){
                            e[a].val = d;
                            e[a].i = i; e[a].j = j; e[a].k = k;
                        }
                    }
                }
            }
        }

        std::cout.precision(6);
        std::cout.setf(std::ios::fixed);
        std::cout << std::endl << std::endl;
        std::cout << "      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        std::cout << "      3D " << Planet::planet_tag() << " iterational process" << std::endl;
        std::cout << "      max total iteration number nm = " << m.nm << std::endl;
        std::cout << m.steady_iter_line();
        std::cout << std::endl;
        std::cout << std::endl << std::endl << Planet::steady_heading() << std::endl << std::endl;

        // The continuity residual first, then the table — this is the original case order.
        print_steady_line(" residuum: continuity equation ", minimum, m.i_res, m.j_res, m.k_res);
        for(int a = 0; a < ne; a++)
            print_steady_line(e[a].name, e[a].val, e[a].i, e[a].j, e[a].k);

        std::cout << std::endl << std::endl;
        return;
    }

private:

    void print_steady_line(const std::string &name_Value, double Value,
                           int i_loc, int j_loc, int k_loc) const
    {
        // The *1.e-3 treated L_atm as METRES. It is kilometres, and i_loc_level is an int, so
        // every level in this report truncated to 0 and the whole column read "0km" whatever
        // the cell was — which is exactly how it looked as though every field changed most at
        // the bottom boundary. searchMinMax_3D has always had it right.
        const int i_loc_level = i_loc * (int)m.L_atm/(m.im-1);
        int j_loc_deg = 0, k_loc_deg = 0;
        std::string deg_lat, deg_lon;
        if(j_loc <= 90){ j_loc_deg = 90 - j_loc; deg_lat = "°N"; }
        if(j_loc > 90) { j_loc_deg = j_loc - 90; deg_lat = "°S"; }
        if(k_loc <= 180){ k_loc_deg = k_loc;       deg_lon = "°E"; }
        if(k_loc > 180) { k_loc_deg = 360 - k_loc; deg_lon = "°W"; }
        std::cout << std::setiosflags(std::ios::left) << std::setw(36) << std::setfill('.')
            << name_Value << " = " << std::resetiosflags(std::ios::left) << std::setw(12)
            << std::fixed << std::setfill(' ') << Value << std::setw(5) << j_loc_deg
            << std::setw(3) << deg_lat << std::setw(4) << k_loc_deg << std::setw(3) << deg_lon
            << std::setw(6) << i_loc_level << std::setw(2) << "km" << std::endl;
    }
};

#include "cUranusModel.h"
#include "Reporting.h"

using namespace std;

namespace{
    string heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
    string heading_2 = " results based on three dimensional considerations of the problem";
    string level = "km";

    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    };
    HemisphereCoords convert_coords(double lon, double lat){
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
}
/*
*
*/
void cUranusModel::printMinMax(){
    // The 39 fields every model shares are the SHARED Reporting<Planet>; see
    // print_minmax_common() there for why they could not be shared until the units were
    // settled. Below is only what is ATURAN's own.
    Reporting<cUranusModel>(*this).print_minmax_common();

    cout << endl;

    cout << endl << " Velocities " << endl;
    cout << endl;

    cout << endl << " Pressures and Mixture-Density" << endl;
    searchMinMax_3D(" max 3D rho_mix  ", " min 3D rho_mix  ", "kg/m³", rho_mix, 1.0);
    cout << endl;

    cout << endl;

    cout << endl << " Hydrogen Sulfide " << endl;
    searchMinMax_3D(" max 3D h2s_cloud ", " min 3D h2s_cloud ", " kg/m3", h2s_cloud, 1.0);
    searchMinMax_3D(" max 3D h2s_ice ", " min 3D h2s_ice ", " kg/m3", h2s_ice, 1.0);

    // The turbulence closure's own fields. Zero unless ATURAN_TURB is set, and printed regardless
    // so that switching the knob on produces something visible: nue* is what the closure exists to
    // compute, and the question it must answer on Uranus is whether it is larger or smaller than
    // the molecular background 1/re = 1e-3. ATSAT's turned out ~77x SMALLER — the opposite of
    // ATJUP's situation and of what ATJUP's comment claimed — and ATURAN's re is 1000 as well.
    searchMinMax_3D(" max 3D tke ", " min 3D tke ", "/", tke, 1.0);
    searchMinMax_3D(" max 3D dis ", " min 3D dis ", "/", dis, 1.0);
    searchMinMax_3D(" max 3D nue ", " min 3D nue ", "/", nue, 1.0);
    // Per-column friction velocity u_tau, the quantity the closure's k* seed goes as the SQUARE
    // of. Printed because k* seeding to zero is otherwise unattributable.
    searchMinMax_2D(" max 2D vel_star ", " min 2D vel_star ", " m/s", vel_star, 1.0);

    // The precipitation scheme's own fields. Zero unless ATURAN_PRECIP is set, and printed
    // regardless — these arrays were computed and read by nothing, which is the position Q_rad and
    // nue* were in before their rows were added, and it makes "the scheme changes no output"
    // unanswerable. mm/day for the surface maps so they can be read against an energy budget.
    searchMinMax_3D(" max 3D P_rain ", " min 3D P_rain ", " kg/m2/s", P_rain, 1.0);
    searchMinMax_3D(" max 3D P_snow ", " min 3D P_snow ", " kg/m2/s", P_snow, 1.0);
    searchMinMax_3D(" max 3D P_graupel ", " min 3D P_graupel ", " kg/m2/s", P_graupel, 1.0);
    searchMinMax_3D(" max 3D P_nh3_rain ", " min 3D P_nh3_rain ", " kg/m2/s", P_nh3_rain, 1.0);
    searchMinMax_3D(" max 3D P_ch4_rain ", " min 3D P_ch4_rain ", " kg/m2/s", P_ch4_rain, 1.0);
    searchMinMax_3D(" max 3D P_nh4sh ", " min 3D P_nh4sh ", " kg/m2/s", P_nh4sh, 1.0);
    searchMinMax_3D(" max 3D Q_precip ", " min 3D Q_precip ", " W/m3", Q_precip, 1.0);
    searchMinMax_2D(" max 2D precip srf total ", " min 2D precip srf total ", " mm/d", precip_srf_total, 86400.0);
    searchMinMax_2D(" max 2D precip srf H2O ", " min 2D precip srf H2O ", " mm/d", precip_srf_h2o, 86400.0);
    searchMinMax_2D(" max 2D precip srf NH3 ", " min 2D precip srf NH3 ", " mm/d", precip_srf_nh3, 86400.0);
    cout << endl;

    cout << endl;

    // The radiation scheme's own fields, on exactly the argument the precipitation block above
    // makes: zero unless ATURAN_RADIATION is set, printed regardless, so that switching the knob
    // on produces something visible. Radiation.h has been filling these three arrays on this model
    // since the port and NOTHING read them — not this file, not ParaView — which is why every
    // number known about ATURAN's radiation so far came from the one stdout line in the shared
    // header. net radiation is the interface flux the scheme balances; emissivity is the only
    // direct view of the OPACITY, and so the first place a mis-tuned kappa would show.
    cout << endl << " Radiation " << endl;
    searchMinMax_3D(" max 3D net radiation ", " min 3D net radiation ", " W/m2", radiation, 1.0);
    searchMinMax_3D(" max 3D Q_rad ", " min 3D Q_rad ", " W/m3", Q_rad, 1.0);
    searchMinMax_3D(" max 3D emissivity ", " min 3D emissivity ", " /", epsilon, 1.0);

    // Equatorial column profile (j=jm/2, k=km/2), top -> bottom, for a direct check of the
    // radiation / Q_rad fields against the actual T(p). Ported from ATJUP, which was the only
    // model that had it.
    //
    // IT IS THE TABLE THAT WOULD HAVE CAUGHT ATNEPT'S 607x PRESSURE DEFECT ON THE FIRST RUN.
    // Neptune's photosphere was reported at 17.2287 bar and read as an opacity failure for two
    // commits; what was actually wrong was that init_PressureStatic put the TOP of the domain at
    // 15 bar, so tau = 1 was reached in the first layer below the ceiling. A single column of
    // p[bar] from top to bottom shows that immediately, and no amount of column-MEAN diagnostics
    // does. The p[bar] column earns its place here even with the radiation knob off.
    {
        const int j0 = jm / 2, k0 = km / 2;
        cout << endl << " Equatorial column  (j=" << j0 << ", k=" << k0
             << ")   top -> bottom" << endl;
        printf("   %3s  %10s  %8s  %8s  %12s  %14s\n",
               "i", "p[bar]", "T[K]", "eps", "netRad[W/m2]", "Q_rad[W/m3]");
        for(int i = im - 1; i >= 0; i--){
            printf("   %3d  %10.4f  %8.2f  %8.4f  %12.4f  %14.4e\n",
                   i, p_stat.x[i][j0][k0], t.x[i][j0][k0] * t_ref,
                   epsilon.x[i][j0][k0], radiation.x[i][j0][k0], Q_rad.x[i][j0][k0]);
        }
    }

    cout << endl << endl;
}
/*
*
*/
/*
*
*/
// Forwarders to the SHARED Reporting<cUranusModel>. The bodies used to be here in full; the
// machinery is identical in all three models. See Reporting.h.
void cUranusModel::searchMinMax_3D(string name_maxValue, string name_minValue,
    string name_unitValue, Array &value_D, double coeff,
    std::function< double(double) > lambda, bool print_heading){
    Reporting<cUranusModel>(*this).searchMinMax_3D(name_maxValue, name_minValue,
        name_unitValue, value_D, coeff, lambda, print_heading);
}
/*
*
*/
void cUranusModel::searchMinMax_2D(string name_maxValue, string name_minValue,
    string name_unitValue, Array_2D &value, double coeff){
    Reporting<cUranusModel>(*this).searchMinMax_2D(name_maxValue, name_minValue,
        name_unitValue, value, coeff);
}

double cUranusModel::out_maxValue() const{
    return maxValue;
}
/*
*
*/
double cUranusModel::out_minValue() const{
    return minValue;
}
/*
*
*/
void cUranusModel::print_welcome_msg(){
    if(verbose){
        cout << endl << endl << endl;
        cout << "***** Atmosphere Jupiter General Circulation Model(ATUran) applied to laminar flow" << endl;
        cout << "***** program for the computation of Jupiter-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 6 additional transport equations to describe the water vapour, cloud water, cloud ice and nh3 vapour, nh3 cloud and nh3 ice" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonally constant" << endl;
        cout << "***** water and nh3 vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code(German Weather Forecast)" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)" << endl << endl;
        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_welcome_msg_printed = true;
    }
    return;
}
/*
*
*/
void cUranusModel::initMsg(){
    cout << "  present state of the computation " << endl << "  current number of iterations " << endl << endl << "    iter_n = " << iter_n << endl;
    return;
}
/*
*
*/
void cUranusModel::print_final_msg(){
    cout << endl << "***** end of the JupiterAtmosphere General Circulation Modell(ATURAN) *****" << endl << endl;
    if(n == nm)   cout <<  "***** number of artificial time steps      n = " << iter_n << ", end of program reached because of limit of maximum artificial time steps ***** \n\n" << endl;
}
/*
*
*/
/*
*
*/
/*
 * The four defects ATURAN's own copy carried, all of which the shared version does not:
 *   - the continuity residual tested fabs(residuum) but stored the SIGNED value, so one negative
 *     residual made every later cell compare true and the reported location became the last cell
 *     scanned;
 *   - min_nh4sh/max_nh4sh were never initialised, so case 13 read the stack;
 *   - the pressure case was commented out with its `break` left OUTSIDE the comment, so it fell
 *     through and printed case 1's continuity residual under the label "dp: pressure Poisson
 *     equation" — p_dynn did not exist to difference against, which is why it was commented out
 *     in the first place. p_dynn is a real array now, maintained by restoreVar with the other
 *     n-copies;
 *   - i_loc_level multiplied by 1.e-3, treating L_atm as metres, so every level printed 0km.
 * None of this had ever been seen, because nothing called the routine.
 */
void cUranusModel::steadyQuery(){
    Reporting<cUranusModel>(*this).steadyQuery();
}

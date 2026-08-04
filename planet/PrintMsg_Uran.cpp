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

    // The turbulence and precipitation rows the other three models print are NOT here: ATURAN
    // has taken neither template, so tke, dis, nue, vel_star, P_rain, P_snow, P_graupel,
    // P_nh3_rain, P_ch4_rain, P_nh4sh, Q_precip and the three surface precipitation maps do not
    // exist in this model. They belong with those ports, not with the reporting one.
    cout << endl;

    cout << endl;

    cout << endl;

    cout << endl;

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

#include "cUranusModel.h"

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

    cout << endl << endl << " Courant time step   dt = " << dt << endl << endl;

    cout << endl << endl << " Temperatures " << endl;
    searchMinMax_3D(" max 3D temperature ", " min 3D temperature ",
        "degC", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max 3D thermalflux ", " min 3D thermalflux ", "kg/m3", thermalmassflux, r_mix);
    cout << endl;

    cout << endl << " Velocities " << endl;
    searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", "m/s", w, u_0);
    cout << endl;

    cout << endl << " Pressures and Mixture-Density " << endl;
    searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", "mbar", p_dyn, 1e3);
    searchMinMax_3D(" max 3D pressure static  ", " min 3D pressure static  ", "bar", p_stat, 1.0);
    searchMinMax_3D(" max 3D rho_mix  ", " min 3D rho_mix  ", "kg/m³", rho_mix, 1.0);
    cout << endl;

    cout << endl << " Water " << endl;
    searchMinMax_3D(" max 3D h2o       ", " min 3D h2o       ", "kg/m3", h2o,       r_mix);
    searchMinMax_3D(" max 3D h2o_cloud ", " min 3D h2o_cloud ", "kg/m3", h2o_cloud, r_mix);
    searchMinMax_3D(" max 3D h2o_ice   ", " min 3D h2o_ice   ", "kg/m3", h2o_ice,   r_mix);
    cout << endl;

    cout << endl << " Methane " << endl;
    searchMinMax_3D(" max 3D ch4       ", " min 3D ch4       ", "kg/m3", ch4,       r_mix);
    searchMinMax_3D(" max 3D ch4_cloud ", " min 3D ch4_cloud ", "kg/m3", ch4_cloud, r_mix);
    searchMinMax_3D(" max 3D ch4_ice   ", " min 3D ch4_ice   ", "kg/m3", ch4_ice,   r_mix);
    cout << endl;

    cout << endl << " Hydrogen Sulfide " << endl;
    searchMinMax_3D(" max 3D h2s       ", " min 3D h2s       ", "kg/m3",  h2s,         r_mix);
    searchMinMax_3D(" max 3D h2s_cloud ", " min 3D h2s_cloud ", "kg/m3",  h2s_cloud,   r_mix);
    searchMinMax_3D(" max 3D h2s_ice   ", " min 3D h2s_ice   ", "kg/m3",  h2s_ice,     r_mix);
    searchMinMax_3D(" max 3D w_h2s     ", " min 3D w_h2s     ", "kg/m3s", w_h2s,       r_mix);
    searchMinMax_3D(" max 3D j_h2s     ", " min 3D j_h2s     ", "kg/m3",  j_h2s,       r_mix);
    searchMinMax_3D(" max 3D jT_h2s    ", " min 3D jT_h2s    ", "kg/m3",  jT_h2s,      r_mix);
    searchMinMax_3D(" max 3D massflux_h2s  ", " min 3D massflux_h2s  ", "kg/m3", massflux_h2s, r_mix);
    searchMinMax_3D(" max 3D diff_h2s  ", " min 3D diff_h2s  ", "kg/m3",  difflux_h2s, r_mix);
    cout << endl;

    cout << endl << " Ammonia " << endl;
    searchMinMax_3D(" max 3D nh3       ", " min 3D nh3       ", "kg/m3",  nh3,         r_mix);
    searchMinMax_3D(" max 3D nh3_cloud ", " min 3D nh3_cloud ", "kg/m3",  nh3_cloud,   r_mix);
    searchMinMax_3D(" max 3D nh3_ice   ", " min 3D nh3_ice   ", "kg/m3",  nh3_ice,     r_mix);
    searchMinMax_3D(" max 3D w_nh3     ", " min 3D w_nh3     ", "kg/m3s", w_nh3,       r_mix);
    searchMinMax_3D(" max 3D j_nh3     ", " min 3D j_nh3     ", "kg/m3",  j_nh3,       r_mix);
    searchMinMax_3D(" max 3D jT_nh3    ", " min 3D jT_nh3    ", "kg/m3",  jT_nh3,      r_mix);
    searchMinMax_3D(" max 3D massflux_nh3  ", " min 3D massflux_nh3  ", "kg/m3", massflux_nh3, r_mix);
    searchMinMax_3D(" max 3D diff_nh3  ", " min 3D diff_nh3  ", "kg/m3",  difflux_nh3, r_mix);
    cout << endl;

    cout << endl << " Ammonia Hydrosulfide " << endl;
    searchMinMax_3D(" max 3D nh4sh         ", " min 3D nh4sh         ", "mg/m3",  nh4sh,          1e6 * r_mix);
    searchMinMax_3D(" max 3D w_nh4sh       ", " min 3D w_nh4sh       ", "mg/m3s", w_nh4sh,        1e6 * r_mix);
    searchMinMax_3D(" max 3D massflux_nh4sh", " min 3D massflux_nh4sh", "mg/m3",  massflux_nh4sh, 1e6 * r_mix);
    searchMinMax_3D(" max 3D j_nh4sh       ", " min 3D j_nh4sh       ", "kg/m3",  j_nh4sh,        1e6 * r_mix);
    searchMinMax_3D(" max 3D jT_nh4sh      ", " min 3D jT_nh4sh      ", "kg/m3",  jT_nh4sh,       1e6 * r_mix);
    searchMinMax_3D(" max 3D diff_nh4sh    ", " min 3D diff_nh4sh    ", "mg/m3",  difflux_nh4sh,  1e6 * r_mix);
    cout << endl;

    cout << endl << " Forces " << endl;
    searchMinMax_3D(" max 3D Coriolis force  ", " min 3D Coriolis force  ", "mN/m3",  CoriolisForce,   1e3);
    searchMinMax_3D(" max 3D centrifugal frc ", " min 3D centrifugal frc ", "mN/m3", CentrifugalForce, 1e3);
    searchMinMax_3D(" max 3D buoyancy force  ", " min 3D buoyancy force  ", "N/m3",  BuoyancyForce,   1.0);
    searchMinMax_3D(" max 3D presgrad force  ", " min 3D presgrad force  ", "N/m3",  PresGradForce,   1.0);
    cout << endl;

    cout << endl << " Energies " << endl;
    searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", "W/m3", Q_Sensible, 1.0);
    searchMinMax_3D(" max 3D latent heat   ", " min 3D latent heat   ", "W/m3", Q_Latent,   1.0);
    cout << endl << endl;
}
/*
*
*/
void cUranusModel::searchMinMax_3D(string name_maxValue, string name_minValue,
    string name_unitValue, Array &value_D, double coeff,
    std::function< double(double) > lambda, bool print_heading){
    double maxValue = value_D.x[0][0][0];
    double minValue = value_D.x[0][0][0];
    int imax = 0;
    int jmax = 0;
    int kmax = 0;
    int imin = 0;
    int jmin = 0;
    int kmin = 0;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(value_D.x[i][j][k] > maxValue){
                    maxValue = value_D.x[i][j][k];
                    imax = i;
                    jmax = j;
                    kmax = k;
                }else if(value_D.x[i][j][k] < minValue){
                    minValue = value_D.x[i][j][k];
                    imin = i;
                    jmin = j;
                    kmin = k;
                }
            }
        }
    }
    int imax_level = imax * (int)L_atm/(im-1);
    int imin_level = imin * (int)L_atm/(im-1);
    //  maximum latitude and longitude units recalculated
    HemisphereCoords coords = convert_coords(kmax, jmax);
    int jmax_deg = coords.lat;
    string deg_lat_max = coords.north_or_south;
    int kmax_deg = coords.lon;
    string deg_lon_max = coords.east_or_west;
    //  minimum latitude and longitude units recalculated
    coords = convert_coords(kmin, jmin);
    int jmin_deg = coords.lat;
    string deg_lat_min= coords.north_or_south;
    int kmin_deg = coords.lon;
    string deg_lon_min = coords.east_or_west;
    cout.precision(6);
    if(print_heading){
        cout << endl << heading_1 << endl << heading_2 << endl << endl;
    }
    maxValue = lambda(maxValue * coeff);
    minValue = lambda(minValue * coeff);
    cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " << 
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(6) << 
        name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg << 
        setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "          " << 
        setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<< 
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(6) << 
        name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg << 
        setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
}
/*
*
*/
void cUranusModel::searchMinMax_2D(string name_maxValue, string name_minValue, 
    string name_unitValue, Array_2D &value, double coeff){
    double minValue = value.y[0][0];
    double maxValue = value.y[0][0];
    int jmax = 0;
    int kmax = 0;
    int jmin = 0;
    int kmin = 0;  
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
            if(value.y[j][k] > maxValue){
                maxValue = value.y[j][k];
                jmax = j;
                kmax = k;
            }else if(value.y[j][k] < minValue){
                minValue = value.y[j][k];
                jmin = j;
                kmin = k;
            }
        }
    }
    int imax_level = 0;
    int imin_level = 0;
    //  maximum latitude and longitude units recalculated
    HemisphereCoords coords = convert_coords(kmax, jmax);
    int jmax_deg = coords.lat;
    string deg_lat_max = coords.north_or_south;
    int kmax_deg = coords.lon;
    string deg_lon_max = coords.east_or_west;
    //  minimum latitude and longitude units recalculated
    coords = convert_coords(kmin, jmin);
    int jmin_deg = coords.lat;
    string deg_lat_min= coords.north_or_south;
    int kmin_deg = coords.lon;
    string deg_lon_min = coords.east_or_west;
    cout.precision(6);
    maxValue = maxValue * coeff;
    minValue = minValue * coeff;
    cout << setiosflags(ios::left) << setw(26) << setfill('.') << name_maxValue << " = " << 
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << maxValue << setw(6) << 
        name_unitValue << setw(5) << jmax_deg << setw(3) << deg_lat_max << setw(4) << kmax_deg << 
        setw(3) << deg_lon_max << setw(6) << imax_level << setw(2) << level << "          " << 
        setiosflags(ios::left) << setw(26) << setfill('.') << name_minValue << " = "<< 
        resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << minValue << setw(6) << 
        name_unitValue << setw(5)  << jmin_deg << setw(3) << deg_lat_min << setw(4) << kmin_deg << 
        setw(3) << deg_lon_min  << setw(6) << imin_level << setw(2) << level << endl;
}
/*
*
*/
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
void cUranusModel::steadyQuery(){
    int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, 
        j_c, k_c, i_cloud, j_cloud, k_cloud, i_ice, j_ice, k_ice, i_nh3, 
        j_nh3, k_nh3, i_nh3_cloud, j_nh3_cloud, k_nh3_cloud, i_nh3_ice, 
        j_nh3_ice, k_nh3_ice, i_nh4sh, j_nh4sh, k_nh4sh, i_p, j_p, k_p;
    int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;
    double max_u, max_v, max_w, max_t, max_c, max_cloud, max_ice, max_p, 
        max_nh3, max_nh3_cloud, max_nh3_ice, max_nh4sh;
    double min_u, min_v, min_w, min_t, min_c, min_cloud, min_ice, 
        min_nh3, min_nh3_cloud, min_nh3_ice, min_nh4sh, min_p;
    double Value;
    string name_Value;
    string level, deg_north, deg_south, deg_west, deg_east, deg_lat, 
        deg_lon, heading;
    min_u = max_u = 0.;
    min_v = max_v = 0.;
    min_w = max_w = 0.;
    min_t = max_t = 0.;
    min_c = max_c = 0.;
    min_p = max_p = 0.;
    min_cloud = max_cloud = 0.;
    min_ice = max_ice = 0.;
    min_nh3 = max_nh3 = 0.;
    min_nh3_cloud = max_nh3_cloud = 0.;
    min_nh3_ice = max_nh3_ice = 0.;
    double sinthe = 0., costhe = 0., rmsinthe = 0.;
    double dudr = 0., dvdthe = 0., dwdphi = 0.;
    double residuum = 0.;
    double minimum = 0.;
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            rmsinthe = rad.z[i] * sinthe;
            for(int k = 1; k < km-1; k++){
                dudr = (u.x[i+1][j][k] - u.x[i-1][j][k])/(2. * dr);
                dvdthe = (v.x[i][j+1][k] - v.x[i][j-1][k])/(2. * dthe);
                dwdphi = (w.x[i][j][k+1] - w.x[i][j][k-1])/(2. * dphi);
                residuum = dudr + 2. * u.x[i][j][k]/rad.z[i] + dvdthe/rad.z[i]
                    + costhe/rmsinthe * v.x[i][j][k] + dwdphi/rmsinthe;
                if(fabs(residuum) >= minimum){
                    minimum = residuum;
                    i_res = i;
                    j_res = j;
                    k_res = k;
                }
            }
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                p_dynn.x[i][j][k] = p_dyn.x[i][j][k];
                max_p = fabs(p_dyn.x[i][j][k] - p_dynn.x[i][j][k]);
                if(max_p >= min_p){
                    min_p = max_p;
                    i_p = i;
                    j_p = j;
                    k_p = k;
                }
                max_u = fabs(u.x[i][j][k] - un.x[i][j][k]);
                if(max_u >= min_u){
                    min_u = max_u;
                    i_u = i;
                    j_u = j;
                    k_u = k;
                }
                max_v = fabs(v.x[i][j][k] - vn.x[i][j][k]);
                if(max_v >= min_v){
                    min_v = max_v;
                    i_v = i;
                    j_v = j;
                    k_v = k;
                }
                max_w = fabs(w.x[i][j][k] - wn.x[i][j][k]);
                if(max_w >= min_w){
                    min_w = max_w;
                    i_w = i;
                    j_w = j;
                    k_w = k;
                }
                max_t = fabs(t.x[i][j][k] - tn.x[i][j][k]);
                if(max_t >= min_t){
                    min_t = max_t;
                    i_t = i;
                    j_t = j;
                    k_t = k;
                }
                max_c = fabs(ch4.x[i][j][k] - ch4n.x[i][j][k]);
                if(max_c >= min_c){
                    min_c = max_c;
                    i_c = i;
                    j_c = j;
                    k_c = k;
                }
                max_cloud = fabs(ch4_cloud.x[i][j][k] - ch4_cloudn.x[i][j][k]);
                if(max_cloud >= min_cloud){
                    min_cloud = max_cloud;
                    i_cloud = i;
                    j_cloud = j;
                    k_cloud = k;
                }
                max_ice = fabs(ch4_ice.x[i][j][k] - ch4_icen.x[i][j][k]);
                if(max_ice >= min_ice){
                    min_ice = max_ice;
                    i_ice = i;
                    j_ice = j;
                    k_ice = k;
                }
                max_c = fabs(h2o.x[i][j][k] - h2on.x[i][j][k]);
                if(max_c >= min_c){
                    min_c = max_c;
                    i_c = i;
                    j_c = j;
                    k_c = k;
                }
                max_cloud = fabs(h2o_cloud.x[i][j][k] - h2o_cloudn.x[i][j][k]);
                if(max_cloud >= min_cloud){
                    min_cloud = max_cloud;
                    i_cloud = i;
                    j_cloud = j;
                    k_cloud = k;
                }
                max_ice = fabs(h2o_ice.x[i][j][k] - h2o_icen.x[i][j][k]);
                if(max_ice >= min_ice){
                    min_ice = max_ice;
                    i_ice = i;
                    j_ice = j;
                    k_ice = k;
                }
                max_nh3 = fabs(nh3.x[i][j][k] - nh3n.x[i][j][k]);
                if(max_nh3 >= min_nh3){
                    min_nh3= max_nh3;
                    i_nh3 = i;
                    j_nh3 = j;
                    k_nh3 = k;
                }
                max_nh3_cloud = fabs(nh3_cloud.x[i][j][k] - nh3_cloudn.x[i][j][k]);
                if(max_nh3_cloud >= min_nh3_cloud){
                    min_nh3_cloud= max_nh3_cloud;
                    i_nh3_cloud = i;
                    j_nh3_cloud = j;
                    k_nh3_cloud = k;
                }
                max_nh3_ice = fabs(nh3_ice.x[i][j][k] - nh3_icen.x[i][j][k]);
                if(max_nh3_ice >= min_nh3_ice){
                    min_nh3_ice= max_nh3_ice;
                    i_nh3_ice = i;
                    j_nh3_ice = j;
                    k_nh3_ice = k;
                }
                max_nh4sh = fabs(nh4sh.x[i][j][k] - nh4shn.x[i][j][k]);
                if(max_nh4sh >= min_nh4sh){
                    min_nh4sh= max_nh4sh;
                    i_nh4sh = i;
                    j_nh4sh = j;
                    k_nh4sh = k;
                }
            }
        }
    }
cout << "   " << i_p << "   " << j_p << "   " << k_p << "   " << max_p << "   " << min_p << endl;
cout << "   " << i_t << "   " << j_t << "   " << k_t << "   " << max_t << "   " << min_t << endl;
    cout.precision(6);
    cout.setf(ios::fixed);
    cout << endl << endl;
    cout << "      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    cout << "      3D ATURAN iterational process" << endl;
    cout << "      max total iteration number nm = " << nm << endl;
    cout << "      n = " << n << "      velocity_iter = " << velocity_iter << "     " << "pressure_iter = " << pressure_iter << endl;
    cout << endl;
    heading = " 3D iterational process for the surface boundary conditions\n printout of maximum and minimum absolute and relative errors of the computed values at their locations: level, latitude, longitude";
    cout << endl << endl << heading << endl << endl;
    level = "km";
    deg_north = "°N";
    deg_south = "°S";
    deg_west = "°W";
    deg_east = "°E";
    int choice = {1};
    preparation:
    switch(choice){
        case 1 :    name_Value = " residuum: continuity equation ";
                        Value = minimum;
                        i_loc = i_res;
                        j_loc = j_res;
                        k_loc = k_res;
                        break;
        case 2 :    name_Value = " dp: pressure Poisson equation ";

                        Value = min_p;
                        i_loc = i_p;
                        j_loc = j_p;
                        k_loc = k_p;
                        break;
        case 3 :    name_Value = " du: Navier Stokes equation ";
                        Value = min_u;
                        i_loc = i_u;
                        j_loc = j_u;
                        k_loc = k_u;
                        break;
        case 4 :    name_Value = " dv: Navier Stokes equation ";
                        Value = min_v;
                        i_loc = i_v;
                        j_loc = j_v;
                        k_loc = k_v;
                        break;
        case 5 :    name_Value = " dw: Navier Stokes equation ";
                        Value = min_w;
                        i_loc = i_w;
                        j_loc = j_w;
                        k_loc = k_w;
                        break;
        case 6 :    name_Value = " dt: energy transport equation ";
                        Value = min_t;
                        i_loc = i_t;
                        j_loc = j_t;
                        k_loc = k_t;
                        break;
        case 7 :    name_Value = " dh2o: h2o vap transport equation ";
                        Value = min_c;
                        i_loc = i_c;
                        j_loc = j_c;
                        k_loc = k_c;
                        break;
        case 8 :    name_Value = " dh2oc: h2o cl transport equation ";
                        Value = min_cloud;
                        i_loc = i_cloud;
                        j_loc = j_cloud;
                        k_loc = k_cloud;
                        break;
        case 9 :    name_Value = " dh2oi: h2o ic transport equation ";
                        Value = min_ice;
                        i_loc = i_ice;
                        j_loc = j_ice;
                        k_loc = k_ice;
                        break;
        case 10 :    name_Value = " dnh3: nh3 transport equation ";
                        Value = min_nh3;
                        i_loc = i_nh3;
                        j_loc = j_nh3;
                        k_loc = k_nh3;
                        break;
        case 11 :    name_Value = " dnh3: nh3 cl transport equation ";
                        Value = min_nh3_cloud;
                        i_loc = i_nh3_cloud;
                        j_loc = j_nh3_cloud;
                        k_loc = k_nh3_cloud;
                        break;
        case 12 :    name_Value = " dnh3i: nh3 i transport equation ";
                        Value = min_nh3_ice;
                        i_loc = i_nh3_ice;
                        j_loc = j_nh3_ice;
                        k_loc = k_nh3_ice;
                        break;
        case 13 :    name_Value = " dnh4sh: nh4sh transport equation ";
                        Value = min_nh4sh;
                        i_loc = i_nh4sh;
                        j_loc = j_nh4sh;
                        k_loc = k_nh4sh;
                        break;
        default :     cout << choice << "error in iterationPrintout_3D member function in class Accuracy" << endl;
    }
    i_loc_level = i_loc*int(L_atm)/(im-1)*1.e-3;
    if(j_loc <= 90){
        j_loc_deg = 90 - j_loc;
        deg_lat = deg_north;
    }
    if(j_loc > 90){
        j_loc_deg = j_loc - 90;
        deg_lat = deg_south;
    }
    if(k_loc <= 180){
        k_loc_deg = k_loc;
        deg_lon = deg_east;
    }
    if(k_loc > 180){
        k_loc_deg = 360 - k_loc;
        deg_lon = deg_west;
    }
    cout << setiosflags(ios::left) << setw(36) << setfill('.') << name_Value << " = " << resetiosflags(ios::left) << setw(12) << fixed << setfill(' ') << Value << setw(5) << j_loc_deg << setw(3) << deg_lat << setw(4) << k_loc_deg << setw(3) << deg_lon << setw(6) << i_loc_level << setw(2) << level << endl;
    choice++;
    if(choice <= 13) goto preparation;
    cout << endl << endl;
    return;
}


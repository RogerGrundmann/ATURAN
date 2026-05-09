/*
 * Atmosphere General Circulation Modell(ATURAN) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cUranusModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

void cUranusModel::writeData(){
    cout << endl << "      ATURAN: writeData" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

//    int i_radial = 0;
//    int i_radial = 3;
    int i_radial = 20;
//    int i_radial = 36;
//    int i_radial = 39;
//    int i_radial = 40;
    paraview_vtk_radial(iter_n, i_radial);
    int j_longal = 112;
    paraview_vtk_longal(iter_n, j_longal);
    int k_zonal = 180;
//    int k_zonal = 0;
    paraview_vtk_zonal(iter_n, k_zonal);

    if(paraview_flag && (iter_n % panorama_print == 0)){
        paraview_panorama_vts(iter_n);
//        paraview_sphere_vts(iter_n);
    }

    UranusPlotData();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for writeData\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: writeData ended" << endl;

    return;
}
/*
*
*/
void cUranusModel::writeResults(){
    cout << endl << "      ATURAN: writeResults" << endl;
//    double coeff_mmWS = r_mix/r_h2o;  // coeff_mmWS = 1.2041/0.0094 [kg/m³/kg/m³] = 128,0827 [/]
//    double coeff_lv = lv_h2o /(cp_h2o * t_0);  // coefficient for the specific latent Evaporation heat(Condensation heat), coeff_lv = 9.1069 in [/]
//    double coeff_ls = ls /(cp_h2o * t_0);  // coefficient for the specific latent Evaporation heat(Condensation heat), coeff_ls = 10.3091 in [/]
//    double a, e;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            precipitable_water.y[j][k] = 0.;  // precipitable water
        }
    }
    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Q_Latent.x[0][j][k] = c43 * Q_Latent.x[1][j][k] - c13 * Q_Latent.x[2][j][k];
            Q_Latent.x[im-1][j][k] = c43 * Q_Latent.x[im-2][j][k] - c13 * Q_Latent.x[im-3][j][k];
            Q_Sensible.x[0][j][k] = c43 * Q_Sensible.x[1][j][k] - c13 * Q_Sensible.x[2][j][k];
            Q_Sensible.x[im-1][j][k] = c43 * Q_Sensible.x[im-2][j][k] - c13 * Q_Sensible.x[im-3][j][k];
            BuoyancyForce.x[0][j][k] = c43 * BuoyancyForce.x[1][j][k] - c13 * BuoyancyForce.x[2][j][k];
            BuoyancyForce.x[im-1][j][k] = c43 * BuoyancyForce.x[im-2][j][k] - c13 * BuoyancyForce.x[im-3][j][k];
/*
            nh3.x[0][j][k] = c43 * nh3.x[1][j][k] - c13 * nh3.x[2][j][k];
            nh3.x[im-1][j][k] = c43 * nh3.x[im-2][j][k] - c13 * nh3.x[im-3][j][k];
            nh3_cloud.x[0][j][k] = c43 * nh3_cloud.x[1][j][k] - c13 * nh3_cloud.x[2][j][k];
            nh3_cloud.x[im-1][j][k] = c43 * nh3_cloud.x[im-2][j][k] - c13 * nh3_cloud.x[im-3][j][k];
            nh3_ice.x[0][j][k] = c43 * nh3_ice.x[1][j][k] - c13 * nh3_ice.x[2][j][k];
            nh3_ice.x[im-1][j][k] = c43 * nh3_ice.x[im-2][j][k] - c13 * nh3_ice.x[im-3][j][k];
*/
        }
    }
    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            t.x[i][0][k] = c43 * t.x[i][1][k] - c13 * t.x[i][2][k];
            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] - c13 * t.x[i][jm-3][k];
            Q_Latent.x[i][0][k] = c43 * Q_Latent.x[i][1][k] - c13 * Q_Latent.x[i][2][k];
            Q_Latent.x[i][jm-1][k] = c43 * Q_Latent.x[i][jm-2][k] - c13 * Q_Latent.x[i][jm-3][k];
            Q_Sensible.x[i][0][k] = c43 * Q_Sensible.x[i][1][k] - c13 * Q_Sensible.x[i][2][k];
            Q_Sensible.x[i][jm-1][k] = c43 * Q_Sensible.x[i][jm-2][k] - c13 * Q_Sensible.x[i][jm-3][k];
            BuoyancyForce.x[i][0][k] = c43 * BuoyancyForce.x[i][1][k] - c13 * BuoyancyForce.x[i][2][k];
            BuoyancyForce.x[i][jm-1][k] = c43 * BuoyancyForce.x[i][jm-2][k] - c13 * BuoyancyForce.x[i][jm-3][k];
/*
            nh3.x[i][0][k] = c43 * nh3.x[i][1][k] - c13 * nh3.x[i][2][k];
            nh3.x[i][jm-1][k] = c43 * nh3.x[i][jm-2][k] - c13 * nh3.x[i][jm-3][k];
            nh3_cloud.x[i][0][k] = c43 * nh3_cloud.x[i][1][k] - c13 * nh3_cloud.x[i][2][k];
            nh3_cloud.x[i][jm-1][k] = c43 * nh3_cloud.x[i][jm-2][k] - c13 * nh3_cloud.x[i][jm-3][k];
            nh3_ice.x[i][0][k] = c43 * nh3_ice.x[i][1][k] - c13 * nh3_ice.x[i][2][k];
            nh3_ice.x[i][jm-1][k] = c43 * nh3_ice.x[i][jm-2][k] - c13 * nh3_ice.x[i][jm-3][k];
            h2o.x[i][0][k] = c43 * h2o.x[i][1][k] - c13 * h2o.x[i][2][k];
            h2o.x[i][jm-1][k] = c43 * h2o.x[i][jm-2][k] - c13 * h2o.x[i][jm-3][k];
            h2o_cloud.x[i][0][k] = c43 * h2o_cloud.x[i][1][k] - c13 * h2o_cloud.x[i][2][k];
            h2o_cloud.x[i][jm-1][k] = c43 * h2o_cloud.x[i][jm-2][k] - c13 * h2o_cloud.x[i][jm-3][k];
            h2o_ice.x[i][0][k] = c43 * h2o_ice.x[i][1][k] - c13 * h2o_ice.x[i][2][k];
            h2o_ice.x[i][jm-1][k] = c43 * h2o_ice.x[i][jm-2][k] - c13 * h2o_ice.x[i][jm-3][k];
*/
/*
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] - c13 * P_snow.x[i][jm-3][k];
*/
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] - c13 * t.x[i][j][2];
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] - c13 * t.x[i][j][km-3];
            t.x[i][j][0] = t.x[i][j][km-1] =(t.x[i][j][0] + t.x[i][j][km-1])/ 2.;
            Q_Latent.x[i][j][0] = c43 * Q_Latent.x[i][j][1] - c13 * Q_Latent.x[i][j][2];
            Q_Latent.x[i][j][km-1] = c43 * Q_Latent.x[i][j][km-2] - c13 * Q_Latent.x[i][j][km-3];
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] =(Q_Latent.x[i][j][0] + Q_Latent.x[i][j][km-1])/ 2.;
            Q_Sensible.x[i][j][0] = c43 * Q_Sensible.x[i][j][1] - c13 * Q_Sensible.x[i][j][2];
            Q_Sensible.x[i][j][km-1] = c43 * Q_Sensible.x[i][j][km-2] - c13 * Q_Sensible.x[i][j][km-3];
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] =(Q_Sensible.x[i][j][0] + Q_Sensible.x[i][j][km-1])/ 2.;
            BuoyancyForce.x[i][j][0] = c43 * BuoyancyForce.x[i][j][1] - c13 * BuoyancyForce.x[i][j][2];
            BuoyancyForce.x[i][j][km-1] = c43 * BuoyancyForce.x[i][j][km-2] - c13 * BuoyancyForce.x[i][j][km-3];
            BuoyancyForce.x[i][j][0] = BuoyancyForce.x[i][j][km-1] =(BuoyancyForce.x[i][j][0] + BuoyancyForce.x[i][j][km-1])/ 2.;
/*
            nh3.x[i][j][0] = c43 * nh3.x[i][j][1] - c13 * nh3.x[i][j][2];
            nh3.x[i][j][km-1] = c43 * nh3.x[i][j][km-2] - c13 * nh3.x[i][j][km-3];
            nh3.x[i][j][0] = nh3.x[i][j][km-1] =(nh3.x[i][j][0] + nh3.x[i][j][km-1])/ 2.;
            nh3_cloud.x[i][j][0] = c43 * nh3_cloud.x[i][j][1] - c13 * nh3_cloud.x[i][j][2];
            nh3_cloud.x[i][j][km-1] = c43 * nh3_cloud.x[i][j][km-2] - c13 * nh3_cloud.x[i][j][km-3];
            nh3_cloud.x[i][j][0] = nh3_cloud.x[i][j][km-1] =(nh3_cloud.x[i][j][0] + nh3_cloud.x[i][j][km-1])/ 2.;
            nh3_ice.x[i][j][0] = c43 * nh3_ice.x[i][j][1] - c13 * nh3_ice.x[i][j][2];
            nh3_ice.x[i][j][km-1] = c43 * nh3_ice.x[i][j][km-2] - c13 * nh3_ice.x[i][j][km-3];
            nh3_ice.x[i][j][0] = nh3_ice.x[i][j][km-1] =(nh3_ice.x[i][j][0] + nh3_ice.x[i][j][km-1])/ 2.;
            h2o_cloud.x[i][j][0] = c43 * h2o_cloud.x[i][j][1] - c13 * h2o_cloud.x[i][j][2];
            h2o_cloud.x[i][j][km-1] = c43 * h2o_cloud.x[i][j][km-2] - c13 * h2o_cloud.x[i][j][km-3];
            h2o_cloud.x[i][j][0] = h2o_cloud.x[i][j][km-1] =(h2o_cloud.x[i][j][0] + h2o_cloud.x[i][j][km-1])/ 2.;
            h2o_ice.x[i][j][0] = c43 * h2o_ice.x[i][j][1] - c13 * h2o_ice.x[i][j][2];
            h2o_ice.x[i][j][km-1] = c43 * h2o_ice.x[i][j][km-2] - c13 * h2o_ice.x[i][j][km-3];
            h2o_ice.x[i][j][0] = h2o_ice.x[i][j][km-1] =(h2o_ice.x[i][j][0] + h2o_ice.x[i][j][km-1])/ 2.;
*/
/*
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] - c13 * P_rain.x[i][j][2];
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] - c13 * P_rain.x[i][j][km-3];
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] =(P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/ 2.;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] - c13 * P_snow.x[i][j][2];
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] - c13 * P_snow.x[i][j][km-3];
            P_snow.x[i][j][0] = P_snow.x[i][j][km-1] =(P_snow.x[i][j][0] + P_snow.x[i][j][km-1])/ 2.;
*/
        }
    }
/*
    precipitable_water.y[0][0] = 0.;
    precipitable_water.y[0][0] = 0.;
    nh3_total.y[0][0] = 0.;
    nh3_cloud_total.y[0][0] = 0.;
    nh3_ice_total.y[0][0] = 0.;
    double precipitablewater_average = 0.;
    double precipitation_average = 0.;
    double h2_vegetation_average = 0.;
    double he_vegetation_average = 0.;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            nh3_total.y[j][k] = nh3.x[0][j][k];
            nh3_cloud_total.y[j][k] = nh3_cloud.x[0][j][k];
            nh3_ice_total.y[j][k] = nh3_ice.x[0][j][k];
            for(int i = 0; i < im; i++){
                e = h2o.x[i][j][k] * p_stat.x[i][j][k]/ep_h2o;  // water vapour pressure in hPa
                a = 216.6 * e /(t.x[i][j][k] * t_0);  // absolute humidity in kg/m3
                precipitable_water.y[j][k] += a * L_atm /(double)(im - 1);  // kg/m³ * m
            }
        }
    }
    double coeff_prec = 86400.;  // dimensions see below
// surface values of precipitation and precipitable water
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Precipitation.y[j][k] = coeff_prec * (P_rain.x[0][j][k] + P_snow.x[0][j][k]);
            // 60 s * 60 min * 24 h = 86400 s == 1 d
            // Precipitation, P_rain and P_snow in kg/ ( m² * s ) = mm/s
            // Precipitation in 86400. * kg/ ( m² * d ) = 86400 mm/d
            // kg/ ( m² * s ) == mm/s ( Kraus, p. 94 )
            if(Precipitation.y[j][k] >= 25.)  Precipitation.y[j][k] = 25.;
            if(Precipitation.y[j][k] <= 0)  Precipitation.y[j][k] = 0.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            precipitablewater_average += precipitable_water.y[j][k];
            precipitation_average += Precipitation.y[j][k];
            h2_vegetation_average += nh3_total.y[j][k];
            he_vegetation_average += nh3_cloud_total.y[j][k];
        }
    }
    h2_vegetation_average = h2_vegetation_average /(double)((jm-1)*(km-1));
    he_vegetation_average = he_vegetation_average /(double)((jm-1)*(km-1));
    precipitablewater_average = precipitablewater_average /(double)((jm-1)*(km-1));
    precipitation_average = 365. * precipitation_average /(double)((jm-1)*(km-1));
    cout.precision(2);
    string level = "m";
    string deg_north = "°N";
    string deg_south = "°S";
    string deg_west = "°W";
    string deg_east = "°E";
    string name_Value_7 = " precipitable water average ";
    string name_Value_8 = " precipitation average per year ";
    string name_Value_9 = " precipitation average per day ";
    string name_Value_22 = " H2_average ";
    string name_Value_25 = " He_average ";
    string name_unit_mmd = " mm/d";
    string name_unit_mm = " mm";
    string name_unit_mma = " mm/a";
    string name_unit_ppm = " kg/kg";
    cout << endl;
    double Value_7 = precipitablewater_average;
    double Value_8 = precipitation_average;
    cout << setw(6)<< setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_7 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_7 << setw(6)<< name_unit_mm << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_8 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_8 << setw(6)<< name_unit_mma << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_9 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_8/365. << setw(6)<< name_unit_mmd << endl;
    double Value_9 = h2_vegetation_average;
    double Value_25 = he_vegetation_average;
    cout << setw(6)<< setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_22 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_9 << setw(6)<< name_unit_ppm << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_25 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_25 << setw(6)<< name_unit_ppm << endl << endl << endl;
    return;
*/

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for writeResults\n", elapsed.count() * 1e-9);

    cout << "      ATURAN: writeResults ended" << endl;

}
/*
*
*/

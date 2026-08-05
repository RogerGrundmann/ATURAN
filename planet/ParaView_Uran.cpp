/*
 * Uranus General Circulation Modell(ATURAN)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include "cUranusModel.h"
#include "ParaViewWriter.h"

using namespace std;
//using namespace AtomUtils;

// The five dumpers moved to the SHARED ParaViewWriter.h, as namespace ParaViewIO. ATURAN's
// copies (namespace ParaViewUranus) were IDENTICAL code to ATSAT's, ATJUP's and ATNEPT's — the
// whole difference was `for (` against `for(` and the /* */ separators between them. Four models
// independently carried the same file-format code.
/*
 *
*/
void cUranusModel::paraview_panorama_vts(int n){
    using namespace ParaViewIO;
    // Header, coordinates, the Velocity array and the Temperature array are the
    // SHARED ParaViewWriter.h. What stays here is the field list below and the
    // scalars string that has to agree with it.
    // Display units, settled across all four models: species and their fluxes in g/m3 (1e3), the
    // NH4SH family in ug/m3 (1e9), thermalmassflux in W/m3 (1.0). These are the same units the
    // shared Reporting.h prints, so an array name means one thing in the .vtk and in the log.
    // See the note above the species block in Reporting.h for why NH4SH needs its own.
    ParaViewWriter<cUranusModel> pv(*this);
    ofstream Uranus_panorama_vts_File = pv.open_panorama(n,
        "Temperature PressureDynamic PressureStatic NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce Q_Latent Q_Sensible BuoyancyForce Q_rad_mW_m3 Radiation ");
    pv.panorama_velocity(Uranus_panorama_vts_File);
    pv.panorama_temperature(Uranus_panorama_vts_File);
    dump_array("u-component", u, u_0, Uranus_panorama_vts_File);
    dump_array("v-component", v, u_0, Uranus_panorama_vts_File);
    dump_array("w-component", w, u_0, Uranus_panorama_vts_File);

    dump_array("PressureDyn", p_dyn, 1e3, Uranus_panorama_vts_File);
    dump_array("PressureStat", p_stat, 1.0, Uranus_panorama_vts_File);
    dump_array("rho_mix", rho_mix, 1.0, Uranus_panorama_vts_File);

    dump_array("H2O", h2o, 1e3, Uranus_panorama_vts_File);
    dump_array("H2OCloud", h2o_cloud, 1e3, Uranus_panorama_vts_File);
    dump_array("H2OIce", h2o_ice, 1e3, Uranus_panorama_vts_File);

    dump_array("CH4", ch4, 1e3, Uranus_panorama_vts_File);
    dump_array("CH4Cloud", ch4_cloud, 1e3, Uranus_panorama_vts_File);
    dump_array("CH4Ice", ch4_ice, 1e3, Uranus_panorama_vts_File);

    dump_array("H2S", h2s, 1e3, Uranus_panorama_vts_File);
    dump_array("H2SCloud", h2s_cloud, 1e3, Uranus_panorama_vts_File);
    dump_array("H2SIce", h2s_ice, 1e3, Uranus_panorama_vts_File);
    dump_array("w_h2s", w_h2s, 1e3, Uranus_panorama_vts_File);
//    dump_array("j_h2s", j_h2s, 1e3, Uranus_panorama_vts_File);
//    dump_array("jT_h2s", jT_h2s, 1e3, Uranus_panorama_vts_File);

    dump_array("NH3", nh3, 1e3, Uranus_panorama_vts_File);
    dump_array("NH3Cloud", nh3_cloud, 1e3, Uranus_panorama_vts_File);
    dump_array("NH3Ice", nh3_ice, 1e3, Uranus_panorama_vts_File);
    dump_array("w_nh3", w_nh3, 1e3, Uranus_panorama_vts_File);
//    dump_array("j_nh3", j_nh3, 1e3, Uranus_panorama_vts_File);
//    dump_array("jT_nh3", jT_nh3, 1e3, Uranus_panorama_vts_File);

    dump_array("NH4SH", nh4sh, 1e9, Uranus_panorama_vts_File);
    dump_array("w_nh4sh", w_nh4sh, 1e9, Uranus_panorama_vts_File);

//    dump_array("Q_Latent", Q_Latent, 1.0, Uranus_panorama_vts_File);
//    dump_array("Q_Sensible", Q_Sensible, 1.0, Uranus_panorama_vts_File);

    // Radiation pair, matching ATJUP and ATSAT: Q_rad scaled to mW/m3 because the raw W/m3 values
    // are ~1e-5 and plot as a flat zero field, radiation left in W/m2. Both are identically zero
    // unless ATURAN_RADIATION is set.
    dump_array("Q_rad_mW_m3", Q_rad, 1.0e3, Uranus_panorama_vts_File);
    dump_array("Radiation", radiation, 1.0, Uranus_panorama_vts_File);

    pv.close_panorama(Uranus_panorama_vts_File, n);
    return;
}
/*
 * 
*/
void cUranusModel::paraview_vtk_radial(int n, int i_radial){
    using namespace ParaViewIO;
    // File name, header, DIMENSIONS/POINTS and the coordinate block are the shared writer.
    // Display units, settled across all four models: species and their fluxes in g/m3 (1e3), the
    // NH4SH family in ug/m3 (1e9), thermalmassflux in W/m3 (1.0). These are the same units the
    // shared Reporting.h prints, so an array name means one thing in the .vtk and in the log.
    // See the note above the species block in Reporting.h for why NH4SH needs its own.
    ofstream Uranus_vtk_radial_File = ParaViewWriter<cUranusModel>(*this)
        .open_slice("radial", "Radial", i_radial, n, km, jm, 0.1, false);
    const double z = 0.0;   // out-of-plane component of the in-plane vector below
    dump_radial("u-Component", u, u_0, i_radial, Uranus_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Uranus_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Uranus_vtk_radial_File);
    Uranus_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Uranus_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Uranus_vtk_radial_File << t.x[i_radial][j][k] * t_ref - 273.15 << endl;
        }
    }

    dump_radial("thermalmassflux", thermalmassflux, 1.0, i_radial, Uranus_vtk_radial_File);

    dump_radial("H2O", h2o, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("H2OCloud", h2o_cloud, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("H2OIce", h2o_ice, 1e3, i_radial, Uranus_vtk_radial_File);

    dump_radial("CH4", ch4, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("CH4Cloud", ch4_cloud, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("CH4Ice", ch4_ice, 1e3, i_radial, Uranus_vtk_radial_File);

    dump_radial("H2S", h2s, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("H2SCloud", h2s_cloud, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("H2SIce", h2s_ice, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("w_h2s", w_h2s, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("j_h2s", j_h2s, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("jT_h2s", jT_h2s, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("massflux_h2s", massflux_h2s, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("difflux_h2s", difflux_h2s, 1e3, i_radial, Uranus_vtk_radial_File);

    dump_radial("NH3", nh3, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("NH3Cloud", nh3_cloud, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("NH3Ice", nh3_ice, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("w_nh3", w_nh3, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("j_nh3", j_nh3, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("jT_nh3", jT_nh3, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("massflux_nh3", massflux_nh3, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("difflux_nh3", difflux_nh3, 1e3, i_radial, Uranus_vtk_radial_File);

    dump_radial("NH4SH", nh4sh, 1e9, i_radial, Uranus_vtk_radial_File);
    dump_radial("w_nh4sh", w_nh4sh, 1e9, i_radial, Uranus_vtk_radial_File);
    dump_radial("massflux_nh4sh", massflux_nh4sh, 1e9, i_radial, Uranus_vtk_radial_File);
    dump_radial("j_nh4sh", j_nh4sh, 1e9, i_radial, Uranus_vtk_radial_File);
    dump_radial("jT_nh4sh", jT_nh4sh, 1e9, i_radial, Uranus_vtk_radial_File);
    dump_radial("difflux_nh4sh", difflux_nh4sh, 1e9, i_radial, Uranus_vtk_radial_File);


    dump_radial("PressureDyn", p_dyn, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("PressureStat", p_stat, 1.0, i_radial, Uranus_vtk_radial_File);
    dump_radial("rho_mix", rho_mix, 1.0, i_radial, Uranus_vtk_radial_File);

    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Uranus_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Uranus_vtk_radial_File);
    dump_radial("PresGradForce", PresGradForce, 1.0, i_radial, Uranus_vtk_radial_File);

    dump_radial("Q_Latent", Q_Latent, 1.0, i_radial, Uranus_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1.0, i_radial, Uranus_vtk_radial_File);

    // Radiation pair, matching ATJUP and ATSAT. Q_rad in mW/m3, radiation in W/m2;
    // both identically zero unless ATURAN_RADIATION is set.
    dump_radial("Q_rad_mW_m3", Q_rad, 1.0e3, i_radial, Uranus_vtk_radial_File);
    dump_radial("Radiation", radiation, 1.0, i_radial, Uranus_vtk_radial_File);

    Uranus_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Uranus_vtk_radial_File << v.x[i_radial][j][k] << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }
    ParaViewWriter<cUranusModel>(*this).close_slice(Uranus_vtk_radial_File, "radial", i_radial, n);
    return;
}
/*
 * 
*/
void cUranusModel::paraview_vtk_zonal(int n, int k_zonal){
    using namespace ParaViewIO;
    // File name, header, DIMENSIONS/POINTS and the coordinate block are the shared writer.
    // Display units, settled across all four models: species and their fluxes in g/m3 (1e3), the
    // NH4SH family in ug/m3 (1e9), thermalmassflux in W/m3 (1.0). These are the same units the
    // shared Reporting.h prints, so an array name means one thing in the .vtk and in the log.
    // See the note above the species block in Reporting.h for why NH4SH needs its own.
    ofstream Uranus_vtk_zonal_File = ParaViewWriter<cUranusModel>(*this)
        .open_slice("zonal", "Zonal", k_zonal, n, jm, im, 0.05, false);
    const double z = 0.0;   // out-of-plane component of the in-plane vector below
    dump_zonal("u-Component", u, u_0, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Uranus_vtk_zonal_File);
    Uranus_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Uranus_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Uranus_vtk_zonal_File << t.x[i][j][k_zonal] * t_ref - 273.15 << endl;
            aux.x[i][j][k_zonal] = get_layer_height(i);
        }
    }
    dump_zonal("thermalmassflux", thermalmassflux, 1.0, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("height", aux, 1.0, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("H2O", h2o, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("H2OCloud", h2o_cloud, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("H2OIce", h2o_ice, 1e3, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("CH4", ch4, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("CH4Cloud", ch4_cloud, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("CH4Ice", ch4_ice, 1e3, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("H2S", h2s, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("H2SCloud", h2s_cloud, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("H2SIce", h2s_ice, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("w_h2s", w_h2s, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("j_h2s", j_h2s, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("jT_h2s", jT_h2s, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("massflux_h2s", massflux_h2s, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("difflux_h2s", difflux_h2s, 1e3, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("NH3", nh3, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("NH3Cloud", nh3_cloud, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("NH3Ice", nh3_ice, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("w_nh3", w_nh3, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("j_nh3", j_nh3, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("jT_nh3", jT_nh3, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("massflux_nh3", massflux_nh3, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("difflux_nh3", difflux_nh3, 1e3, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("NH4SH", nh4sh, 1e9, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("w_nh4sh", w_nh4sh, 1e9, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("massflux_nh4sh", massflux_nh4sh, 1e9, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("j_nh4sh", j_nh4sh, 1e9, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("jT_nh4sh", jT_nh4sh, 1e9, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("difflux_nh4sh", difflux_nh4sh, 1e9, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("PressureDyn", p_dyn, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("PressureStat", p_stat, 1.0, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("rho_mix", rho_mix, 1.0, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("CentrifugalForce", CentrifugalForce, 1e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("PresGradForce", PresGradForce, 1.0, k_zonal, Uranus_vtk_zonal_File);

    dump_zonal("Q_Latent", Q_Latent, 1.0, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1.0, k_zonal, Uranus_vtk_zonal_File);

    // Radiation pair, matching ATJUP and ATSAT. Q_rad in mW/m3, radiation in W/m2;
    // both identically zero unless ATURAN_RADIATION is set.
    dump_zonal("Q_rad_mW_m3", Q_rad, 1.0e3, k_zonal, Uranus_vtk_zonal_File);
    dump_zonal("Radiation", radiation, 1.0, k_zonal, Uranus_vtk_zonal_File);

    Uranus_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Uranus_vtk_zonal_File << u.x[i][j][k_zonal] << " " << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }
    ParaViewWriter<cUranusModel>(*this).close_slice(Uranus_vtk_zonal_File, "zonal", k_zonal, n);
    return;
}
/*
 * 
*/
void cUranusModel::paraview_vtk_longal(int n, int j_longal){
    using namespace ParaViewIO;
    // File name, header, DIMENSIONS/POINTS and the coordinate block are the shared writer.
    // Display units, settled across all four models: species and their fluxes in g/m3 (1e3), the
    // NH4SH family in ug/m3 (1e9), thermalmassflux in W/m3 (1.0). These are the same units the
    // shared Reporting.h prints, so an array name means one thing in the .vtk and in the log.
    // See the note above the species block in Reporting.h for why NH4SH needs its own.
    ofstream Uranus_vtk_longal_File = ParaViewWriter<cUranusModel>(*this)
        .open_slice("longal", "Longitudinal", j_longal, n, km, im, 0.025, true);
    const double y = 0.0;   // out-of-plane component of the in-plane vector below
    dump_longal("u-Component", u, u_0, j_longal, Uranus_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Uranus_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Uranus_vtk_longal_File);
    Uranus_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Uranus_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Uranus_vtk_longal_File << t.x[i][j_longal][k] * t_ref - 273.15 << endl;
            aux.x[i][j_longal][k] = get_layer_height(i);
        }
    }

    dump_longal("thermalmassflux", thermalmassflux, 1.0, j_longal, Uranus_vtk_longal_File);

    dump_longal("height", aux, 1.0, j_longal, Uranus_vtk_longal_File);

    dump_longal("H2O", h2o, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("H2OCloud", h2o_cloud, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("H2OIce", h2o_ice, 1e3, j_longal, Uranus_vtk_longal_File);

    dump_longal("CH4", ch4, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("CH4Cloud", ch4_cloud, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("CH4Ice", ch4_ice, 1e3, j_longal, Uranus_vtk_longal_File);

    dump_longal("H2S", h2s, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("H2SCloud", h2s_cloud, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("H2SIce", h2s_ice, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("w_h2s", w_h2s, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("j_h2s", j_h2s, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("jT_h2s", jT_h2s, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("massflux_h2s", massflux_h2s, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("difflux_h2s", difflux_h2s, 1e3, j_longal, Uranus_vtk_longal_File);

    dump_longal("NH3", nh3, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("NH3Cloud", nh3_cloud, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("NH3Ice", nh3_ice, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("w_nh3", w_nh3, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("j_nh3", j_nh3, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("jT_nh3", jT_nh3, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("massflux_nh3", massflux_nh3, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("difflux_nh3", difflux_nh3, 1e3, j_longal, Uranus_vtk_longal_File);

    dump_longal("NH4SH", nh4sh, 1e9, j_longal, Uranus_vtk_longal_File);
    dump_longal("w_nh4sh", w_nh4sh, 1e9, j_longal, Uranus_vtk_longal_File);
    dump_longal("massflux_nh4sh", massflux_nh4sh, 1e9, j_longal, Uranus_vtk_longal_File);
    dump_longal("j_nh4sh", j_nh4sh, 1e9, j_longal, Uranus_vtk_longal_File);
    dump_longal("jT_nh4sh", jT_nh4sh, 1e9, j_longal, Uranus_vtk_longal_File);
    dump_longal("difflux_nh4sh", difflux_nh4sh, 1e9, j_longal, Uranus_vtk_longal_File);


    dump_longal("PressureDyn", p_dyn, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("PressureStat", p_stat, 1.0, j_longal, Uranus_vtk_longal_File);
    dump_longal("rho_mix", rho_mix, 1.0, j_longal, Uranus_vtk_longal_File);

    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Uranus_vtk_longal_File);
    dump_longal("CentrifugalForce", CentrifugalForce, 1e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1.0, j_longal, Uranus_vtk_longal_File);
    dump_longal("PresGradForce", PresGradForce, 1.0, j_longal, Uranus_vtk_longal_File);

    dump_longal("Q_Latent", Q_Latent, 1.0, j_longal, Uranus_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1.0, j_longal, Uranus_vtk_longal_File);

    // Radiation pair, matching ATJUP and ATSAT. Q_rad in mW/m3, radiation in W/m2;
    // both identically zero unless ATURAN_RADIATION is set.
    dump_longal("Q_rad_mW_m3", Q_rad, 1.0e3, j_longal, Uranus_vtk_longal_File);
    dump_longal("Radiation", radiation, 1.0, j_longal, Uranus_vtk_longal_File);

    Uranus_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Uranus_vtk_longal_File << u.x[i][j_longal][k] << " " 
                << y << " " << w.x[i][j_longal][k] << endl;
        }
    }
    ParaViewWriter<cUranusModel>(*this).close_slice(Uranus_vtk_longal_File, "longal", j_longal, n);
    return;
}
/*
 * 
*/
void cUranusModel::paraview_sphere_vts(int n){
    using namespace ParaViewIO;
    double x, y, z, sinthe, sinphi, costhe, cosphi;
    // Display units, settled across all four models: species and their fluxes in g/m3 (1e3), the
    // NH4SH family in ug/m3 (1e9), thermalmassflux in W/m3 (1.0). These are the same units the
    // shared Reporting.h prints, so an array name means one thing in the .vtk and in the log.
    // See the note above the species block in Reporting.h for why NH4SH needs its own.
    string Uranus_sphere_vts_File_Name = output_path + "/Uranus_sphere_" 
        + std::to_string(n) + ".vts";
    ofstream Uranus_sphere_vts_File;
    Uranus_sphere_vts_File.precision(4);
    Uranus_sphere_vts_File.setf(ios::fixed);
    Uranus_sphere_vts_File.open(Uranus_sphere_vts_File_Name);
    if (!Uranus_sphere_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Uranus_sphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Uranus_sphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Uranus_sphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Uranus_sphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Uranus_sphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDyn PressureStat  H2S H2SCloud H2SIce NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce \">\n"  << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        sinphi = sin( phi.z[k]);
        cosphi = cos( phi.z[k]);
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            for(int i = 0; i < im; i++){
                aux_u.x[i][j][k] = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi * v.x[i][j][k] - sinphi * w.x[i][j][k];
                aux_v.x[i][j][k] = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe * v.x[i][j][k] + cosphi * w.x[i][j][k];
                aux_w.x[i][j][k] = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];
                Uranus_sphere_vts_File << aux_u.x[i][j][k] << " " << aux_v.x[i][j][k] << " " << aux_w.x[i][j][k]  << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << t.x[i][j][k] * t_ref - 273.15 << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDyn\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * p_dyn.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureStat\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e-3 * p_stat.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2O\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2o.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2S\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2s.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * nh3.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH4SH\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e9 * nh4sh.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OCloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2o_cloud.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2SCloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2s_cloud.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Cloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * nh3_cloud.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Uranus_sphere_vts_File <<  "\n"  << endl;

    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OIce\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2o_ice.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2SIce\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * h2s_ice.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Ice\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * nh3_ice.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
/*
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"HE\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << 1e3 * he.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
*/
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << u_0 * aux_u.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << u_0 * aux_v.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << u_0 * aux_w.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Seamount\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Uranus_sphere_vts_File << SeaMount.x[i][j][k] << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Uranus_sphere_vts_File <<  "\n"  << endl;
    Uranus_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Uranus_sphere_vts_File <<  "   </PointData>\n" << endl;
    Uranus_sphere_vts_File <<  "   <Points>\n"  << endl;
    Uranus_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                x = rad.z[i] * sin( the.z[j])* cos(phi.z[k]);
                y = rad.z[i] * sin( the.z[j])* sin(phi.z[k]);
                z = rad.z[i] * cos( the.z[j]);
                Uranus_sphere_vts_File << x << " " << y << " " << z  << endl;
            }
            Uranus_sphere_vts_File <<  "\n"  << endl;
        }
        Uranus_sphere_vts_File <<  "\n"  << endl;
    }
    Uranus_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Uranus_sphere_vts_File <<  "   </Points>\n"  << endl;
    Uranus_sphere_vts_File <<  "  </Piece>\n"  << endl;
    Uranus_sphere_vts_File <<  " </StructuredGrid>\n"  << endl;
    Uranus_sphere_vts_File <<  "</VTKFile>\n"  << endl;
    Uranus_sphere_vts_File.close();
    cout << "   File:  " << "Uran_sphere_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
}
/*
 * 
*/
void cUranusModel::UranusPlotData(){
    string Name_PlotData_File = output_path + "/PlotData_Uranus.xyz";
    ofstream PlotData_File;
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    PlotData_File.open(Name_PlotData_File);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File << "lons(deg)" << ", " << "lats(deg)" << ", " 
        << "topography" << ", " << "v-velocity(m/s)" << ", " 
        << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", " 
        << "temperature(Celsius)" << ", " << "water_vapour(g/kg)" 
        << ", " << "precipitation(mm)" << ", " 
        <<  "precipitable water(mm)" << endl;
    double vel_mag;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            vel_mag = sqrt(pow(v.x[0][j][k] * u_0, 2) + pow(w.x[0][j][k] * u_0, 2));
//            PlotData_File << k << " " << j << " " << SeaMount.x[0][j][k] << " " 
            PlotData_File << k << " " << j << " "
                << v.x[0][j][k] * u_0 << " " << w.x[0][j][k] * u_0 << " " 
                << vel_mag << " " << t.x[0][j][k] * t_ref - t_ref << " " 
                << h2o.x[0][j][k] << " "<< nh3.x[0][j][k] <<  endl;
        }
    }
    PlotData_File.close();
    return;
}



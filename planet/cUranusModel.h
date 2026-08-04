#ifndef _CJUPITERMODEL_H
#define _CJUPITERMODEL_H

#include <fenv.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"
#include "PythonStream.h"
#include "Utils.h"
#include "Config.h"


#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

namespace{
    std::function<double(double)> default_lambda=[](double i)->double{return i;};
}

class cUranusModel{

    // Shared physics/output templates (SHARED.md5, `make check-shared`). ATURAN is the fourth
    // model to take these; see ParaViewWriter.h for what it provides and what it does not.
    template<class M> friend class ParaViewWriter;
    template<class M> friend class ConvectiveAdjustment;
    template<class M> friend class FluxLimiter;
    friend class ChemistryUran;
    friend class PressureSolverUran;
    friend class SaturationAdjustmentUran;
    friend class BC_Uran;
    friend class VelocityInitializerUran;

public:

    // ---- Hooks for the shared ParaViewWriter<Planet> (ParaViewWriter.h) ----
    // planet_name() is the word in an output FILE name ("Uranus_radial_20_1.vtk");
    // planet_short() is the abbreviation inside a .vtk title line
    // ("Radial_Data_Uran_Circulation"). ATURAN carried both spellings by hand, and like ATNEPT
    // its "has been written" line used the SHORT one — announcing Uran_radial_20_1.vtk for a
    // file actually called Uranus_radial_20_1.vtk. The shared writer names the file it wrote.
    static const char* planet_name(){ return "Uranus"; }
    static const char* planet_short(){ return "Uran"; }

    // The model's own name in log lines written by SHARED code — "ATURAN: ..." — so a shared
    // header can say which planet it is running on without knowing anything else about it. It is
    // also the prefix the shared modules build their environment-variable names from. ATSAT,
    // ATJUP and ATNEPT carry the same accessor.
    static const char* planet_tag(){ return "ATURAN"; }

    // ---- Hooks for the shared FluxLimiter<Planet> (FluxLimiter.h) ----
    //
    // metricRadius() is the established hook for the one place the models genuinely differ in the
    // limiter: ATJUP shifts rad.z itself at initialisation and so returns rm unchanged, while
    // ATSAT shifts the metric factors here instead. ATURAN is in ATJUP's position for the simpler
    // reason that it has no metric radius at all, so this is the identity unless
    // ATURAN_METRIC_RADIUS is set, and the shared limiter reproduces the m.rad.z[i] the
    // hand-written copy used, exactly.
    double metricRadius(double rm){
        static const double R_km = [](){
            const char* e = getenv("ATURAN_METRIC_RADIUS"); return e ? atof(e) : 0.0; }();
        if(!(R_km > 0.0)) return rm;
        return rm + (R_km / L_atm - 1.0);
    }

    // The model's own floor on sin(theta) in the METRIC. ATURAN declares none, as ATSAT and
    // ATNEPT do not; it exists so ATPhys::polar_divisor_floor<Planet>() compiles. With
    // ATURAN_SINTHE_TRACK unset that function returns the literal 0.4 the hand-written limiter
    // used, so this value is not reached by default.
    static double sinthe_min(){
        static const double v = [](){
            const char* e = getenv("ATURAN_SINTHE_MIN");
            const double x = e ? atof(e) : 0.0;
            return (x >= 0.0 && x < 1.0) ? x : 0.0;
        }();
        return v;
    }

    // What the panorama .vts prints in its "Temperature" array — degrees Celsius, as ATJUP and
    // ATNEPT do; ATSAT writes kelvin/10, and one array name carrying two quantities across four
    // models remains unsettled.
    //
    // ATURAN DID NOT COMPUTE THIS AT ALL. Four of its five ParaView temperatures read
    //     t.x[i][j][k] * 273.15 - 273.15
    // with 273.15 standing where t_ref belongs. t_ref is 76.4 on Uranus, so the constant was
    // wrong by a factor of 273.15/76.4 = 3.5752, and those temperatures were meaningless — the
    // initial deep equator, 421.0 K = t_nd 5.5105, printed as 1232.0 degC where it should read
    // 147.9. The telling part is that the ZONAL writer already used t_ref correctly, so the model
    // disagreed with ITSELF: one slice was in degrees Celsius and the other three plus the
    // panorama were in nothing at all. Fixed here and at the three remaining call sites.
    double paraview_temperature(double t_nd) const { return t_nd * t_ref - 273.15; }

    const char *filename;

    cUranusModel();
    ~cUranusModel();

    cUranusModel(const cUranusModel&) = delete;
    cUranusModel& operator=(const cUranusModel&) = delete;

    static cUranusModel* get_model(){
        if(!m_model){
            m_model = new cUranusModel();
        }
        return m_model;
    }

    void LoadConfig(const char *filename);
    void Run();

    #include "UranusParams.h.inc"

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0;


    int n, n_print, n_paraview, panorama, panorama_step;
    int iter, iter_max, j_max;
    int velocity_iter, pressure_iter;
    int iter_n, panorama_cnt;
    int i_res, j_res, k_res;

    std::vector<double> tropopause_layers; // keep the tropopause layer index
    std::vector<std::vector<int> > i_topography;
    std::vector<double> u_trans;
    std::vector<double> v_trans;
    std::vector<double> w_trans;

    double maxValue, minValue;
    /*
     * Given a latitude, return the layer index of tropopause
    */
    int get_tropopause_layer(int j){
        assert(j>=0);
        assert(j<jm);
        //refer to  BC_Thermo::TropopauseLocation and BC_Thermo::GetTropopauseHightAdd
        //tropopause height is proportional to the mean tropospheric temperature.
        //higher near the equator - warm troposphere
        //lower at the poles - cold troposphere

//        tropopause_layers[j] = 30;
//        tropopause_layers[j] = 35;
//        init_tropopause_layers();
        tropopause_layers[j] = im_tropopause[j];
        return tropopause_layers[j];
    }
    /*
     *
    */
    int get_surface_layer(int j, int k){
        return i_topography[j][k];
    }
    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    std::vector<float> get_layer_heights(){
        return m_layer_heights;
    }   
    /*
    * Given a altitude, return the layer index
    */
    int get_layer_index(float height){
        std::size_t i = 0;
        for(; i<m_layer_heights.size(); i++){
            if(height<m_layer_heights[i])
                return i-1;
        }
        return i;
    }



private:

    static cUranusModel* m_model;

    PythonStream ps;
    std::streambuf *backup;

    const double c43 = 4.0/3.0, c13 = 1.0/3.0, c32 = 3.0/2.0, c42 = 4.0/2.0, c12 = 1.0/2.0;
    static const int im = 41, jm = 181, km = 361;

    int i_max = 40;  // corresponds to about 125 km above 10e6 Pa pressure level, maximum hight of the tropopause at equator
    int i_beg = 30;  // corresponds to about 100 km above 10e6 Pa pressure level, maximum hight of the tropopause at poles
//    int i_beg = 40;  // corresponds to about 100 km above 10e6 Pa pressure level, maximum hight of the tropopause at poles

    double gam = 0.0;    // dry adiabatic lapse rate [K/km], computed as g*1e3/cp_mix after ThermalPropertiesUran()
    double mue_mix, k_mix, cp_mix, rho_cond_mix, r_mix, R_mix, c_mix, x_mix;

// at 230K NH3 and H2S condense via a heterogenious reaction: NH3 + H2S -> NH4SH ( Planetary Scienses, de Pater, Lissauer)


// temperatures at triple point and ice formation
    double t_0_h2o = 273.15;  // in K == 0°C, triple point
    double t_00_h2o = 210.15;  // in K == -67°C (Planetary Siences)
//    double t_00_h2o = 241.15;  // in K == -32°C (COSMO)
    double t_000 = 235.15;  // in K == -20°C (precipitation module)

    double t_0_h2s = 187.65;  // in K == -85.5°C, h2s-ice cloud formation (Planetary Sciences, p. 96)
    double t_00_h2s = 220.0;  // in K == -133.15°C, h2s-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_ch4 = 90.69;  // in K == -182.456°C, triple point
    double t_00_ch4 = 190.56;  // in K == -85.5°C, h2s-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh3 = 195.5;  // in K == -77.65°C, triple point, gas and liquid pase
    double t_00_nh3 = 220.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)
//    double t_00_nh3 = 140.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)
//    double t_00_nh3 = 120.0;  // in K == -153.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh4sh = 230.0;  // in K == -43.15°C, nh4sh formation onset (Planetary Sciences, p. 96)
    double t_00_nh4sh = 200.0;  // in K == -73.15°C, nh4sh formation end (Planetary Sciences, p. 96)


// pressures take from  Planetary Sciences, p. 96
    double p_0_h2o = 21.0;  // in bar
    double p_00_h2o = 50.0;  // in bar

    double p_0_ch4 = 0.1;  // in bar
    double p_00_ch4 = 1.1;  // in bar

    double p_0_h2s = 2.8;  // in bar
    double p_00_h2s = 7.0;  // in bar

    double p_0_nh3 = 20.2;  // in bar
    double p_00_nh3 = 38.0;  // in bar

    double p_0_nh4sh = 20.5;  // in bar
    double p_00_nh4sh = 38.0;  // in bar


// constants for Clausius-Clapeyron law
    double coeff_h2_A = -2000.0; // invented
    double coeff_h2_B = 8.0; // invented

    double coeff_he_A = -3000.0; // invented
    double coeff_he_B = 10.0; // invented


    double coeff_h2o_A = -4961.04;  // from triple and critical point values for water
    double coeff_h2o_B = 13.0662;  // from triple and critical point values for water

    double coeff_h2o_A_i = -4961.04; // invented
    double coeff_h2o_B_i = 13.0662; // invented


    double coeff_ch4_A = -1033.3;  // from triple and critical point values for methane
    double coeff_ch4_B = 6.3910;  // from triple and critical point values for methane

    double coeff_ch4_A_i = -1033.3; // invented
    double coeff_ch4_B_i = 6.3910; // invented

    double coeff_h2s_A = -2251.66;  // from triple and critical point values for hydrogen sulfide
    double coeff_h2s_B = 10.5253;  // from triple and critical point values for hydrogen sulfide

    double coeff_h2s_A_i = -2251.66; // invented
    double coeff_h2s_B_i = 10.5253; // invented


    double coeff_nh3_A = -2836.56;  // from triple and critical point values for ammonia
    double coeff_nh3_B = 11.7271;  // from triple and critical point values for ammonia

    double coeff_nh3_A_i = -2836.56; // invented
    double coeff_nh3_B_i = 11.7271; // invented


    double coeff_nh4sh_A = -2836.56;  // invented
    double coeff_nh4sh_B = 11.7271; // invented


// condensate/crystal densities in kg/m³ (Planetary Sciences p. 90, 2010).
// These are phase-change densities for cloud microphysics — NOT for gas-phase diffusion.
    double rho_cond_h2    = 0.408;    // liquid hydrogen density in kg/m³
    double rho_cond_he    = 0.1752;   // liquid helium density in kg/m³
    double rho_cond_nh3   = 0.7623;   // liquid ammonia density in kg/m³
    double rho_cond_nh4sh = 1170.0;   // solid ammonium hydrosulfide density in kg/m³
    double rho_cond_h2s   = 1.5357;   // liquid hydrogen sulfide density in kg/m³
    double rho_cond_h2o   = 1000.0;    // liquid water density in kg/m³
    double rho_cond_ch4   = 0.657;    // liquid methane density in kg/m³

// molecular weights
    double m_h2 = 2.016;  // molecular weight of hydrogen in kg/Kmol (molar mass)
    double m_he = 4.02602;  // molecular weight of helium in kg/Kmol
    double m_nh3 = 17.03052;  // molecular weight of ammonia in kg/Kmol
    double m_nh4sh = 51.1114;  // molecular weight of ammonium hydrosulfide in kg/Kmol
    double m_h2s = 34.08088;  // molecular weight of hydrogen sulfide in kg/Kmol
    double m_h2o = 18.01588;  // molecular weight of water in kg/Kmol
    double m_ch4 = 16.042;  // molecular weight of water in kg/Kmol

// vapour mass densities of gases                Planetary sciences p. 90 2010
    double r_h2 = 0.864;  // density of hydrogen vapour in kg/m³
    double r_he = 0.136;  // density of helium vapour  in kg/m³
    double r_ch4 = 0.19;  // density of methane vapour in kg/m³
    double r_nh3 = 0.1;  // density of ammonia vapour in kg/m³
    double r_h2s = 0.3;  // density of hydrogen sulfide vapour in kg/m³
    double r_h2o = 0.08;  // density of water vapour in kg/m³
    double r_nh4sh = 0.007;  // density of ammonium hydrosulfide vapour in kg/m³  assumption
    double r_nh3_add = 0.09;  // density of ammonia vapour in kg/m³

// vapour mass densities of clouds and ices               Planetary sciences p. 90 2010
    double r_nh3_cloud = 0.004;  // density of ammonia cloud in kg/m³
    double r_h2o_cloud = 0.009;  // density of water cloud in kg/m³
    double r_nh3_ice = 0.0003;  // density of ammonia ice in kg/m³
    double r_h2o_ice = 0.0038;  // density of water ice in kg/m³
    double r_h2s_ice = 0.007;  // density of water ice in kg/m³
    double r_ch4_ice = 0.0082;  // density of water ice in kg/m³

// NH4SH particle sedimentation
    double r_p_nh4sh = 1.0e-5;  // NH4SH crystal radius in m

// vapour molar densities of gases
    double c_h2 = r_h2/m_h2;  // density of hydrogen vapour in kmol/m³
    double c_he = r_he/m_he;  // density of helium vapour  in kmol/m³
    double c_nh3 = r_nh3/m_nh3;  // density of ammonia vapour in kmol/m³
    double c_h2s = r_h2s/m_h2s;  // density of hydrogen sulfide vapour in kmol/m³
    double c_h2o = r_h2o/m_h2o;  // density of water vapour in kmol/m³
    double c_ch4 = r_ch4/m_ch4;  // density of water vapour in kmol/m³
    double c_nh4sh = r_nh4sh/m_nh4sh;  // density of ammonium hydrosulfide vapour in kmol/m³  assumption

 // ratio of vapour molecular weight to mean molecular weight
    double X_h2 = 0.864;
    double X_he = 0.136;
    double X_h2o = 1.7e-3;
    double X_nh3 = 2.0e-4;
    double X_h2s = 7.7e-5;
    double X_nh4sh = 3.6e-5;
    double X_ch4 = 3.6e-5;

// gas constants
    double R_h2 = 4124.2; // gas constant of hydrogen in J/(kg*K)
    double R_he = 2077.1; // gas constant of helium in J/(kg*K)
    double R_nh3 = 488.21; // gas constant of ammoinia in J/(kg*K)
    double R_nh4sh = 261.0; // gas constant of ammoinium hydrosulfide in J/(kg*K)   invented for the initial distribution of nh4sh
    double R_h2s = 243.96; // gas constant of hydrogen sulfide inJ/(kg*K) 
    double R_h2o = 461.52; // gas constant of water inJ/(kg*K)
    double R_ch4 = 518.28; // gas constant of methane inJ/(kg*K)

// dynamic viscosities
    double mue_h2 = 0.84e-5; // dynamic viscosity of hydrogen in Ns/m²
    double mue_he = 1.87e-5; // dynamic viscosity of helium in Ns/m²
    double mue_nh3 = 0.92e-5; // dynamic viscosity of ammonia in Ns/m²
    double mue_nh4sh = 0.99e-5; // dynamic viscosity of ammonium sulfide in Ns/m²
    double mue_h2s = 1.3e-5; // dynamic viscosity of hydrogen sulfid in Ns/m²
    double mue_h2o = 1.308e-3; // dynamic viscosity of water in Ns/m²
    double mue_ch4 = 1.107e-2; // dynamic viscosity of methane in Ns/m²

// thermal conductivities
    double k_h2 = 0.1317; // thermal conductivity of hydrogen in W/(m*K)
    double k_he = 0.1193; // thermal conductivity of helium in W/(m*K)
    double k_nh3 = 0.02102; // thermal conductivity of ammonia in W/(m*K)
    double k_nh4sh = 0.02102; // thermal conductivity of ammonium hydrosulfide in W/(m*K)
    double k_h2s = 0.013; // thermal conductivity of hydrogen sulfide in W/(m*K)
    double k_h2o = 0.0187; // thermal conductivity of water in W/(m*K)
    double k_ch4 = 0.0; // thermal conductivity of methane in W/(m*K)

// specific heat capacities
    double cp_h2 = 14.32e3;  // specific heat capacity of hydrogen in J/(kg*K)
    double cp_he = 5.19e3;  // specific heat capacity of helium in J/(kg*K)
    double cp_nh3 = 2.19e3;  // specific heat capacity of ammonia in J/(kg*K)
    double cp_nh4sh = 2.00e3;  // specific heat capacity of ammonium hydrosulfide in J/(kg*K)
    double cp_h2s = 2.24e3;  // specific heat capacity of hydrogen sulfid in J/(kg*K)
    double cp_h2o = 1.93e3;  // specific heat capacity of water in J/(kg*K)
    double cp_ch4 = 2.232e3;  // specific heat capacity of methane in J/(kg*K)

// ratios of gas constants of dry gas to vapour or vapour molecular weight to mean atmospheric molecular weight or m/m_mix
    double ep_h2 = 0.8572;  // ratio of the gas constants of dry air to h2 non-dimensional or m/m_mix
    double ep_he = 1.7152;  // ratio of the gas constants of dry he to h2 non-dimensional
    double ep_h2o = 8.1253;  // ratio of the gas constants of dry air to h2 non-dimensional
    double ep_h2s = 14.5192;  // ratio of the gas constants of dry hydrogen sulfide to h2 non-dimensional
    double ep_nh3 = 7.6752;  // ratio of the gas constants of dry ammonia to h2 non-dimensional
    double ep_ch4 = 7.6752;  // ratio of the gas constants of dry methane to h2 non-dimensional
    double ep_nh4sh = 21.7745;  // ratio of the gas constants of dry methane to h2 non-dimensional       invented

// latent heat of evaporation
    double lv_h2o = 2.5009e6;  // latent heat of h2o evaporation at 0°C in J/kg
    double lv_h2s = 3.5340e6;  // latent heat of h2s evaporation at -73°C in J/kg
    double lv_nh3 = 1.3720e6;  // latent heat of nh3 evaporation at -33.33 in J/kg
    double lv_ch4 = 5.11e5;  // latent heat of nh3 evaporation at -33.33 in J/kg

// latent heat of sublimation
    double ls_h2o = 2.8339e6;  // latent heat of h2o sublimation at 0°C in J/kg
    double ls_h2s = 7.4500e5;  // latent heat of h2s sublimation at -98°C in J/kg
    double ls_nh3 = 1.8320e6;  // latent heat of nh3 sublimation at -93.15 in J/kg
    double ls_ch4 = 5.11e5;  // latent heat of ch4 sublimation at -93.15 in J/kg

// Schmidt number
    double sc_h2 = 0.20;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_he = 0.22;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_ch4 = 0.99;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2o = 0.61;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2s = 0.94;  // Schmidt number of h2s, Sc = nue/D 
    double sc_nh3 = 0.61;  // Schmidt number of nh3, Sc = nue/D 
    double sc_nh4sh = 0.7;  // Schmidt number of nh4sh, Sc = nue/D

// Sponge layer (top quarter of domain, quadratic Rayleigh damping on u)
    double alpha_sponge = 10.0;  // peak damping rate [1/time_nd] at i=im-1

// Prantl numbers
    double Pr = 0.72;  // Prandtl-number 

// diffusion coefficients
    double D_nh3 = 1.5e-9; // ordinary diffusion coefficient of ammonia in m*m/s 
    double D_h2s = 1.36e-9; // ordinary diffusion coefficient of hydrogen sulfid in m*m/s 
    double D_nh4sh = 1.45e-9; // ordinary diffusion coefficient of ammonium hydrosulfide in m*m/s

// thermal diffusion coefficients
    double DT_nh3 = 1.54e-9; // thermal diffusion coefficient of ammonia in kg/(s*m)                           unklar
    double DT_h2s = 1.36e-9; // thermal diffusion coefficient of hydrogen sulfid in kg/(s*m)
    double DT_nh4sh = 1.45e-9; // thermal diffusion coefficient of ammonium hydrosulfide in kg/(s*m)

// constants for saturation vapour pressure and latent heat
    double C_h2o = 25.096;  //  in bar
    double C_ch4 = 1.627;  //  in bar
    double C_nh3 = 27.863;  //  in bar
    double C_h2s = 17.064;  //  in bar
    double C_nh4sh = 75.678;  //  in bar

    double L0_h2o = 3148.2;  //  in J/g
    double L0_ch4 = 553.1;  //  in J/g
    double L0_nh3 = 2016.0;  //  in J/g
    double L0_h2s = 747.0;  //  in J/g
    double L0_nh4sh = 2915.7;  //  in J/g

    double del_alf_h2o = 0.0;
    double del_bet_h2o = - 8.7e-3;

    double del_alf_ch4 = 1.002;
    double del_bet_ch4 = -4.1e-3;

    double del_alf_nh3 = - 0.888;
    double del_bet_nh3 = 0.0;

    double del_alf_h2s = 0.0;
    double del_bet_h2s = - 2.9e-3;

    double del_alf_nh4sh = - 1.760;
    double del_bet_nh4sh = 7.8e-4;
 
    std::vector<std::vector<int> > j_ellipse;
    bool has_welcome_msg_printed;
    double out_maxValue() const;
    double out_minValue() const;

    void init_layer_heights(){
        float h = L_atm/(im-1);
        for(int i=0; i<im; i++){
            m_layer_heights.push_back(i * h);
        } 
        return;
    }

    struct CellGeometry {
        double rm, rm2;
        double sinthe, sinthe2, costhe, cotanthe;
        double inv_rm, inv_rm2;
        double inv_rmsinthe, inv_rm2sinthe, inv_rm2sinthe2;
        double costhe_inv_rm2sinthe;
        double inv_2dr, inv_2dthe, inv_2dphi;
        double inv_dr2, inv_dthe2, inv_dphi2;
    };

    void SetDefaultConfig();
    void RHSUran(int i, int j, int k, const CellGeometry& geo);
    void RungeKuttaUran();
    void UranusPlotData();
    void paraview_vtk_longal(int n, int j_longal);
    void paraview_vtk_radial(int n, int i_radial);
    void paraview_vtk_zonal(int n, int k_zonal);
    void paraview_panorama_vts(int n);
    void paraview_sphere_vts(int n);

    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1., 
        std::function< double(double) > lambda = default_lambda,
        bool print_heading=false);

    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.0);

    void print_welcome_msg();
    void print_final_msg();
    void printMinMax();
    void initMsg();
    void writeResults();
    void writeData();

    void BC_phi();
    void BC_radius();
    void BC_theta();
    void BC_seamount();
    void BC_solidground();
    void resetArrays();
    void TropopauseLocation();
    void UranusCellStructure();
    void UranusCellStructure_new();
    void Uraniter_PlotData();

    void init_tropopause_layers();
    void init_u(Array &u, int j);
    void form_diagonals(Array &a, int start, int end);
    void init_v_or_w(Array &v_or_w, int j, double coeff_trop, double coeff_sl);
    void init_v_or_w_above_tropopause(Array &v_or_w, int j, double coeff);

    void computePressure();
    void init_temperature();
    void init_PressureStatic();
    void init_PressureDynamic();

    void init_vapour(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &X,
        Array &c, Array &cloud, Array &ice, Array &cloudiness);

    void init_vapour_cloud_ice(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet,
        Array &c, Array &cloud, Array &ice, Array &cloudiness);

    void init_h2s(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, Array &c);

    void Saturation_Adjustment(std::string gas, 
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00, 
        double &ep, double &lv, double &ls, double &cp, double &r,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &m,
        Array &c, Array &cloud, Array &ice);

    void OneCategoryIceScheme();

    void ChemMassRateUran();
    void DiffMassFluxUran();
    void ThermalPropertiesUran();
    void Latent_Heat();
    void Forces();

    void steadyQuery();
    void restoreVar(double coeff);

    std::vector<int> im_tropopause; // keep the tropopause layer index
    std::vector<float> m_layer_heights;
    std::vector<double> cloud_loc; // lateral cloudwater distribution
    std::vector<double> r_max; // lateral r_max distribution
    std::vector<double> r_max_add; // lateral r_max distribution
    std::vector<double> t_add; // lateral r_max distribution

    Array_1D rad;
    Array_1D the;
    Array_1D phi;

    Array_2D Topography; // topography
    Array_2D LatentHeat;        // areas of higher latent heat
    Array_2D Precipitation;        // areas of higher precipitation
    Array_2D precipitable_water;// areas of precipitable water in the air
    Array_2D nh3_total;            // areas of higher nh3 concentration
    Array_2D nh3_cloud_total;    // areas of higher nh3_cloud concentration
    Array_2D nh3_ice_total;        // areas of higher nh3_ice concentration
    Array_2D aux_2D_v;            // auxilliar field v
    Array_2D aux_2D_w;            // auxilliar field w
    Array_2D tropopause_height; // local height of the tropopause

    Array t;                    // temperature
    Array u;                    // u-component velocity component in r-direction
    Array v;                    // v-component velocity component in theta-direction
    Array w;                    // w-component velocity component in phi-direction

    Array ch4;                    // methane vapour
    Array ch4_cloud;                // methane cloud
    Array ch4_ice;                    // methane ice
    Array h2o;                    // water vapour
    Array h2o_cloud;                // cloud water
    Array h2o_ice;                    // cloud ice
    Array h2s;                    // water vapour
    Array h2s_cloud;                // cloud water
    Array h2s_ice;                    // cloud ice
    Array nh3;                    // nh3-vapour
    Array nh3_cloud;            // nh3-cloud
    Array nh3_ice;                // nh3-ice
    Array nh4sh;                    // nh4sh-vapour

    Array tn;                    // temperature new
    Array un;                    // u-velocity component in r-direction new
    Array vn;                    // v-velocity component in theta-direction new
    Array wn;                    // w-velocity component in phi-direction new
    Array ch4n;                    // water vapour new
    Array ch4_cloudn;                    // water vapour new
    Array ch4_icen;                // cloud water new
    Array h2on;                    // water vapour new
    Array h2o_cloudn;                    // water vapour new
    Array h2o_icen;                // cloud water new
    Array h2sn;                    // water vapour new
    Array h2s_cloudn;                    // water vapour new
    Array h2s_icen;                // cloud water new
    Array nh3n;                    // nh3 new
    Array nh3_cloudn;            // nh3_cloud new
    Array nh3_icen;                 // nh3_ice new
    Array nh4shn;                    // nh4sh new

    Array massflux_h2s;    // mass flux h2s
    Array massflux_nh3;    // mass flux nh3
    Array massflux_nh4sh;  // mass flux nh4sh
    Array fluxlim_nh4sh;   // TVD flux-limiter correction for nh4sh advection

    Array difflux_h2s;   // diffusive flux h2s
    Array difflux_nh3;   // diffusive flux nh3
    Array difflux_nh4sh; // diffusive flux nh4sh

    Array thermalmassflux;   // thermal massflux_h2s

    Array cloudiness_ch4; // cloudiness, N in literature
    Array cloudiness_h2o; // cloudiness, N in literature
    Array cloudiness_h2s; // cloudiness, N in literature
    Array cloudiness_nh3; // cloudiness, N in literature

    Array p_dyn;                // dynamic pressure
    Array p_dynn;                // dynamic pressure
    Array p_stat;                // static pressure
    Array rho_mix;              // local mixture density from ideal gas: p_stat/(R_mix*T)

    Array rhs_t;                // auxilliar field RHS temperature
    Array rhs_u;                // auxilliar field RHS u-velocity component
    Array rhs_v;                // auxilliar field RHS v-velocity component
    Array rhs_w;                // auxilliar field RHS w-velocity component

    Array rhs_ch4;                // auxilliar field RHS water vapour
    Array rhs_ch4_cloud;            // auxilliar field RHS cloud water
    Array rhs_ch4_ice;            // auxilliar field RHS cloud water
    Array rhs_h2o;                // auxilliar field RHS water vapour
    Array rhs_h2o_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2o_ice;                // auxilliar field RHS cloud ice
    Array rhs_h2s;                // auxilliar field RHS water vapour
    Array rhs_h2s_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2s_ice;                // auxilliar field RHS cloud ice
    Array rhs_nh3;                // auxilliar field RHS nh3
    Array rhs_nh3_cloud;        // auxilliar field RHS nh3_cloud
    Array rhs_nh3_ice;            // auxilliar field RHS nh3_ice
    Array rhs_nh4sh;                // auxilliar field RHS nh4sh

    Array aux;                // auxilliar field u-velocity component
    Array aux_u;                // auxilliar field u-velocity component
    Array aux_v;                // auxilliar field v-velocity component
    Array aux_w;                // auxilliar field w-velocity component

    Array Q_Latent;                // latent heat
    Array Q_Sensible;            // sensible heat
    Array CoriolisForce;        // coriolis force
    Array CentrifugalForce;             // centrifugal force
    Array BuoyancyForce;        // buoyancy force, Boussinesque approximation
    Array PresGradForce;// pressure gradient force
    Array SeaMount;             // sea mount contour

    Array w_nh3;                // reaction rate nh3
    Array w_h2s;                // reaction rate h2s
    Array w_nh4sh;                // reaction rate nh4sh

    Array j_nh3;                // ordinary-diffusion mass flux of nh3
    Array j_h2s;                // ordinary-diffusion mass flux of h2s
    Array j_nh4sh;              // ordinary-diffusion mass flux of nh4sh

    Array jT_nh3;               // thermo-diffusion mass flux of nh3
    Array jT_h2s;               // thermo-diffusion mass flux of h2s
    Array jT_nh4sh;             // thermo-diffusion mass flux of nh4sh
};
#endif

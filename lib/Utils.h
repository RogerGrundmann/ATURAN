#ifndef _UTILS_
#define _UTILS_

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <vector>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#define logger() \
if (false) ; \
else get_logger()

namespace AtomUtils{
    using namespace std;
    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    }; 

    HemisphereCoords convert_coords(double lon, double lat);

    int lon_2_index(double lon);
    int lat_2_index(double lat);

    std::ofstream& get_logger();

    inline bool is_land(const Array& SeaMount, int i, int j, int k){
        return fabs(SeaMount.x[i][j][k] - 1) < std::numeric_limits<double>::epsilon();
    }

    inline bool is_air(const Array& SeaMount, int i, int j, int k){
        return !is_land(SeaMount,i,j,k);
    }

    inline bool is_water(const Array& SeaMount, int i, int j, int k){
        return !is_land(SeaMount,i,j,k);
    }

    inline bool is_ocean_surface(const Array& SeaMount, int i, int j, int k){
        return i==0 && !is_land(SeaMount, i, j, k);
    }

    inline bool is_land_surface(const Array& SeaMount, int i, int j, int k){
        return is_land(SeaMount, i, j, k) && !is_land(SeaMount, i+1, j, k);
    }

    inline double parabola(double x){
        return x*x - 2*x;
    }

    inline double Agnesi(double x, double equator){
        return pow(equator, 3.) // Versiera di Agnesi
            / ( pow(equator, 2.) + pow(x, 2.));
    }

    inline double parabola_interp(double lower_bound, double upper_bound, double x){
        return upper_bound + (upper_bound-lower_bound)*(x*x - 2*x);
    }

    inline bool is_east_coast(const Array& SeaMount, int j, int k){
        if(k == SeaMount.get_km()-1 ) return false;//on grid boundary
        return is_land(SeaMount, 0, j, k) && !is_land(SeaMount, 0, j, k+1);
    }

    inline bool is_west_coast(const Array& SeaMount, int j, int k){
        if( k == 0 ) return false;//on grid boundary
        return is_land(SeaMount, 0, j, k) && !is_land(SeaMount, 0, j, k-1);
    }

    inline bool is_north_coast(const Array& SeaMount, int j, int k){
        if( j == 0 ) return false;//on grid boundary
        return is_land(SeaMount, 0, j, k) && !is_land(SeaMount, 0, j-1, k);
    }

    inline bool is_south_coast(const Array& SeaMount, int j, int k){
        if(j == SeaMount.get_jm()-1 ) return false;//on grid boundary
        return is_land(SeaMount, 0, j, k) && !is_land(SeaMount, 0, j+1, k);
    }

    //change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
    void move_data(double* data, int len);
    void move_data(std::vector<int>& data, int len);

    void move_data_to_new_arrays(int im, int jm, int km, double coeff, std::vector<Array*>& arrays, 
                        std::vector<Array*>& new_arrays);

    void move_data_to_new_arrays( int jm, int km, double coeff, std::vector<Array*>& arrays, 
                        std::vector<Array*>& new_arrays, int i = 0);

    std::tuple<double, int, int, int>
    max_diff(int i, int j, int k, const Array &a1, const Array &a2);

    void load_map_from_file(const std::string& fn, std::map<float, float>& m);

    // integration tools
    double simpson(int n1, int n2, double dstep, Array_1D &value);
    double trapezoidal(int n1, int n2, double dstep, Array_1D &value);
    double rectangular(int n1, int n2, double dstep, Array_1D &value);

    double C_Dalton (double coeff_Dalton, double u_0, double v, double w);

    void read_IC(const string& fn, double** a, int jm, int km);

    void smooth_steps(int k, int im, int jm, Array &value);

    void smooth_tropopause(int jm, std::vector<double> &value);

    void fft_gaussian_filter(Array& data, int sigma);

    void fft_gaussian_filter_3d(Array& data, int sigma);

    void mirror_padding(double* data, size_t i_len, size_t p_len);

    void fft_gaussian_filter(double* _data, double* kernel, size_t len);

    template<class T>
    T clamp(T value, T lo, T hi){
        return value < lo ? lo : (value > hi ? hi : value);
    }

    // Shapiro (1-2-1) wiggle-damping filter.
    // f_new[n] = f[n] + strength*0.25*(f[n-1] - 2*f[n] + f[n+1])
    // k (longitude): periodic wrap; j (latitude): clamped at poles; i (vertical): clamped.
    // i_surface: first fluid cell per (j,k); pass nullptr to filter all cells.
    void damp_wiggles(Array& field,
                      const std::vector<std::vector<int>>* i_surface,
                      bool along_i, bool along_j, bool along_k,
                      double strength = 1.0,
                      int passes = 1);

    // 2-D overload (Array_2D, axes j and k only).
    void damp_wiggles(Array_2D& field,
                      bool along_j, bool along_k,
                      double strength = 1.0,
                      int passes = 1);

    template<class T>
    void set_values(T* a, T value, int len){
        for(int i = 0 ; i < len; i++){
            a[i] = value;
        }
    }

    template<class T>
    void chessboard_grid(T** a, int j, int k, int jm, int km)
    {
        T value_1 = 0.9, value_2 = 1.1;
        T* line_1 = new T[km];
        T* line_2 = new T[km];
        int i = 0;
        bool flip = true;
        while(i+k < km){
            if(flip){
                set_values(line_1+i, value_1, k);
                set_values(line_2+i, value_2, k);
            }else{
                set_values(line_1+i, value_2, k);
                set_values(line_2+i, value_1, k);
            }
            i+=k;
            flip = !flip;
        }
        if(i<km-1){
            if(flip){
                set_values(line_1+i, value_1, km-1-i);
                set_values(line_2+i, value_2, km-1-i);
            }else{
                set_values(line_1+i, value_2, km-1-i);
                set_values(line_2+i, value_1, km-1-i);
            }
        }
        int cnt=0;
        flip = true;
        for(i = 0; i<jm; i++){
            if(flip)
                memcpy(a[i], line_1, km*sizeof(T));
            else
                memcpy(a[i], line_2, km*sizeof(T));
            cnt++;
            if(cnt==j){
                flip = !flip;
                cnt = 0;
            }
        }

        delete[] line_1;
        delete[] line_2;
        return;
    }
}

#endif

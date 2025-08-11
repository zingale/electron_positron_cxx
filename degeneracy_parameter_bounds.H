#ifndef DEGENERACY_PARAMETER_BOUNDS_H
#define DEGENERACY_PARAMETER_BOUNDS_H

#include <algorithm>
#include <array>

#include "real_type.H"
#include "util.H"

using namespace literals;

namespace bounds {

    constexpr std::array<real_t, 10> Ts{
        {1.e+03_rt, 1.e+04_rt, 1.e+05_rt, 1.e+06_rt, 1.e+07_rt,
         1.e+08_rt, 1.e+09_rt, 1.e+10_rt, 1.e+11_rt, 1.e+12_rt}};

    constexpr std::array<real_t, 24> rhoYes{
        {1.e-12_rt, 1.e-11_rt, 1.e-10_rt, 1.e-09_rt, 1.e-08_rt, 1.e-07_rt,
         1.e-06_rt, 1.e-05_rt, 1.e-04_rt, 1.e-03_rt, 1.e-02_rt, 1.e-01_rt,
         1.e+00_rt, 1.e+01_rt, 1.e+02_rt, 1.e+03_rt, 1.e+04_rt, 1.e+05_rt,
         1.e+06_rt, 1.e+07_rt, 1.e+08_rt, 1.e+09_rt, 1.e+10_rt, 1.e+11_rt}};

    // degeneracy parameters for electrons + positrons

    constexpr std::array<std::array<real_t, 10>, 24> eta_pos{
        {{   -19.35_rt,    -22.81_rt,    -26.26_rt,    -29.71_rt,    -33.17_rt,    -36.65_rt,    -40.37_rt,    -45.29_rt,    -51.75_rt,     -58.6_rt},
         {   -17.05_rt,     -20.5_rt,    -23.96_rt,    -27.41_rt,    -30.87_rt,    -34.35_rt,    -38.06_rt,    -42.99_rt,    -49.44_rt,     -56.3_rt},
         {   -14.75_rt,     -18.2_rt,    -21.65_rt,    -25.11_rt,    -28.56_rt,    -32.05_rt,    -35.76_rt,    -40.69_rt,    -47.14_rt,       -54_rt},
         {   -12.44_rt,     -15.9_rt,    -19.35_rt,    -22.81_rt,    -26.26_rt,    -29.74_rt,    -33.46_rt,    -38.39_rt,    -44.84_rt,    -51.69_rt},
         {   -10.14_rt,    -13.59_rt,    -17.05_rt,     -20.5_rt,    -23.96_rt,    -27.44_rt,    -31.16_rt,    -36.08_rt,    -42.54_rt,    -49.39_rt},
         {   -7.838_rt,    -11.29_rt,    -14.75_rt,     -18.2_rt,    -21.66_rt,    -25.14_rt,    -28.85_rt,    -33.78_rt,    -40.23_rt,    -47.09_rt},
         {   -5.534_rt,     -8.99_rt,    -12.44_rt,     -15.9_rt,    -19.35_rt,    -22.84_rt,    -26.55_rt,    -31.48_rt,    -37.93_rt,    -44.79_rt},
         {   -3.219_rt,    -6.687_rt,    -10.14_rt,     -13.6_rt,    -17.05_rt,    -20.53_rt,    -24.25_rt,    -29.17_rt,    -35.63_rt,    -42.48_rt},
         {  -0.7919_rt,     -4.38_rt,    -7.838_rt,    -11.29_rt,    -14.75_rt,    -18.23_rt,    -21.95_rt,    -26.87_rt,    -33.33_rt,    -40.18_rt},
         {    2.697_rt,    -2.038_rt,    -5.534_rt,     -8.99_rt,    -12.45_rt,    -15.93_rt,    -19.64_rt,    -24.57_rt,    -31.02_rt,    -37.88_rt},
         {    13.95_rt,    0.6542_rt,    -3.219_rt,    -6.687_rt,    -10.14_rt,    -13.63_rt,    -17.34_rt,    -22.27_rt,    -28.72_rt,    -35.58_rt},
         {       65_rt,     6.369_rt,   -0.7919_rt,     -4.38_rt,    -7.841_rt,    -11.32_rt,    -15.04_rt,    -19.96_rt,    -26.42_rt,    -33.27_rt},
         {    301.8_rt,     30.15_rt,     2.697_rt,    -2.038_rt,    -5.538_rt,    -9.021_rt,    -12.73_rt,    -17.66_rt,    -24.11_rt,    -30.97_rt},
         {     1400_rt,       140_rt,     13.95_rt,    0.6536_rt,    -3.222_rt,    -6.718_rt,    -10.43_rt,    -15.36_rt,    -21.81_rt,    -28.67_rt},
         {     6498_rt,     649.8_rt,     64.96_rt,     6.365_rt,   -0.7958_rt,    -4.412_rt,     -8.13_rt,    -13.06_rt,    -19.51_rt,    -26.36_rt},
         { 3.01e+04_rt,      3010_rt,       301_rt,     30.07_rt,     2.686_rt,    -2.071_rt,    -5.826_rt,    -10.75_rt,    -17.21_rt,    -24.06_rt},
         {1.384e+05_rt, 1.384e+04_rt,      1384_rt,     138.4_rt,     13.78_rt,    0.6006_rt,    -3.515_rt,    -8.451_rt,     -14.9_rt,    -21.76_rt},
         {6.179e+05_rt, 6.179e+04_rt,      6179_rt,     617.9_rt,     61.78_rt,     6.009_rt,    -1.133_rt,    -6.148_rt,     -12.6_rt,    -19.46_rt},
         {2.493e+06_rt, 2.493e+05_rt, 2.493e+04_rt,      2493_rt,     249.3_rt,     24.88_rt,     1.869_rt,    -3.842_rt,     -10.3_rt,    -17.15_rt},
         {8.257e+06_rt, 8.257e+05_rt, 8.257e+04_rt,      8257_rt,     825.7_rt,     82.55_rt,         8_rt,    -1.505_rt,    -7.997_rt,    -14.85_rt},
         {2.246e+07_rt, 2.246e+06_rt, 2.246e+05_rt, 2.246e+04_rt,      2246_rt,     224.6_rt,     22.35_rt,     1.083_rt,    -5.694_rt,    -12.55_rt},
         {5.419e+07_rt, 5.419e+06_rt, 5.419e+05_rt, 5.419e+04_rt,      5419_rt,     541.9_rt,     54.13_rt,      4.87_rt,    -3.387_rt,    -10.25_rt},
         {1.231e+08_rt, 1.231e+07_rt, 1.231e+06_rt, 1.231e+05_rt, 1.231e+04_rt,      1231_rt,     123.1_rt,     12.05_rt,    -1.046_rt,    -7.944_rt},
         {2.718e+08_rt, 2.718e+07_rt, 2.718e+06_rt, 2.718e+05_rt, 2.718e+04_rt,      2718_rt,     271.8_rt,     27.06_rt,     1.562_rt,    -5.641_rt}}
    };

    // degeneracy parameters for just electrons

    constexpr std::array<std::array<real_t, 10>, 24> eta_ele{
        {{   -19.35_rt,    -22.81_rt,    -26.26_rt,    -29.71_rt,    -33.17_rt,    -36.65_rt,    -40.37_rt,    -45.29_rt,    -51.75_rt,     -58.6_rt},
         {   -17.05_rt,     -20.5_rt,    -23.96_rt,    -27.41_rt,    -30.87_rt,    -34.35_rt,    -38.06_rt,    -42.99_rt,    -49.44_rt,     -56.3_rt},
         {   -14.75_rt,     -18.2_rt,    -21.65_rt,    -25.11_rt,    -28.56_rt,    -32.05_rt,    -35.76_rt,    -40.69_rt,    -47.14_rt,       -54_rt},
         {   -12.44_rt,     -15.9_rt,    -19.35_rt,    -22.81_rt,    -26.26_rt,    -29.74_rt,    -33.46_rt,    -38.39_rt,    -44.84_rt,    -51.69_rt},
         {   -10.14_rt,    -13.59_rt,    -17.05_rt,     -20.5_rt,    -23.96_rt,    -27.44_rt,    -31.16_rt,    -36.08_rt,    -42.54_rt,    -49.39_rt},
         {   -7.838_rt,    -11.29_rt,    -14.75_rt,     -18.2_rt,    -21.66_rt,    -25.14_rt,    -28.85_rt,    -33.78_rt,    -40.23_rt,    -47.09_rt},
         {   -5.534_rt,     -8.99_rt,    -12.44_rt,     -15.9_rt,    -19.35_rt,    -22.84_rt,    -26.55_rt,    -31.48_rt,    -37.93_rt,    -44.79_rt},
         {   -3.219_rt,    -6.687_rt,    -10.14_rt,     -13.6_rt,    -17.05_rt,    -20.53_rt,    -24.25_rt,    -29.17_rt,    -35.63_rt,    -42.48_rt},
         {  -0.7919_rt,     -4.38_rt,    -7.838_rt,    -11.29_rt,    -14.75_rt,    -18.23_rt,    -21.95_rt,    -26.87_rt,    -33.33_rt,    -40.18_rt},
         {    2.697_rt,    -2.038_rt,    -5.534_rt,     -8.99_rt,    -12.45_rt,    -15.93_rt,    -19.64_rt,    -24.57_rt,    -31.02_rt,    -37.88_rt},
         {    13.95_rt,    0.6542_rt,    -3.219_rt,    -6.687_rt,    -10.14_rt,    -13.63_rt,    -17.34_rt,    -22.27_rt,    -28.72_rt,    -35.58_rt},
         {       65_rt,     6.369_rt,   -0.7919_rt,     -4.38_rt,    -7.841_rt,    -11.32_rt,    -15.04_rt,    -19.96_rt,    -26.42_rt,    -33.27_rt},
         {    301.8_rt,     30.15_rt,     2.697_rt,    -2.038_rt,    -5.538_rt,    -9.021_rt,    -12.73_rt,    -17.66_rt,    -24.11_rt,    -30.97_rt},
         {     1400_rt,       140_rt,     13.95_rt,    0.6536_rt,    -3.222_rt,    -6.718_rt,    -10.43_rt,    -15.36_rt,    -21.81_rt,    -28.67_rt},
         {     6498_rt,     649.8_rt,     64.96_rt,     6.365_rt,   -0.7958_rt,    -4.412_rt,     -8.13_rt,    -13.06_rt,    -19.51_rt,    -26.36_rt},
         { 3.01e+04_rt,      3010_rt,       301_rt,     30.07_rt,     2.686_rt,    -2.071_rt,    -5.826_rt,    -10.75_rt,    -17.21_rt,    -24.06_rt},
         {1.384e+05_rt, 1.384e+04_rt,      1384_rt,     138.4_rt,     13.78_rt,    0.6006_rt,    -3.515_rt,    -8.451_rt,     -14.9_rt,    -21.76_rt},
         {6.179e+05_rt, 6.179e+04_rt,      6179_rt,     617.9_rt,     61.78_rt,     6.009_rt,    -1.133_rt,    -6.148_rt,     -12.6_rt,    -19.46_rt},
         {2.493e+06_rt, 2.493e+05_rt, 2.493e+04_rt,      2493_rt,     249.3_rt,     24.88_rt,     1.869_rt,    -3.842_rt,     -10.3_rt,    -17.15_rt},
         {8.257e+06_rt, 8.257e+05_rt, 8.257e+04_rt,      8257_rt,     825.7_rt,     82.55_rt,         8_rt,    -1.505_rt,    -7.997_rt,    -14.85_rt},
         {2.246e+07_rt, 2.246e+06_rt, 2.246e+05_rt, 2.246e+04_rt,      2246_rt,     224.6_rt,     22.35_rt,     1.083_rt,    -5.694_rt,    -12.55_rt},
         {5.419e+07_rt, 5.419e+06_rt, 5.419e+05_rt, 5.419e+04_rt,      5419_rt,     541.9_rt,     54.13_rt,      4.87_rt,    -3.387_rt,    -10.25_rt},
         {1.231e+08_rt, 1.231e+07_rt, 1.231e+06_rt, 1.231e+05_rt, 1.231e+04_rt,      1231_rt,     123.1_rt,     12.05_rt,    -1.046_rt,    -7.944_rt},
         {2.718e+08_rt, 2.718e+07_rt, 2.718e+06_rt, 2.718e+05_rt, 2.718e+04_rt,      2718_rt,     271.8_rt,     27.06_rt,     1.562_rt,    -5.641_rt}}
    };

    template <typename T>
    inline
    std::pair<T, T> get_eta_bounds(T rhoYe0, T temp0, bool include_positrons=true) {

        // find the index into the table such that
        //   rhoYes[ir] <= rhoYe0 < rhoYe[ir+1]
        // and likewise for T

        const auto *const ip = std::upper_bound(rhoYes.cbegin(), rhoYes.cend(), rhoYe0);
        int ir{-1};
        if (ip != rhoYes.cend()) {
            ir = std::clamp(static_cast<int>(std::distance(rhoYes.cbegin(), ip)) - 1,
                              1, static_cast<int>(rhoYes.size())-3);
        }

        const auto *const ip2 = std::upper_bound(Ts.cbegin(), Ts.cend(), temp0);
        int it{-1};
        if (ip2 != Ts.cend()) {
            it = std::clamp(static_cast<int>(std::distance(Ts.cbegin(), ip2)) - 1,
                               1, static_cast<int>(Ts.size())-3);
        }

        if (it < 0 || ir < 0) {
            util::safe_println("thermodynamic state (T = {}, rhoYe = {}) out of bounds", temp0, rhoYe0);
            abort();
        }

        T eta_min{};
        T eta_max{};
        if (include_positrons) {
            eta_min = std::min(
                {eta_pos[ir-1][it-1], eta_pos[ir-1][it], eta_pos[ir-1][it+1], eta_pos[ir-1][it+2],
                 eta_pos[ir][it-1],   eta_pos[ir][it],   eta_pos[ir][it+1],   eta_pos[ir][it+2],
                 eta_pos[ir+1][it-1], eta_pos[ir+1][it], eta_pos[ir+1][it+1], eta_pos[ir+1][it+2],
                 eta_pos[ir+2][it-1], eta_pos[ir+2][it], eta_pos[ir+2][it+1], eta_pos[ir+2][it+2]});
            eta_max = std::max(
                {eta_pos[ir-1][it-1], eta_pos[ir-1][it], eta_pos[ir-1][it+1], eta_pos[ir-1][it+2],
                 eta_pos[ir][it-1],   eta_pos[ir][it],   eta_pos[ir][it+1],   eta_pos[ir][it+2],
                 eta_pos[ir+1][it-1], eta_pos[ir+1][it], eta_pos[ir+1][it+1], eta_pos[ir+1][it+2],
                 eta_pos[ir+2][it-1], eta_pos[ir+2][it], eta_pos[ir+2][it+1], eta_pos[ir+2][it+2]});

        } else {
            eta_min = std::min(
                {eta_ele[ir-1][it-1], eta_ele[ir-1][it], eta_ele[ir-1][it+1], eta_ele[ir-1][it+2],
                 eta_ele[ir][it-1],   eta_ele[ir][it],   eta_ele[ir][it+1],   eta_ele[ir][it+2],
                 eta_ele[ir+1][it-1], eta_ele[ir+1][it], eta_ele[ir+1][it+1], eta_ele[ir+1][it+2],
                 eta_ele[ir+2][it-1], eta_ele[ir+2][it], eta_ele[ir+2][it+1], eta_ele[ir+2][it+2]});
            eta_max = std::max(
                {eta_ele[ir-1][it-1], eta_ele[ir-1][it], eta_ele[ir-1][it+1], eta_ele[ir-1][it+2],
                 eta_ele[ir][it-1],   eta_ele[ir][it],   eta_ele[ir][it+1],   eta_ele[ir][it+2],
                 eta_ele[ir+1][it-1], eta_ele[ir+1][it], eta_ele[ir+1][it+1], eta_ele[ir+1][it+2],
                 eta_ele[ir+2][it-1], eta_ele[ir+2][it], eta_ele[ir+2][it+1], eta_ele[ir+2][it+2]});
        }

        return {eta_min, eta_max};
    }

}

#endif

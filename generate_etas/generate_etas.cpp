#include <array>
#include <ranges>
#include <iostream>
#include <iomanip>

#include "real_type.H"
#include "electron_positron.H"
#include "util.H"

using namespace literals;

const std::array<real_t, 19> Ts_v{
    {1.e3_rt, 3.e3_rt, 1.e4_rt, 3.e4_rt, 1.e5_rt, 3.e5_rt, 1.e6_rt, 3.e6_rt, 1.e7_rt, 3.e7_rt,
     1.e8_rt, 3.e8_rt, 1.e9_rt, 3.e9_rt, 1.e10_rt, 3.e10_rt, 1.e11_rt, 3.e11_rt, 1.e12_rt}};

const std::array<real_t, 25> rhoYes_v{
    {1.e-12_rt, 1.e-11_rt, 1.e-10_rt, 1.e-09_rt, 1.e-08_rt, 1.e-07_rt,
     1.e-06_rt, 1.e-05_rt, 1.e-04_rt, 1.e-03_rt, 1.e-02_rt, 1.e-01_rt,
     1.e+00_rt, 1.e+01_rt, 1.e+02_rt, 1.e+03_rt, 1.e+04_rt, 1.e+05_rt,
     1.e+06_rt, 1.e+07_rt, 1.e+08_rt, 1.e+09_rt, 1.e+10_rt, 1.e+11_rt,
     1.e+12_rt}};

auto main() -> int
{

    const real_t Ye(1.0);

    for (auto [ir, rho] : std::views::enumerate(rhoYes_v)) {
        for (auto [it, T] : std::views::enumerate(Ts_v)) {

            real_t eta{};
            real_t beta = C::dbeta_dT * T;
            real_t n_e_net = rho * Ye * C::N_A;

            try {
                eta = brent<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto n_e = n_e_constraint(_eta, beta);
                        auto n_pos = n_p_constraint(_eta, beta);
                        return n_e_net - (n_e - n_pos);
                    }, -100.0_rt, 1.e12_rt);
            } catch (const std::out_of_range& e) {
                util::red_println("bounds failed for rho = {:8.3g}, T = {:8.3g}", rho, T);
            }

            if (it == 0) {
                if (ir == 0) {
                    std::cout << "{{";
                } else {
                    std::cout << " {";
                }
            }

            std::cout << util::format("{:12.8g}", eta);
            if (it < Ts_v.size()-1) {
                std::cout << ", ";
            }

            if (it == Ts_v.size()-1) {
                if (ir == rhoYes_v.size()-1) {
                    std::cout << "}}\n";
                } else {
                    std::cout << "},\n";
                }
            }
        }
    }

}

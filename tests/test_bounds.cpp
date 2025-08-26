#include "degeneracy_parameter_bounds.H"
#include "real_type.H"
#include "util.H"
#include "electron_positron.H"

const std::array<real_t, 5> Ts{1.e4_rt, 1.e5_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

using namespace literals;

auto main() -> int
{

    const real_t Ye(0.5);

    for (auto T : Ts) {
        for (auto rho : rhos) {

            auto [eta_min, eta_max] = bounds::get_eta_bounds(rho * Ye, T);

            real_t eta{};
            real_t beta = C::dbeta_dT * T;
            real_t n_e_net = rho * Ye * C::N_A;

            try {
                eta = brent<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto n_e = n_e_constraint(_eta, beta);
                        auto n_pos = n_p_constraint(_eta, beta);
                        return n_e_net - (n_e - n_pos);
                    }, eta_min, eta_max);
                util::green_println("bounds passed for rho = {:8.3g}, T = {:8.3g}; eta = {:9.2f}", rho, T, eta);

            } catch (const std::out_of_range& e) {
                util::red_println("bounds failed for rho = {:8.3g}, T = {:8.3g}", rho, T);
            }
        }
    }
}


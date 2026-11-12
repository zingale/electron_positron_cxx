#include <array>
#include <stdexcept>

#include "brent.H"
#include "degeneracy_parameter_bounds.H"
#include "real_type.H"
#include "util.H"
#include "electron_positron.H"
#include "fundamental_constants.H"

const std::array<real_t, 6> Ts{2.e3_rt, 1.e4_rt, 1.e5_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 6> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt, 1.e11_rt};

using namespace literals;

auto main() -> int
{

    const real_t Ye(0.5);

    // Exercise the first and last temperature intervals at a density near
    // the upper end of the table, where using the wrong interval misses eta.
    const real_t edge_rhoYe = 0.5_rt * bounds::rhoYes.back();
    const std::array<real_t, 3> edge_temps{
        bounds::Ts.front(),
        0.5_rt * (bounds::Ts[0] + bounds::Ts[1]),
        0.5_rt * (bounds::Ts[bounds::Ts.size()-2] + bounds::Ts.back())};

    for (auto temp : edge_temps) {
        const auto [eta_min, eta_max] = bounds::get_eta_bounds(edge_rhoYe, temp);
        const real_t beta = C::dbeta_dT * temp;
        const real_t n_e_net = edge_rhoYe * C::N_A;
        const auto residual = [=] (real_t eta) -> real_t {
            return n_e_net - (n_e_constraint(eta, beta) - n_p_constraint(eta, beta));
        };

        // Charge neutrality requires n_e - n_pos = rhoYe * N_A, so the
        // residual is zero at the desired eta.  Increasing eta increases
        // n_e - n_pos: the lower bound must give too few net electrons
        // (residual >= 0), and the upper bound too many (residual <= 0).
        // Check the table's bracket directly because get_eta() could hide
        // a bad bracket by successfully solving with its wider fallback.
        if (!(residual(eta_min) >= 0.0_rt && residual(eta_max) <= 0.0_rt)) {
            throw std::runtime_error("Temperature edge bracket misses eta");
        }

        const auto eta = get_eta(edge_rhoYe, temp);
        if (!(mp::abs(residual(eta) / n_e_net) < 1.e-10_rt)) {
            throw std::runtime_error("Temperature edge solve violates charge neutrality");
        }
    }

    for (auto T : Ts) {
        for (auto rho : rhos) {

            auto [eta_min, eta_max] = bounds::get_eta_bounds(rho * Ye, T);

            real_t eta{};
            const real_t beta = C::dbeta_dT * T;
            const real_t n_e_net = rho * Ye * C::N_A;

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


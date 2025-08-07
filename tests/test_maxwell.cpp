#include <print>

#include <array>
#include <cmath>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"

constexpr std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
constexpr std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

// number density

void
test_maxwell() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    {

        std::println("testing p = ρ² ∂e/∂ρ|ᴛ + T ∂p/∂T|ᵨ");

        for (auto T : Ts) {
            for (auto rho : rhos) {
                auto es = eos.pe_state(rho, T, Ye);

                real_t ptot = es.p_e + es.p_pos;
                real_t de_drho = es.dee_drho + es.dep_drho;
                real_t dp_dT = es.dpe_dT + es.dpp_dT;

                real_t term = std::abs((ptot - (rho * rho * de_drho + T * dp_dT)) / ptot);

                std::println("ρ = {:8.3g} T = {:8.3g}:  p⁻ + p⁺ = {:15.8g},  p = ρ² ∂e/∂ρ|ᴛ + T ∂p/∂T|ᵨ error = {:15.5g}",
                             rho, T, ptot, term);

            }
        }
    }

    {
        std::println("");
        std::println("testing ∂e/∂T|ᵨ = T ∂s/∂T|ᵨ");

        for (auto T : Ts) {
            for (auto rho : rhos) {
                auto es = eos.pe_state(rho, T, Ye);

                real_t de_dT = es.dee_dT + es.dep_dT;
                real_t ds_dT = es.dse_dT + es.dsp_dT;

                real_t term = std::abs((de_dT - T * ds_dT) / de_dT);

                std::println("ρ = {:8.3g} T = {:8.3g}:  ∂e/∂T|ᵨ = {:15.8g},  ∂e/∂T|ᵨ = T ∂s/∂T|ᵨ error = {:15.5g}",
                             rho, T, de_dT, term);

            }
        }
    }

    {
        std::println("");
        std::println("testing -∂s/∂ρ|ᴛ = 1/ρ² ∂p/∂T|ᵨ");

        for (auto T : Ts) {
            for (auto rho : rhos) {
                auto es = eos.pe_state(rho, T, Ye);

                real_t ds_drho = es.dse_drho + es.dsp_drho;
                real_t dp_dT = es.dpe_dT + es.dpp_dT;

                real_t term = std::abs((ds_drho + 1.0_rt / (rho * rho) * dp_dT) / ds_drho);

                std::println("ρ = {:8.3g} T = {:8.3g}:  ∂s/∂ρ|ᴛ = {:15.8g},  -∂s/∂ρ|ᴛ = 1/ρ² ∂p/∂T|ᵨ error = {:15.5g}",
                             rho, T, dp_dT, term);

            }
        }
    }
}


int main() {

    test_maxwell();

}

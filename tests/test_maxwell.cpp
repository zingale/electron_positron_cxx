#include <print>

#include <array>

#include "real_type.H"
#include "maxwell_relations.H"
#include "util.H"

constexpr std::array<real_t, 5> Ts{1.e4_rt, 1.e5_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
constexpr std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

// number density



int main() {

    constexpr real_t Ye{0.5_rt};

    util::green_println("testing p = ρ² ∂e/∂ρ|ᴛ + T ∂p/∂T|ᵨ");

    for (auto T : Ts) {
        for (auto rho : rhos) {

            auto [scale, error] = maxwell_1<real_t>(rho, T, Ye);

            std::println("ρ = {:8.3g} T = {:8.3g}:  p⁻ + p⁺ = {:15.8g},  p = ρ² ∂e/∂ρ|ᴛ + T ∂p/∂T|ᵨ error = {:15.5g}",
                         rho, T, scale, error);

        }
    }

    std::println("");
    util::green_println("testing ∂e/∂T|ᵨ = T ∂s/∂T|ᵨ");

    for (auto T : Ts) {
        for (auto rho : rhos) {

            auto [scale, error] = maxwell_2<real_t>(rho, T, Ye);

            std::println("ρ = {:8.3g} T = {:8.3g}:  ∂e/∂T|ᵨ = {:15.8g},  ∂e/∂T|ᵨ = T ∂s/∂T|ᵨ error = {:15.5g}",
                         rho, T, scale, error);

        }
    }

    std::println("");
    util::green_println("testing -∂s/∂ρ|ᴛ = 1/ρ² ∂p/∂T|ᵨ");

    for (auto T : Ts) {
        for (auto rho : rhos) {

            auto [scale, error] = maxwell_3<real_t>(rho, T, Ye);

            std::println("ρ = {:8.3g} T = {:8.3g}:  ∂s/∂ρ|ᴛ = {:15.8g},  -∂s/∂ρ|ᴛ = 1/ρ² ∂p/∂T|ᵨ error = {:15.5g}",
                         rho, T, scale, error);

        }
    }

}

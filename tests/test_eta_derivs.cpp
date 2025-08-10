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
test_eta_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("testing ∂η/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.eta;
                }, rho, drho);

            real_t err = std::abs(es.deta_drho - deriv) / std::abs(es.deta_drho);
            std::println("ρ = {:10.3g} T = {:10.3g} η = {:15.8}: ∂η/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.eta, es.deta_drho, err);
        }
    }
}

void
test_eta_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂η/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.eta;
                }, T, dT);

            real_t err{};
            err = std::abs(es.deta_dT - deriv) / std::abs(es.deta_dT);
            std::println("ρ = {:10.3g} T = {:10.3g} η = {:15.8}: ∂η/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.eta, es.deta_dT, err);
        }
    }
}

void
test_eta_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("testing ∂²η/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.eta;
                }, rho, drho);

            real_t err = std::abs(es.d2eta_drho2 - deriv) / std::abs(es.d2eta_drho2);
            std::println("ρ = {:10.3g} T = {:10.3g} η = {:15.8}: ∂²η/∂ρ² = {:15.8g}, error = {:15.5g}",
                         rho, T, es.eta, es.d2eta_drho2, err);
        }
    }
}

void
test_eta_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂²η/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.eta;
                }, T, dT);

            real_t err{};
            err = std::abs(es.d2eta_dT2 - deriv) / std::abs(es.d2eta_dT2);
            std::println("ρ = {:10.3g} T = {:10.3g} η = {:15.8}: ∂²η/∂T² = {:15.8g}, error = {:15.5g}",
                         rho, T, es.eta, es.d2eta_dT2, err);
        }
    }
}


int main() {

    test_eta_rho_derivs();
    test_eta_T_derivs();

    test_eta_rho2_derivs();
    test_eta_T2_derivs();

}

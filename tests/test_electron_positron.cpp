#include <print>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"

constexpr std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
constexpr std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

void
test_ne_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("testing ∂n⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_e;
                }, rho, drho);

            real_t err = std::abs(es.dne_drho - deriv) / std::abs(es.dne_drho);
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂n⁻/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dne_drho, err);
        }
    }
}

void
test_ne_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂n⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _T) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.n_e;
                }, T, dT);

            real_t err{};
            if (es.dne_dT == 0.0_rt) {
                real_t scale = es.n_e / T;
                err = std::abs(es.dne_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dne_dT - deriv) / std::abs(es.dne_dT);
            }
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂n⁻/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dne_dT, err);
        }
    }
}

void
test_np_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂n⁺/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dnp_drho == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_pos;
                }, rho, drho);

            real_t err = std::abs(es.dnp_drho - deriv) / std::abs(es.dnp_drho);
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂n⁺/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dnp_drho, err);
        }
    }

}

void
test_np_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂n⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dnp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _T) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.n_pos;
                }, T, dT);

            real_t err{};
            if (es.dnp_dT == 0.0_rt) {
                real_t scale = es.n_pos / T;
                err = std::abs(es.dnp_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dnp_dT - deriv) / std::abs(es.dnp_dT);
            }
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂n⁺/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dnp_dT, err);
        }
    }
}

void
test_pe_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂p⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_e;
                }, rho, drho);

            real_t err = std::abs(es.dpe_drho - deriv) / std::abs(es.dpe_drho);
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂p⁻/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpe_drho, err);
        }
    }
}

void
test_pe_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂p⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _T) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.p_e;
                }, T, dT);

            real_t err{};
            if (es.dpe_dT == 0.0_rt) {
                real_t scale = es.p_e / T;
                err = std::abs(es.dpe_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dpe_dT - deriv) / std::abs(es.dpe_dT);
            }
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂p⁻/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpe_dT, err);
        }
    }
}

void
test_pp_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂p⁺/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dpp_drho == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_pos;
                }, rho, drho);

            real_t err = std::abs(es.dpp_drho - deriv) / std::abs(es.dpp_drho);
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂p⁺/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpp_drho, err);
        }
    }

}

void
test_pp_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    std::println("testing ∂p⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dpp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _T) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.p_pos;
                }, T, dT);

            real_t err{};
            if (es.dpp_dT == 0.0_rt) {
                real_t scale = es.p_pos / T;
                err = std::abs(es.dpp_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dpp_dT - deriv) / std::abs(es.dpp_dT);
            }
            std::println("ρ = {:10.3g} T = {:10.3g}, ∂p⁺/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpp_dT, err);
        }
    }
}


int main() {

    test_ne_rho_derivs();
    test_ne_T_derivs();

    test_np_rho_derivs();
    test_np_T_derivs();

    test_pe_rho_derivs();
    test_pe_T_derivs();

    test_pp_rho_derivs();
    test_pp_T_derivs();

}

#include <iostream>
#include <print>

#include "electron_positron.H"
#include "difference_utils.H"

using RealT = long double;


constexpr std::array<RealT, 4> Ts{1.e4L, 1.e6L, 1.e8L, 5.e9L};
constexpr std::array<RealT, 5> rhos{1.e-2L, 1.e2L, 1.e5L, 1.e7L, 5.e9L};

void
test_ne_rho_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("testing ∂n⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _rho) -> RealT
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_e;
                }, rho, drho);

            RealT err = std::abs(es.dne_drho - deriv) / std::abs(es.dne_drho);
            std::println("ρ = {:10} T = {:10}, ∂n⁻/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dne_drho, err);
        }
    }
}

void
test_ne_T_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂n⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _T) -> RealT
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.n_e;
                }, T, dT);

            RealT err{};
            if (es.dne_dT == 0.0L) {
                RealT scale = es.n_e / T;
                err = std::abs(es.dne_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dne_dT - deriv) / std::abs(es.dne_dT);
            }
            std::println("ρ = {:10} T = {:10}, ∂n⁻/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dne_dT, err);
        }
    }
}

void
test_np_rho_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂n⁺/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dnp_drho == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _rho) -> RealT
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_pos;
                }, rho, drho);

            RealT err = std::abs(es.dnp_drho - deriv) / std::abs(es.dnp_drho);
            std::println("ρ = {:10} T = {:10}, ∂n⁺/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dnp_drho, err);
        }
    }

}

void
test_np_T_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂n⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dnp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _T) -> RealT
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.n_pos;
                }, T, dT);

            RealT err{};
            if (es.dnp_dT == 0.0L) {
                RealT scale = es.n_pos / T;
                err = std::abs(es.dnp_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dnp_dT - deriv) / std::abs(es.dnp_dT);
            }
            std::println("ρ = {:10} T = {:10}, ∂n⁺/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dnp_dT, err);
        }
    }
}

void
test_pe_rho_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂p⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _rho) -> RealT
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_e;
                }, rho, drho);

            RealT err = std::abs(es.dpe_drho - deriv) / std::abs(es.dpe_drho);
            std::println("ρ = {:10} T = {:10}, ∂p⁻/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpe_drho, err);
        }
    }
}

void
test_pe_T_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂p⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _T) -> RealT
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.p_e;
                }, T, dT);

            RealT err{};
            if (es.dpe_dT == 0.0L) {
                RealT scale = es.p_e / T;
                err = std::abs(es.dpe_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dpe_dT - deriv) / std::abs(es.dpe_dT);
            }
            std::println("ρ = {:10} T = {:10}, ∂p⁻/∂T = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpe_dT, err);
        }
    }
}

void
test_pp_rho_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂p⁺/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dpp_drho == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _rho) -> RealT
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_pos;
                }, rho, drho);

            RealT err = std::abs(es.dpp_drho - deriv) / std::abs(es.dpp_drho);
            std::println("ρ = {:10} T = {:10}, ∂p⁺/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dpp_drho, err);
        }
    }

}

void
test_pp_T_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    std::println("");
    std::println("testing ∂p⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dpp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<RealT>([&] (RealT _T) -> RealT
                {
                    auto es_eps = eos.pe_state(rho, _T, Ye);
                    return es_eps.p_pos;
                }, T, dT);

            RealT err{};
            if (es.dpp_dT == 0.0L) {
                RealT scale = es.p_pos / T;
                err = std::abs(es.dpp_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dpp_dT - deriv) / std::abs(es.dpp_dT);
            }
            std::println("ρ = {:10} T = {:10}, ∂p⁺/∂T = {:15.8g}, error = {:15.5g}",
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

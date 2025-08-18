#include <print>

#include <array>
#include <cmath>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"
#include "util.H"

constexpr std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
constexpr std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

// ∂e/∂ρ

void
test_ee_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂e⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.e_e;
                }, rho, drho);

            real_t err = std::abs(es.dee_drho - deriv) / std::abs(es.dee_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂e⁻/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dee_drho, err);
        }
    }
}


void
test_ep_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂e⁺/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dep_drho == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.e_pos;
                }, rho, drho);

            real_t err = std::abs(es.dep_drho - deriv) / std::abs(es.dep_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂e⁺/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dep_drho, err);
        }
    }
}


// ∂e/∂T

void
test_ee_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂e⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.e_e;
                }, T, dT);

            real_t err{};
            if (es.dee_dT == 0.0_rt) {
                const real_t scale = es.e_e / T;
                err = std::abs(es.dee_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dee_dT - deriv) / std::abs(es.dee_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂e⁻/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dee_dT, err);
        }
    }
}


void
test_ep_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂e⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dep_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.e_pos;
                }, T, dT);

            real_t err{};
            if (es.dep_dT == 0.0_rt) {
                const real_t scale = es.e_pos / T;
                err = std::abs(es.dep_dT - deriv / scale) ;
            } else {
                err = std::abs(es.dep_dT - deriv) / std::abs(es.dep_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂e⁺/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dep_dT, err);
        }
    }
}


// ∂²e/∂ρ²

void
test_ee_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²e⁻/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.e_e;
                }, rho, drho);

            real_t err = std::abs(es.d2ee_drho2 - deriv) / std::abs(es.d2ee_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²e⁻/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ee_drho2, err);
        }
    }
}


void
test_ep_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²e⁺/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2ep_drho2 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.e_pos;
                }, rho, drho);

            real_t err = std::abs(es.d2ep_drho2 - deriv) / std::abs(es.d2ep_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²e⁺/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ep_drho2, err);
        }
    }

}


// ∂²e/∂T²

void
test_ee_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²e⁻/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.e_e;
                }, T, dT);

            real_t err{};
            if (es.d2ee_dT2 == 0.0_rt) {
                const real_t scale = es.e_e / T / T;
                err = std::abs(es.d2ee_dT2 - deriv / scale) ;
            } else {
                err = std::abs(es.d2ee_dT2 - deriv) / std::abs(es.d2ee_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²e⁻/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ee_dT2, err);
        }
    }
}


void
test_ep_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²e⁺/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2ep_dT2 == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.e_pos;
                }, T, dT);

            real_t err{};
            if (es.d2ep_dT2 == 0.0_rt) {
                const real_t scale = es.e_pos / T / T;
                err = std::abs(es.d2ep_dT2 - deriv / scale) ;
            } else {
                err = std::abs(es.d2ep_dT2 - deriv) / std::abs(es.d2ep_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²e⁺/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ep_dT2, err);
        }
    }
}


// ∂²e/∂ρ∂T

void
test_ee_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²e⁻/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dee_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2ee_drhodT == 0.0_rt) {
                const real_t scale = es.e_e / rho / T;
                err = std::abs(es.d2ee_drhodT - deriv / scale) ;
            } else {
                err = std::abs(es.d2ee_drhodT - deriv) / std::abs(es.d2ee_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²e⁻/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ee_drhodT, err);
        }
    }
}


void
test_ep_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    constexpr real_t Ye{0.5_rt};

    constexpr real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²e⁺/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2ep_drhodT == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dep_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2ep_drhodT == 0.0_rt) {
                const real_t scale = es.e_pos / rho / T;
                err = std::abs(es.d2ep_drhodT - deriv / scale) ;
            } else {
                err = std::abs(es.d2ep_drhodT - deriv) / std::abs(es.d2ep_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²e⁺/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ep_drhodT, err);
        }
    }
}


int main() {

    test_ee_rho_derivs();
    test_ep_rho_derivs();

    test_ee_T_derivs();
    test_ep_T_derivs();

    test_ee_rho2_derivs();
    test_ep_rho2_derivs();

    test_ee_T2_derivs();
    test_ep_T2_derivs();

    test_ee_rhoT_derivs();
    test_ep_rhoT_derivs();

}

#include <print>

#include <array>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"

const std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};


// ∂p/∂ρ

void
test_pe_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂p⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_e;
                }, rho, drho);

            real_t err = mp::abs(es.dpe_drho - deriv) / mp::abs(es.dpe_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂p⁻/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dpe_drho, err);
        }
    }
}


void
test_pp_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂p⁺/∂ρ via differencing");

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

            real_t err = mp::abs(es.dpp_drho - deriv) / mp::abs(es.dpp_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂p⁺/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dpp_drho, err);
        }
    }

}


// ∂p/∂T

void
test_pe_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂p⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.p_e;
                }, T, dT);

            real_t err{};
            if (es.dpe_dT == 0.0_rt) {
                const real_t scale = es.p_e / T;
                err = mp::abs(es.dpe_dT - deriv / scale) ;
            } else {
                err = mp::abs(es.dpe_dT - deriv) / mp::abs(es.dpe_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂p⁻/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dpe_dT, err);
        }
    }
}


void
test_pp_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂p⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dpp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.p_pos;
                }, T, dT);

            real_t err{};
            if (es.dpp_dT == 0.0_rt) {
                const real_t scale = es.p_pos / T;
                err = mp::abs(es.dpp_dT - deriv / scale) ;
            } else {
                err = mp::abs(es.dpp_dT - deriv) / mp::abs(es.dpp_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂p⁺/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dpp_dT, err);
        }
    }
}


// ∂²p/∂ρ²

void
test_pe_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²p⁻/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_e;
                }, rho, drho);

            real_t err = mp::abs(es.d2pe_drho2 - deriv) / mp::abs(es.d2pe_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²p⁻/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2pe_drho2, err);
        }
    }
}


void
test_pp_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²p⁺/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2pp_drho2 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.p_pos;
                }, rho, drho);

            real_t err = mp::abs(es.d2pp_drho2 - deriv) / mp::abs(es.d2pp_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²p⁺/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2pp_drho2, err);
        }
    }

}


// ∂²p/∂T²

void
test_pe_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²p⁻/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.p_e;
                }, T, dT);

            real_t err{};
            if (es.d2pe_dT2 == 0.0_rt) {
                const real_t scale = es.p_e / T / T;
                err = mp::abs(es.d2pe_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2pe_dT2 - deriv) / mp::abs(es.d2pe_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²p⁻/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2pe_dT2, err);
        }
    }
}


void
test_pp_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²p⁺/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2pp_dT2 == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.p_pos;
                }, T, dT);

            real_t err{};
            if (es.d2pp_dT2 == 0.0_rt) {
                const real_t scale = es.p_pos / T / T;
                err = mp::abs(es.d2pp_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2pp_dT2 - deriv) / mp::abs(es.d2pp_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²p⁺/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2pp_dT2, err);
        }
    }
}


// ∂²p/∂ρ∂T

void
test_pe_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²p⁻/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dpe_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2pe_drhodT == 0.0_rt) {
                const real_t scale = es.p_e / rho / T;
                err = mp::abs(es.d2pe_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2pe_drhodT - deriv) / mp::abs(es.d2pe_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²p⁻/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2pe_drhodT, err);
        }
    }
}


void
test_pp_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²p⁺/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2pp_drhodT == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dpp_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2pp_drhodT == 0.0_rt) {
                const real_t scale = es.p_pos / rho / T;
                err = mp::abs(es.d2pp_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2pp_drhodT - deriv) / mp::abs(es.d2pp_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²p⁺/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2pp_drhodT, err);
        }
    }
}


auto main() -> int
{

    test_pe_rho_derivs();
    test_pp_rho_derivs();

    test_pe_T_derivs();
    test_pp_T_derivs();

    test_pe_rho2_derivs();
    test_pp_rho2_derivs();

    test_pe_T2_derivs();
    test_pp_T2_derivs();

    test_pe_rhoT_derivs();
    test_pp_rhoT_derivs();

}

#include <print>

#include <array>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"

const std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

// first derivatives

// ∂n/∂ρ

void
test_ne_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂n⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_e;
                }, rho, drho);

            real_t err = mp::abs(es.dne_drho - deriv) / mp::abs(es.dne_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂n⁻/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dne_drho, err);
        }
    }
}


void
test_np_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂n⁺/∂ρ via differencing");

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

            real_t err = mp::abs(es.dnp_drho - deriv) / mp::abs(es.dnp_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂n⁺/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dnp_drho, err);
        }
    }

}


// ∂n/∂T

void
test_ne_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂n⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.n_e;
                }, T, dT);

            real_t err{};
            if (es.dne_dT == 0.0_rt) {
                const real_t scale = es.n_e / T;
                err = mp::abs(es.dne_dT - deriv / scale) ;
            } else {
                err = mp::abs(es.dne_dT - deriv) / mp::abs(es.dne_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂n⁻/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dne_dT, err);
        }
    }
}


void
test_np_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂n⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dnp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.n_pos;
                }, T, dT);

            real_t err{};
            if (es.dnp_dT == 0.0_rt) {
                const real_t scale = es.n_pos / T;
                err = mp::abs(es.dnp_dT - deriv / scale) ;
            } else {
                err = mp::abs(es.dnp_dT - deriv) / mp::abs(es.dnp_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂n⁺/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dnp_dT, err);
        }
    }
}


// second derivatives

// ∂²n/∂ρ²

void
test_ne_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n⁻/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_e;
                }, rho, drho);

            real_t err{};
            if (es.d2ne_drho2 == 0) {
                err = mp::abs(es.d2ne_drho2 - deriv);
            } else {
                err = mp::abs(es.d2ne_drho2 - deriv) / mp::abs(es.d2ne_drho2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n⁻/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ne_drho2, err);
        }
    }
}


void
test_np_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n⁺/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2np_drho2 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_pos;
                }, rho, drho);

            real_t err{};
            if (es.d2np_drho2 == 0) {
                err = mp::abs(es.d2np_drho2 - deriv);
            } else {
                err = mp::abs(es.d2np_drho2 - deriv) / mp::abs(es.d2np_drho2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n⁺/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2np_drho2, err);
        }
    }

}


// ∂²n/∂T²

void
test_ne_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n⁻/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.n_e;
                }, T, dT);

            real_t err{};
            if (es.d2ne_dT2 == 0.0_rt) {
                const real_t scale = es.n_e / T / T;
                err = mp::abs(es.d2ne_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2ne_dT2 - deriv) / mp::abs(es.d2ne_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n⁻/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ne_dT2, err);
        }
    }
}


void
test_np_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n⁺/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2np_dT2 == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.n_pos;
                }, T, dT);

            real_t err{};
            if (es.d2np_dT2 == 0.0_rt) {
                const real_t scale = es.n_pos / T / T;
                err = mp::abs(es.d2np_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2np_dT2 - deriv) / mp::abs(es.d2np_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n⁺/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2np_dT2, err);
        }
    }
}


// ∂²n/∂ρ∂T

void
test_ne_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n⁻/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dne_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2ne_drhodT == 0.0_rt) {
                const real_t scale = es.n_e / rho / T;
                err = mp::abs(es.d2ne_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2ne_drhodT - deriv) / mp::abs(es.d2ne_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n⁻/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2ne_drhodT, err);
        }
    }
}


void
test_np_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n⁺/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2np_drhodT == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dnp_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2np_drhodT == 0.0_rt) {
                const real_t scale = es.n_pos / rho / T;
                err = mp::abs(es.d2np_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2np_drhodT - deriv) / mp::abs(es.d2np_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n⁺/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2np_drhodT, err);
        }
    }
}


// third derivatives

// ∂³n/∂ρ³

void
test_ne_rho3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁻/∂ρ³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dne_drho;
                }, rho, drho);

            real_t err{};
            if (es.d3ne_drho3 == 0) {
                err = mp::abs(es.d3ne_drho3 - deriv);
            } else {
                err = mp::abs(es.d3ne_drho3 - deriv) / mp::abs(es.d3ne_drho3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁻/∂ρ³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3ne_drho3, err);
        }
    }
}


void
test_np_rho3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁺/∂ρ³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3np_drho3 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dnp_drho;
                }, rho, drho);

            real_t err{};
            if (es.d3np_drho3 == 0) {
                err = mp::abs(es.d3np_drho3 - deriv);
            } else {
                err = mp::abs(es.d3np_drho3 - deriv) / mp::abs(es.d3np_drho3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁺/∂ρ³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3np_drho3, err);
        }
    }

}


// ∂³n/∂ρ²∂T

void
test_ne_rho2T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁻/∂ρ²∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dne_dT;
                }, rho, drho);

            real_t err{};
            if (es.d3ne_drho2dT == 0) {
                err = mp::abs(es.d3ne_drho2dT - deriv);
            } else {
                err = mp::abs(es.d3ne_drho2dT - deriv) / mp::abs(es.d3ne_drho2dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁻/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3ne_drho2dT, err);
        }
    }
}

void
test_np_rho2T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁺/∂ρ²∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3np_drho2dT == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dnp_dT;
                }, rho, drho);

            real_t err{};
            if (es.d3np_drho3 == 0) {
                err = mp::abs(es.d3np_drho2dT - deriv);
            } else {
                err = mp::abs(es.d3np_drho2dT - deriv) / mp::abs(es.d3np_drho2dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁺/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3np_drho2dT, err);
        }
    }

}



auto main() -> int
{

#if 0
    test_ne_rho_derivs();
    test_np_rho_derivs();

    test_ne_T_derivs();
    test_np_T_derivs();

    test_ne_rho2_derivs();
    test_np_rho2_derivs();

    test_ne_T2_derivs();
    test_np_T2_derivs();

    test_ne_rhoT_derivs();
    test_np_rhoT_derivs();
#endif

    test_ne_rho3_derivs();
    test_np_rho3_derivs();

    test_ne_rho2T_derivs();
    test_np_rho2T_derivs();

}

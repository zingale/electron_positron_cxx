#include <print>

#include "real_type.H"
#include "electron_positron.H"
#include "eos_types.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"
#include "derivative_comparison.H"


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

void
test_n_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n;
                }, rho, drho);

            real_t err{};
            if (es.d2n_drho2 == 0) {
                err = mp::abs(es.d2n_drho2 - deriv);
            } else {
                err = mp::abs(es.d2n_drho2 - deriv) / mp::abs(es.d2n_drho2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2n_drho2, err);
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

void
test_n_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.n;
                }, T, dT);

            real_t err{};
            if (es.d2n_dT2 == 0.0_rt) {
                const real_t scale = es.n / T / T;
                err = mp::abs(es.d2n_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2n_dT2 - deriv) / mp::abs(es.d2n_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2n_dT2, err);
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

void
test_n_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²n/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dn_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2n_drhodT == 0.0_rt) {
                const real_t scale = es.n / rho / T;
                err = mp::abs(es.d2n_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2n_drhodT - deriv) / mp::abs(es.d2n_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²n/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2n_drhodT, err);
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

void
test_n_rho3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n/∂ρ³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dn_drho;
                }, rho, drho);

            real_t err{};
            if (es.d3n_drho3 == 0) {
                err = mp::abs(es.d3n_drho3 - deriv);
            } else {
                err = mp::abs(es.d3n_drho3 - deriv) / mp::abs(es.d3n_drho3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n/∂ρ³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3n_drho3, err);
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
            if (es.d3np_drho2dT == 0) {
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

void
test_n_rho2T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n/∂ρ²∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dn_dT;
                }, rho, drho);

            real_t err{};
            if (es.d3n_drho2dT == 0) {
                err = mp::abs(es.d3n_drho2dT - deriv);
            } else {
                err = mp::abs(es.d3n_drho2dT - deriv) / mp::abs(es.d3n_drho2dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3n_drho2dT, err);
        }
    }

}

// ∂³n/∂ρ∂T²

void
test_ne_rhoT2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁻/∂ρ∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.d2ne_dT2;
                }, rho, drho);

            real_t err{};
            if (es.d3ne_drhodT2 == 0) {
                err = mp::abs(es.d3ne_drhodT2 - deriv);
            } else {
                err = mp::abs(es.d3ne_drhodT2 - deriv) / mp::abs(es.d3ne_drhodT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁻/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3ne_drhodT2, err);
        }
    }
}

void
test_np_rhoT2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁺/∂ρ∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3np_drhodT2 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.d2np_dT2;
                }, rho, drho);

            real_t err{};
            if (es.d3np_drhodT2 == 0) {
                err = mp::abs(es.d3np_drhodT2 - deriv);
            } else {
                err = mp::abs(es.d3np_drhodT2 - deriv) / mp::abs(es.d3np_drhodT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁺/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3np_drhodT2, err);
        }
    }
}

void
test_n_rhoT2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n/∂ρ∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.d2n_dT2;
                }, rho, drho);

            real_t err{};
            if (es.d3n_drhodT2 == 0) {
                err = mp::abs(es.d3n_drhodT2 - deriv);
            } else {
                err = mp::abs(es.d3n_drhodT2 - deriv) / mp::abs(es.d3n_drhodT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3n_drhodT2, err);
        }
    }
}

// ∂³n/∂T³

void
test_ne_T3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁻/∂T³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.dne_dT;
                }, T, dT);

            real_t err{};
            if (es.d3ne_dT3 == 0.0_rt) {
                const real_t scale = es.n_e / T / T / T;
                err = mp::abs(es.d3ne_dT3 - deriv / scale) ;
            } else {
                err = mp::abs(es.d3ne_dT3 - deriv) / mp::abs(es.d3ne_dT3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁻/∂T³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3ne_dT3, err);
        }
    }
}

void
test_np_T3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n⁺/∂T³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3np_dT3 == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.dnp_dT;
                }, T, dT);

            real_t err{};
            if (es.d3np_dT3 == 0.0_rt) {
                const real_t scale = es.n_pos / T / T / T;
                err = mp::abs(es.d3np_dT3 - deriv / scale) ;
            } else {
                err = mp::abs(es.d3np_dT3 - deriv) / mp::abs(es.d3np_dT3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n⁺/∂T³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3np_dT3, err);
        }
    }
}

void
test_n_T3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³n/∂T³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.dn_dT;
                }, T, dT);

            real_t err{};
            if (es.d3n_dT3 == 0.0_rt) {
                const real_t scale = es.n / T / T / T;
                err = mp::abs(es.d3n_dT3 - deriv / scale) ;
            } else {
                err = mp::abs(es.d3n_dT3 - deriv) / mp::abs(es.d3n_dT3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³n/∂T³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3n_dT3, err);
        }
    }
}


auto main() -> int
{

    // test ∂n/∂ρ
    test_rho_deriv(&EOSState<real_t>::dne_drho, &EOSState<real_t>::n_e, "n⁻");
    test_rho_deriv(&EOSState<real_t>::dnp_drho, &EOSState<real_t>::n_pos, "n⁺");
    test_rho_deriv(&EOSState<real_t>::dn_drho, &EOSState<real_t>::n, "n");

    test_T_deriv(&EOSState<real_t>::dne_dT, &EOSState<real_t>::n_e, "n⁻");
    test_T_deriv(&EOSState<real_t>::dnp_dT, &EOSState<real_t>::n_pos, "n⁺");
    test_T_deriv(&EOSState<real_t>::dn_dT, &EOSState<real_t>::n, "n");

    test_ne_rho2_derivs();
    test_np_rho2_derivs();
    test_n_rho2_derivs();

    test_ne_T2_derivs();
    test_np_T2_derivs();
    test_n_T2_derivs();

    test_ne_rhoT_derivs();
    test_np_rhoT_derivs();
    test_n_rhoT_derivs();

    test_ne_rho3_derivs();
    test_np_rho3_derivs();
    test_n_rho3_derivs();

    test_ne_rho2T_derivs();
    test_np_rho2T_derivs();
    test_n_rho2T_derivs();

    test_ne_rhoT2_derivs();
    test_np_rhoT2_derivs();
    test_n_rhoT2_derivs();

    test_ne_T3_derivs();
    test_np_T3_derivs();
    test_n_T3_derivs();

}

#include <print>

#include <array>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"

const std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};


// ∂s/∂ρ

void
test_se_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂s⁻/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.s_e;
                }, rho, drho);

            real_t err = mp::abs(es.dse_drho - deriv) / mp::abs(es.dse_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂s⁻/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dse_drho, err);
        }
    }
}


void
test_sp_rho_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂s⁺/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dsp_drho == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.s_pos;
                }, rho, drho);

            real_t err = mp::abs(es.dsp_drho - deriv) / mp::abs(es.dsp_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂s⁺/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dsp_drho, err);
        }
    }
}


// ∂s/∂T

void
test_se_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂s⁻/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.s_e;
                }, T, dT);

            real_t err{};
            if (es.dse_dT == 0.0_rt) {
                const real_t scale = es.s_e / T;
                err = mp::abs(es.dse_dT - deriv / scale) ;
            } else {
                err = mp::abs(es.dse_dT - deriv) / mp::abs(es.dse_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂s⁻/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dse_dT, err);
        }
    }
}


void
test_sp_T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂s⁺/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.dsp_dT == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.s_pos;
                }, T, dT);

            real_t err{};
            if (es.dsp_dT == 0.0_rt) {
                const real_t scale = es.s_pos / T;
                err = mp::abs(es.dsp_dT - deriv / scale) ;
            } else {
                err = mp::abs(es.dsp_dT - deriv) / mp::abs(es.dsp_dT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂s⁺/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.dsp_dT, err);
        }
    }
}


// second derivatives

// ∂²s/∂ρ²

void
test_se_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²s⁻/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.s_e;
                }, rho, drho);

            real_t err = mp::abs(es.d2se_drho2 - deriv) / mp::abs(es.d2se_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²s⁻/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2se_drho2, err);
        }
    }
}


void
test_sp_rho2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²s⁺/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2sp_drho2 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.s_pos;
                }, rho, drho);

            real_t err = mp::abs(es.d2sp_drho2 - deriv) / mp::abs(es.d2sp_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²s⁺/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2sp_drho2, err);
        }
    }

}


// ∂²s/∂T²

void
test_se_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²s⁻/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.s_e;
                }, T, dT);

            real_t err{};
            if (es.d2se_dT2 == 0.0_rt) {
                const real_t scale = es.s_e / T / T;
                err = mp::abs(es.d2se_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2se_dT2 - deriv) / mp::abs(es.d2se_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²s⁻/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2se_dT2, err);
        }
    }
}


void
test_sp_T2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²s⁺/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2sp_dT2 == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.s_pos;
                }, T, dT);

            real_t err{};
            if (es.d2sp_dT2 == 0.0_rt) {
                const real_t scale = es.s_pos / T / T;
                err = mp::abs(es.d2sp_dT2 - deriv / scale) ;
            } else {
                err = mp::abs(es.d2sp_dT2 - deriv) / mp::abs(es.d2sp_dT2);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²s⁺/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2sp_dT2, err);
        }
    }
}


// ∂²s/∂ρ∂T

void
test_se_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²s⁻/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dse_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2se_drhodT == 0.0_rt) {
                const real_t scale = es.s_e / rho / T;
                err = mp::abs(es.d2se_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2se_drhodT - deriv) / mp::abs(es.d2se_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²s⁻/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2se_drhodT, err);
        }
    }
}


void
test_sp_rhoT_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²s⁺/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d2sp_drhodT == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return es_eps.dsp_dT;
                }, rho, drho);

            real_t err{};
            if (es.d2sp_drhodT == 0.0_rt) {
                const real_t scale = es.s_pos / rho / T;
                err = mp::abs(es.d2sp_drhodT - deriv / scale) ;
            } else {
                err = mp::abs(es.d2sp_drhodT - deriv) / mp::abs(es.d2sp_drhodT);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²s⁺/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d2sp_drhodT, err);
        }
    }
}


// third derivatives

// ∂³s/∂ρ³

void
test_se_rho3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁻/∂ρ³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dse_drho;
                }, rho, drho);

            real_t err = mp::abs(es.d3se_drho3 - deriv) / mp::abs(es.d3se_drho3);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁻/∂ρ³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3se_drho3, err);
        }
    }
}


void
test_sp_rho3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁺/∂ρ³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3sp_drho3 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dsp_drho;
                }, rho, drho);

            real_t err = mp::abs(es.d3sp_drho3 - deriv) / mp::abs(es.d3sp_drho3);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁺/∂ρ³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3sp_drho3, err);
        }
    }

}

// ∂³s/∂ρ²∂T

void
test_se_rho2T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁻/∂ρ²∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dse_dT;
                }, rho, drho);

            real_t err = mp::abs(es.d3se_drho2dT - deriv) / mp::abs(es.d3se_drho2dT);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁻/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3se_drho2dT, err);
        }
    }
}

void
test_sp_rho2T_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁺/∂ρ²∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3sp_drho2dT == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.dsp_dT;
                }, rho, drho);

            real_t err = mp::abs(es.d3sp_drho2dT - deriv) / mp::abs(es.d3sp_drho2dT);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁺/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3sp_drho2dT, err);
        }
    }

}

// ∂³s/∂ρ∂T²

void
test_se_rhoT2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁻/∂ρ∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.d2se_dT2;
                }, rho, drho);

            real_t err = mp::abs(es.d3se_drhodT2 - deriv) / mp::abs(es.d3se_drhodT2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁻/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3se_drhodT2, err);
        }
    }
}

void
test_sp_rhoT2_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁺/∂ρ∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3sp_drhodT2 == 0.0) {
                continue;
            }
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.d2sp_dT2;
                }, rho, drho);

            real_t err = mp::abs(es.d3sp_drhodT2 - deriv) / mp::abs(es.d3sp_drhodT2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁺/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3sp_drhodT2, err);
        }
    }

}


// ∂³s/∂T³

void
test_se_T3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁻/∂T³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.dse_dT;
                }, T, dT);

            real_t err{};
            if (es.d3se_dT3 == 0.0_rt) {
                const real_t scale = es.s_e / T / T / T;
                err = mp::abs(es.d3se_dT3 - deriv / scale) ;
            } else {
                err = mp::abs(es.d3se_dT3 - deriv) / mp::abs(es.d3se_dT3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁻/∂T³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3se_dT3, err);
        }
    }
}


void
test_sp_T3_derivs() {

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³s⁺/∂T³ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto es = eos.pe_state(rho, T, Ye);
            if (es.n_pos == 0.0 && es.d3sp_dT3 == 0.0) {
                continue;
            }
            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return es_eps.dsp_dT;
                }, T, dT);

            real_t err{};
            if (es.d3sp_dT3 == 0.0_rt) {
                const real_t scale = es.s_pos / T / T / T;
                err = mp::abs(es.d3sp_dT3 - deriv / scale) ;
            } else {
                err = mp::abs(es.d3sp_dT3 - deriv) / mp::abs(es.d3sp_dT3);
            }
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³s⁺/∂T³ = {:15.8g},  error = {:11.5g}",
                                    rho, T, es.d3sp_dT3, err);
        }
    }
}


auto main() -> int
{

#if 0
    test_se_rho_derivs();
    test_se_T_derivs();

    test_sp_rho_derivs();
    test_sp_T_derivs();

    test_se_rho2_derivs();
    test_sp_rho2_derivs();

    test_se_T2_derivs();
    test_sp_T2_derivs();

    test_se_rhoT_derivs();
    test_sp_rhoT_derivs();
#endif

    test_se_rho3_derivs();
    //test_sp_rho3_derivs();

    test_se_rho2T_derivs();
    //test_sp_rho2T_derivs();

    test_se_rhoT2_derivs();
    //test_sp_rhoT2_derivs();

    test_se_T3_derivs();
    //test_sp_T3_derivs();

}

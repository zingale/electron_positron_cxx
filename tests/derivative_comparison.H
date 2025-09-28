#ifndef DERIVATIVE_HELPERS_H
#define DERIVATIVE_HELPERS_H

#include <array>
#include <vector>
#include <string_view>

#include "real_type.H"
#include "electron_positron.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"


const std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};

struct TestState {
    real_t rho{};
    real_t T{};
    real_t value{};
    real_t error{};
};


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_rho_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                    std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂{}/∂ρ via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto drho{eps * rho};

            auto [deriv, _err] = fd::adaptive_diff<real_t>(
                [&] (real_t _rho) -> real_t {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return std::invoke(diff_member, es_eps);
                }, rho, drho);

            const auto ref = std::invoke(test_member, es);
            const auto err = mp::abs(ref - deriv) / mp::abs(ref);

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂{}/∂ρ = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_T_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                    std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂{}/∂T via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto dT{eps * T};

            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return std::invoke(diff_member, es_eps);
                }, T, dT);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0.0_rt) {
                const real_t scale = std::invoke(diff_member, es) / T;
                err = mp::abs(ref - deriv / scale) ;
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂{}/∂T = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_rho2_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                     std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²{}/∂ρ² via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto drho{eps * rho};

            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return std::invoke(diff_member, es_eps);
                }, rho, drho);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0) {
                err = mp::abs(ref - deriv);
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂²{}/∂ρ² = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_T2_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                   std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²{}/∂T² via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto dT{eps * T};

            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return std::invoke(diff_member, es_eps);
                }, T, dT);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0.0_rt) {
                const real_t scale = std::invoke(diff_member, es) / T / T;
                err = mp::abs(ref - deriv / scale);
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂²{}/∂T² = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_rhoT_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                     std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²{}/∂ρ∂T via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto drho{eps * rho};

            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho_, T, Ye);
                    return std::invoke(diff_member, es_eps);
                }, rho, drho);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0.0_rt) {
                const real_t scale = std::invoke(diff_member, es) / rho / T;
                err = mp::abs(ref - deriv / scale) ;
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂²{}/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_rho3_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                     std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³{}/∂ρ³ via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto drho{eps * rho};

            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return std::invoke(diff_member, es_eps);
                }, rho, drho);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0) {
                err = mp::abs(ref - deriv);
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂³{}/∂ρ³ = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_rho2T_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                      std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³{}/∂ρ²∂T via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto drho{eps * rho};

            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return std::invoke(diff_member, es_eps);
                }, rho, drho);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0) {
                err = mp::abs(ref - deriv);
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂³{}/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_rhoT2_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                      std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³{}/∂ρ∂T² via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto drho{eps * rho};

            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return std::invoke(diff_member, es_eps);
                }, rho, drho);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0) {
                err = mp::abs(ref - deriv);
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂³{}/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


template <typename TestMemberPtr, typename DiffMemberPtr>
void test_T3_deriv(TestMemberPtr test_member, DiffMemberPtr diff_member,
                   std::string_view name) {

    // TestMemberPtr is what we are testing (the reference)
    // DiffMemberPtr is what we are differentiating to compare with the reference

    ElectronPositronEOS<real_t> eos;
    const real_t Ye{0.5_rt};
    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³{}/∂T³ via differencing", name);

    std::vector<TestState> tests(Ts.size() * rhos.size());

    #pragma omp parallel for collapse(2)
    for (int it = 0; it < static_cast<int>(Ts.size()); ++it) {
        for (int ir = 0; ir < static_cast<int>(rhos.size()); ++ir) {
            const auto rho = rhos[ir];
            const auto T = Ts[it];
            const int idx = it * static_cast<int>(rhos.size()) + ir;

            auto es = eos.pe_state(rho, T, Ye);
            if (std::invoke(diff_member, es) == 0.0_rt &&
                std::invoke(test_member, es) == 0.0_rt) {
                continue;
            }
            auto dT{eps * T};

            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto es_eps = eos.pe_state(rho, T_, Ye);
                    return std::invoke(diff_member, es_eps);
                }, T, dT);

            const auto ref = std::invoke(test_member, es);
            real_t err{};
            if (ref == 0.0_rt) {
                const real_t scale = std::invoke(diff_member, es) / T / T / T;
                err = mp::abs(ref - deriv / scale) ;
            } else {
                err = mp::abs(ref - deriv) / mp::abs(ref);
            }

            tests[idx].rho = rho;
            tests[idx].T = T;
            tests[idx].value = ref;
            tests[idx].error = err;
        }
    }

    for (auto const& test : tests) {
        if (test.rho == 0.0_rt) {
            continue;
        }
        util::threshold_println(test.error,
                                "ρ = {:8.3g} T = {:8.3g},  ∂³{}/∂T³ = {:15.8g},  error = {:11.5g}",
                                test.rho, test.T, name, test.value, test.error);
    }
}


#endif

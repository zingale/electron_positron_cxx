#include <print>

#include <array>

#include "real_type.H"
#include "helmholtz.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"

const std::array<real_t, 4> Ts{1.e4_rt, 1.e6_rt, 1.e8_rt, 5.e9_rt};
const std::array<real_t, 5> rhos{1.e-2_rt, 1.e2_rt, 1.e5_rt, 1.e7_rt, 5.e9_rt};


// ∂F/∂ρ

void
test_helm_rho_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂F/∂ρ via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(_rho, T, Ye);
                    return helm_eps.F;
                }, rho, drho);

            real_t err = mp::abs(helm.dF_drho - deriv) / mp::abs(helm.dF_drho);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂F/∂ρ = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.dF_drho, err);
        }
    }
}



// ∂F/∂T

void
test_helm_T_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂F/∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(rho, T_, Ye);
                    return helm_eps.F;
                }, T, dT);

            real_t err = mp::abs(helm.dF_dT - deriv) / mp::abs(helm.dF_dT);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂F/∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.dF_dT, err);
        }
    }
}



// ∂²F/∂ρ²

void
test_helm_rho2_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²F/∂ρ² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t _rho) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(_rho, T, Ye);
                    return helm_eps.F;
                }, rho, drho);

            real_t err = mp::abs(helm.d2F_drho2 - deriv) / mp::abs(helm.d2F_drho2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²F/∂ρ² = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.d2F_drho2, err);
        }
    }
}


// ∂²F/∂T²

void
test_helm_T2_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²F/∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff2<real_t>([&] (real_t T_) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(rho, T_, Ye);
                    return helm_eps.F;
                }, T, dT);

            real_t err = mp::abs(helm.d2F_dT2 - deriv) / mp::abs(helm.d2F_dT2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²F/∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.d2F_dT2, err);
        }
    }
}


// ∂²p/∂ρ∂T

void
test_helm_rhoT_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂²F/∂ρ∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t rho_) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(rho_, T, Ye);
                    return helm_eps.dF_dT;
                }, rho, drho);

            real_t err = mp::abs(helm.d2F_drhodT - deriv) / mp::abs(helm.d2F_drhodT);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂²F/∂ρ∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.d2F_drhodT, err);
        }
    }
}


// ∂³F/∂ρ∂T²

void
test_helm_rhoT2_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³F/∂ρ∂T² via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(_rho, T, Ye);
                    return helm_eps.d2F_dT2;
                }, rho, drho);

            real_t err = mp::abs(helm.d3F_drhodT2 - deriv) / mp::abs(helm.d3F_drhodT2);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³F/∂ρ∂T² = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.d3F_drhodT2, err);
        }
    }
}


// ∂³F/∂ρ²∂T²

void
test_helm_rho2T_derivs() {

    const real_t Ye{0.5_rt};

    const real_t eps{0.01_rt};

    std::println("");
    util::green_println("testing ∂³F/∂ρ²∂T via differencing");

    for (auto T : Ts) {
        for (auto rho : rhos) {
            auto [helm, eos] = get_helmholtz_terms(rho, T, Ye);

            auto dT{eps * T};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t T_) -> real_t
                {
                    auto [helm_eps, eos_eps] = get_helmholtz_terms(rho, T_, Ye);
                    return helm_eps.d2F_drho2;
                }, T, dT);

            real_t err = mp::abs(helm.d3F_drho2dT - deriv) / mp::abs(helm.d3F_drho2dT);
            util::threshold_println(err,
                                    "ρ = {:8.3g} T = {:8.3g},  ∂³F/∂ρ²∂T = {:15.8g},  error = {:11.5g}",
                                    rho, T, helm.d3F_drho2dT, err);
        }
    }
}



auto main() -> int
{

    test_helm_rho_derivs();
    test_helm_T_derivs();

    test_helm_rho2_derivs();
    test_helm_T2_derivs();
    test_helm_rhoT_derivs();

    test_helm_rhoT2_derivs();
    test_helm_rho2T_derivs();

}

#include <cassert>
#include <cmath>
#include <print>

#include "real_type.H"
#include "fermi_integrals.H"
#include "difference_utils.H"
#include "electron_positron.H"
#include "util.H"


int main() {

    // value below which we assume that the positron contribution is zero
    constexpr real_t pos_thresh{1.e-500_rt};

    {

        const real_t h = 0.05_rt;

        util::green_println("testing ∂n⁻/∂η");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check ∂n⁻/∂η
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                        {
                            auto ne = n_e_constraint<real_t>(_eta, beta);
                            return ne;
                        }, eta, _h);

                real_t rel_err = std::abs(dn_e.deta - diff) / std::abs(dn_e.deta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁻/∂η = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_e.deta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂n⁺/∂η");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(dn_pos.deta) < pos_thresh) {
                    continue;
                }

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check ∂n⁺/∂η
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto np = n_p_constraint<real_t>(_eta, beta);
                        return np;
                    }, eta, _h);

                real_t rel_err = std::abs(dn_pos.deta - diff) / std::abs(dn_pos.deta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁺/∂η = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_pos.deta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂n⁻/∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = h * std::abs(beta);

                // check ∂n⁻/∂β
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(eta, _beta);
                        return ne;
                    }, beta, _h);

                real_t rel_err = std::abs(dn_e.dbeta - diff) / std::abs(dn_e.dbeta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁻/∂β = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_e.dbeta, rel_err);
            }
        }
    }

    {
        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂n⁺/∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(dn_pos.dbeta) < pos_thresh) {
                    continue;
                }

                const real_t _h = h * std::abs(beta);

                // check ∂n⁺/∂β
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto np = n_p_constraint<real_t>(eta, _beta);
                        return np;
                    }, beta, _h);

                real_t rel_err = std::abs(dn_pos.dbeta - diff) / std::abs(dn_pos.dbeta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁺/∂β = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_pos.dbeta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁻/∂η²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁻/∂η²
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(_eta, beta);
                        return ne;
                    }, eta, _h);

                real_t rel_err = std::abs(dn_e.deta2 - diff) / std::abs(dn_e.deta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁻/∂η² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_e.deta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁺/∂η²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(dn_pos.deta2) < pos_thresh) {
                    continue;
                }

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁺/∂η²
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto np = n_p_constraint<real_t>(_eta, beta);
                        return np;
                    }, eta, _h);

                real_t rel_err = std::abs(dn_pos.deta2 - diff) / std::abs(dn_pos.deta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁺/∂η² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_pos.deta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁻/∂β²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = h * std::abs(beta);

                // check ∂²n⁻/∂β²
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(eta, _beta);
                        return ne;
                    }, beta, _h);

                real_t rel_err = std::abs(dn_e.dbeta2 - diff) / std::abs(dn_e.dbeta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁻/∂β² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_e.dbeta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁺/∂β²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs
                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(dn_pos.dbeta2) < pos_thresh) {
                    continue;
                }

                const real_t _h = h * std::abs(beta);

                // check ∂²n⁺/∂β²
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto np = n_p_constraint<real_t>(eta, _beta);
                        return np;
                    }, beta, _h);

                real_t rel_err = std::abs(dn_pos.dbeta2 - diff) / std::abs(dn_pos.dbeta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁺/∂β² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, dn_pos.dbeta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁻/∂η∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs

                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h1 = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁻/∂η∂β by differencing ∂n⁻/∂β
                auto [diff_eta, err1] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        const real_t eta_tilde_tmp = -_eta - 2.0_rt / beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, _eta, beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, _eta, beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, beta);
                        f32_pos_tmp.evaluate(1);

                        const auto [dn_e_tmp, dn_pos_tmp] = get_n_derivs<real_t>(beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return dn_e_tmp.dbeta;
                    }, eta, _h1);

                real_t err_eta = std::abs(dn_e.detadbeta - diff_eta) / std::abs(dn_e.detadbeta);

                // check ∂²n⁻/∂η∂β by differencing ∂²n⁻/∂η∂β

                const real_t _h2 = h * std::abs(beta);

                // check ∂²n⁻/∂η∂β by differencing ∂n⁻/∂η
                auto [diff_beta, err2] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        const real_t eta_tilde_tmp = -eta - 2.0_rt / _beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, eta, _beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, eta, _beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, _beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, _beta);
                        f32_pos_tmp.evaluate(1);

                        const auto [dn_e_tmp, dn_pos_tmp] = get_n_derivs<real_t>(_beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return dn_e_tmp.deta;
                    }, beta, _h2);

                real_t err_beta = std::abs(dn_e.detadbeta - diff_beta) / std::abs(dn_e.detadbeta);

                util::threshold_println(err_eta,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁻/∂η∂β = {:12.5g},  D_η(∂n⁻/∂β) error = {:12.5g}  D_β(∂n⁻/∂η) error = {:12.5g}",
                                        eta, beta, dn_e.detadbeta, err_eta, err_beta);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁺/∂η∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                const real_t eta_tilde = -eta - 2.0_rt / beta;

                // construct the needed Fermi integrals
                FermiIntegral<real_t> f12(0.5_rt, eta, beta);
                f12.evaluate(2);

                FermiIntegral<real_t> f32(1.5_rt, eta, beta);
                f32.evaluate(2);

                FermiIntegral<real_t> f12_pos(0.5_rt, eta_tilde, beta);
                f12_pos.evaluate(2);

                FermiIntegral<real_t> f32_pos(1.5_rt, eta_tilde, beta);
                f32_pos.evaluate(2);

                // get the derivs

                const auto [dn_e, dn_pos] = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(dn_pos.detadbeta) < pos_thresh) {
                    continue;
                }

                const real_t _h1 = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁺/∂η∂β by differencing ∂n⁺/∂β
                auto [diff_eta, err1] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        const real_t eta_tilde_tmp = -_eta - 2.0_rt / beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, _eta, beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, _eta, beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, beta);
                        f32_pos_tmp.evaluate(1);

                        const auto [dn_e_tmp, dn_pos_tmp] = get_n_derivs<real_t>(beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return dn_pos_tmp.dbeta;
                    }, eta, _h1);

                real_t err_eta = std::abs(dn_pos.detadbeta - diff_eta) / std::abs(dn_pos.detadbeta);

                // check ∂²n⁺/∂η∂β by differencing ∂²n⁺/∂η∂β

                const real_t _h2 = h * std::abs(beta);

                // check ∂²n⁻/∂η∂β by differencing ∂n⁻/∂η
                auto [diff_beta, err2] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        const real_t eta_tilde_tmp = -eta - 2.0_rt / _beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, eta, _beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, eta, _beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, _beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, _beta);
                        f32_pos_tmp.evaluate(1);

                        const auto [dn_e_tmp, dn_pos_tmp] = get_n_derivs<real_t>(_beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return dn_pos_tmp.deta;
                    }, beta, _h2);

                real_t err_beta = std::abs(dn_pos.detadbeta - diff_beta) / std::abs(dn_pos.detadbeta);

                util::threshold_println(err_eta,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁺/∂η∂β = {:12.5g},  D_η(∂n⁺/∂β) error = {:12.5g}  D_β(∂n⁺/∂η) error = {:12.5g}",
                                        eta, beta, dn_pos.detadbeta, err_eta, err_beta);
            }
        }
    }

}


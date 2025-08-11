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

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check ∂n⁻/∂η
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(_eta, beta);
                        return ne;
                    }, eta, _h);

                real_t rel_err = std::abs(nd.dne_deta - diff) / std::abs(nd.dne_deta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁻/∂η = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.dne_deta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂n⁺/∂η");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(nd.dnp_deta) < pos_thresh) {
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

                real_t rel_err = std::abs(nd.dnp_deta - diff) / std::abs(nd.dnp_deta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁺/∂η = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.dnp_deta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂n⁻/∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = h * std::abs(beta);

                // check ∂n⁻/∂β
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(eta, _beta);
                        return ne;
                    }, beta, _h);

                real_t rel_err = std::abs(nd.dne_dbeta - diff) / std::abs(nd.dne_dbeta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁻/∂β = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.dne_dbeta, rel_err);
            }
        }
    }

    {
        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂n⁺/∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(nd.dnp_dbeta) < pos_thresh) {
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

                real_t rel_err = std::abs(nd.dnp_dbeta - diff) / std::abs(nd.dnp_dbeta);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂n⁺/∂β = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.dnp_dbeta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁻/∂η²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁻/∂η²
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(_eta, beta);
                        return ne;
                    }, eta, _h);

                real_t rel_err = std::abs(nd.d2ne_deta2 - diff) / std::abs(nd.d2ne_deta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁻/∂η² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.d2ne_deta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁺/∂η²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(nd.d2np_deta2) < pos_thresh) {
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

                real_t rel_err = std::abs(nd.d2np_deta2 - diff) / std::abs(nd.d2np_deta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁺/∂η² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.d2np_deta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁻/∂β²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h = h * std::abs(beta);

                // check ∂²n⁻/∂β²
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(eta, _beta);
                        return ne;
                    }, beta, _h);

                real_t rel_err = std::abs(nd.d2ne_dbeta2 - diff) / std::abs(nd.d2ne_dbeta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁻/∂β² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.d2ne_dbeta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁺/∂β²");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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
                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(nd.d2np_dbeta2) < pos_thresh) {
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

                real_t rel_err = std::abs(nd.d2np_dbeta2 - diff) / std::abs(nd.d2np_dbeta2);

                util::threshold_println(rel_err,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁺/∂β² = {:12.5g},  error = {:12.5g}",
                                        eta, beta, nd.d2np_dbeta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁻/∂η∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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

                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                const real_t _h1 = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁻/∂η∂β by differencing ∂n⁻/∂β
                auto [diff_eta, err1] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        real_t eta_tilde_tmp = -_eta - 2.0_rt / beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, _eta, beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, _eta, beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, beta);
                        f32_pos_tmp.evaluate(1);

                        auto nd_tmp = get_n_derivs<real_t>(beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return nd_tmp.dne_dbeta;
                    }, eta, _h1);

                real_t err_eta = std::abs(nd.d2ne_detadbeta - diff_eta) / std::abs(nd.d2ne_detadbeta);

                // check ∂²n⁻/∂η∂β by differencing ∂²n⁻/∂η∂β

                const real_t _h2 = h * std::abs(beta);

                // check ∂²n⁻/∂η∂β by differencing ∂n⁻/∂η
                auto [diff_beta, err2] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        real_t eta_tilde_tmp = -eta - 2.0_rt / _beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, eta, _beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, eta, _beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, _beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, _beta);
                        f32_pos_tmp.evaluate(1);

                        auto nd_tmp = get_n_derivs<real_t>(_beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return nd_tmp.dne_deta;
                    }, beta, _h2);

                real_t err_beta = std::abs(nd.d2ne_detadbeta - diff_beta) / std::abs(nd.d2ne_detadbeta);

                util::threshold_println(err_eta,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁻/∂η∂β = {:12.5g},  D_η(∂n⁻/∂β) error = {:12.5g}  D_β(∂n⁻/∂η) error = {:12.5g}",
                                        eta, beta, nd.d2ne_detadbeta, err_eta, err_beta);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing ∂²n⁺/∂η∂β");

        for (auto eta : {-10.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
            for (auto beta : {1.e-6_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                real_t eta_tilde = -eta - 2.0_rt / beta;

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

                auto nd = get_n_derivs<real_t>(beta, f12, f32, f12_pos, f32_pos);

                if (std::abs(nd.d2np_detadbeta) < pos_thresh) {
                    continue;
                }

                const real_t _h1 = (eta == 0) ? h : h * std::abs(eta);

                // check ∂²n⁺/∂η∂β by differencing ∂n⁺/∂β
                auto [diff_eta, err1] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        real_t eta_tilde_tmp = -_eta - 2.0_rt / beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, _eta, beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, _eta, beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, beta);
                        f32_pos_tmp.evaluate(1);

                        auto nd_tmp = get_n_derivs<real_t>(beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return nd_tmp.dnp_dbeta;
                    }, eta, _h1);

                real_t err_eta = std::abs(nd.d2np_detadbeta - diff_eta) / std::abs(nd.d2np_detadbeta);

                // check ∂²n⁺/∂η∂β by differencing ∂²n⁺/∂η∂β

                const real_t _h2 = h * std::abs(beta);

                // check ∂²n⁻/∂η∂β by differencing ∂n⁻/∂η
                auto [diff_beta, err2] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        real_t eta_tilde_tmp = -eta - 2.0_rt / _beta;

                        FermiIntegral<real_t> f12_tmp(0.5_rt, eta, _beta);
                        f12_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_tmp(1.5_rt, eta, _beta);
                        f32_tmp.evaluate(1);

                        FermiIntegral<real_t> f12_pos_tmp(0.5_rt, eta_tilde_tmp, _beta);
                        f12_pos_tmp.evaluate(1);

                        FermiIntegral<real_t> f32_pos_tmp(1.5_rt, eta_tilde_tmp, _beta);
                        f32_pos_tmp.evaluate(1);

                        auto nd_tmp = get_n_derivs<real_t>(_beta, f12_tmp, f32_tmp, f12_pos_tmp, f32_pos_tmp);
                        return nd_tmp.dnp_deta;
                    }, beta, _h2);

                real_t err_beta = std::abs(nd.d2np_detadbeta - diff_beta) / std::abs(nd.d2np_detadbeta);

                util::threshold_println(err_eta,
                                        "eta = {:9.3f}, beta = {:8.3g},  ∂²n⁺/∂η∂β = {:12.5g},  D_η(∂n⁺/∂β) error = {:12.5g}  D_β(∂n⁺/∂η) error = {:12.5g}",
                                        eta, beta, nd.d2np_detadbeta, err_eta, err_beta);
            }
        }
    }

}


#include <cassert>
#include <cmath>
#include <print>

#include "real_type.H"
#include "fermi_integrals.H"
#include "difference_utils.H"
#include "electron_positron.H"
#include "util.H"

int main() {


    {

        const real_t h = 0.05_rt;

        util::green_println("testing dn⁻/dη");

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

                // check dn⁻/dη
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(_eta, beta);
                        return ne;
                    }, eta, _h);

                    real_t rel_err = std::abs(nd.dne_deta - diff) / std::abs(nd.dne_deta);

                    std::println("eta = {:9.3f}, beta = {:9.3g}, dn⁻/dη = {:15.8g}, error = {:15.8g}",
                                 eta, beta, nd.dne_deta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing dn⁺/dη");

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

                if (nd.dnp_deta == 0.0_rt) {
                    continue;
                }

                const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                // check dn⁻/dη
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto np = n_p_constraint<real_t>(_eta, beta);
                        return np;
                    }, eta, _h);

                    real_t rel_err = std::abs(nd.dnp_deta - diff) / std::abs(nd.dnp_deta);

                    std::println("eta = {:9.3f}, beta = {:9.3g}, dn⁺/dη = {:15.8g}, error = {:15.8g}",
                                 eta, beta, nd.dnp_deta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing dn⁻/dβ");

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

                // check dn⁻/dβ
                auto [diff, err] =
                    fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(eta, _beta);
                        return ne;
                    }, beta, _h);

                    real_t rel_err = std::abs(nd.dne_dbeta - diff) / std::abs(nd.dne_dbeta);

                    std::println("eta = {:9.3f}, beta = {:9.3g}, dn⁻/dβ = {:15.8g}, error = {:15.8g}",
                                 eta, beta, nd.dne_dbeta, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing d²n⁻/dη²");

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

                // check dn⁻/dη
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _eta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(_eta, beta);
                        return ne;
                    }, eta, _h);

                    real_t rel_err = std::abs(nd.d2ne_deta2 - diff) / std::abs(nd.d2ne_deta2);

                    std::println("eta = {:9.3f}, beta = {:9.3g}, d²n⁻/dη² = {:15.8g}, error = {:15.8g}",
                                 eta, beta, nd.d2ne_deta2, rel_err);
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println("");
        util::green_println("testing d²n⁻/dβ²");

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

                // check dn⁻/dβ
                auto [diff, err] =
                    fd::adaptive_diff2<real_t>([=] (real_t _beta) -> real_t
                    {
                        auto ne = n_e_constraint<real_t>(eta, _beta);
                        return ne;
                    }, beta, _h);

                    real_t rel_err = std::abs(nd.d2ne_dbeta2 - diff) / std::abs(nd.d2ne_dbeta2);

                    std::println("eta = {:9.3f}, beta = {:9.3g}, d²n⁻/dβ² = {:15.8g}, error = {:15.8g}",
                                 eta, beta, nd.d2ne_dbeta2, rel_err);
            }
        }
    }

}


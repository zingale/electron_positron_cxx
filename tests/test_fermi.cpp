#include <cassert>
#include <cmath>
#include <print>

#include "real_type.H"
#include "fermi_integrals.H"
#include "difference_utils.H"
#include "util.H"

int main() {

    // compare against the Gong et al. code

    {

        // k = 1/2

        const real_t k = 0.5_rt;

        {
            const real_t eta = -100_rt;
            const real_t beta = 0.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 3.2968314946796124e-044_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 3.2968314946796139e-044_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 1.2363118105048548e-044_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 3.2968314946796129e-044_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 1.2363118105048556e-044_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -7.7269488156553455e-045_rt) < 1.e-15);
        }

        {
            const real_t eta = -50_rt;
            const real_t beta = 1.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 2.2314386397062104e-022_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 2.2314386397062108e-022_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 4.4513616188460852e-023_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 2.2314386397062099e-022_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 4.4513616188460870e-023_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -1.0696566775877673e-023_rt) < 1.e-15);
        }

        {
            const real_t eta = 1.0_rt;
            const real_t beta = 10.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 4.3093121530612724_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 3.0909917562756424_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 0.19002958762134178_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 1.6796585258434078_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 0.13977966489586524_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -8.4816696742261236e-003_rt) < 1.e-15);
        }

        {
            const real_t eta = 500.0_rt;
            const real_t beta = 100.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 883930.45936891437_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 3535.6046159037460_rt) < 1.e-14);
            assert (util::rel_error(f.dF_dbeta, 4419.2988177917523_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 7.0710678132825819_rt) < 1.e-11);
            assert (util::rel_error(f.d2F_detadbeta, 17.677315986879442_rt) < 1.e-14);
            assert (util::rel_error(f.d2F_dbeta2, -22.094727366363866_rt) < 1.e-15);
        }

    }

    {

        // k = -1/2

        const real_t k = -0.5_rt;

        {
            const real_t eta = -100_rt;
            const real_t beta = 100.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 2.7816742731063666e-043_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 2.7816742731063678e-043_rt) <1.e-15);
            assert (util::rel_error(f.dF_dbeta, 1.2653970717385436e-045_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 2.7816742731063674e-043_rt) <1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 1.2653970717385444e-045_rt) <1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -5.9528644489672306e-048_rt) <1.e-15);
        }

        {
            const real_t eta = -50_rt;
            const real_t beta = 0.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 3.4186200954570750e-022_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 3.4186200954570741e-022_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 4.2732751193213414e-023_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 3.4186200954570736e-022_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 4.2732751193213450e-023_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -1.6024781697455028e-023_rt) < 1.e-15);
        }

        {
            const real_t eta = 100.0_rt;
            const real_t beta = 100.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 707.87776608227455_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 7.0717751162111329_rt) < 1.e-14);
            assert (util::rel_error(f.dF_dbeta, 3.5323860528718423_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, -7.0773544184610238e-006_rt) < 1.e-9);
            assert (util::rel_error(f.d2F_detadbeta, 3.5351802891431347e-002_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -1.7633986737239839e-002_rt) < 1.e-15);
        }

    }

    {
        // k = 3/2

        const real_t k = 1.5_rt;

        {
            const real_t eta = -75_rt;
            const real_t beta = 10.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 1.2553979904636453e-032_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 1.2553979904636447e-032_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 5.7230338002235925e-034_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 1.2553979904636453e-032_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 5.7230338002235968e-034_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -2.6212332799199379e-035_rt) < 1.e-15);
        }

        {
            const real_t eta = -20_rt;
            const real_t beta = 100.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 2.9294159670904935e-008_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 2.9294159663320433e-008_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 1.4502712363335996e-010_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 2.9294159648151430e-008_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 1.4502712359617485e-010_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -7.1804917485691793e-013_rt) < 1.e-15);
        }

        {
            const real_t eta = 40_rt;
            const real_t beta = 10000.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 1517805.2872690351_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 113369.99663886119_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 75.889697517711070_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 5656.8613205601887_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 5.6684715477425227_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -3.7944565338813361e-003_rt) < 1.e-15);
        }
    }

    {

        // k = 5/2

        const real_t k = 2.5_rt;

        {
            const real_t eta = -55_rt;
            const real_t beta = 15.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 2.1821367305243183e-023_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 2.1821367305243180e-023_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 6.9671437808138560e-025_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 2.1821367305243186e-023_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 6.9671437808138587e-025_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -2.2261771808312711e-026_rt) < 1.e-15);
        }

        {
            const real_t eta = 10.0_rt;
            const real_t beta = 1.e-4_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 1034.9073854351607_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 335.76592223107093_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 2231.0524130570820_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 80.071134602962601_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_detadbeta, 905.09777777902468_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -5201.1741381184684_rt) < 1.e-15);
        }

        {
            const real_t eta = 100.0_rt;
            const real_t beta = 1.0_rt;

            FermiIntegral<real_t> f(k, eta, beta);
            f.evaluate(2);

            assert (util::rel_error(f.F, 17946771.869476113_rt) < 1.e-15);
            assert (util::rel_error(f.dF_deta, 714843.05556091876_rt) < 1.e-15);
            assert (util::rel_error(f.dF_dbeta, 8740888.8928929009_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_deta2, 21361.250145668106_rt) < 1.e-14);
            assert (util::rel_error(f.d2F_detadbeta, 350417.80107778130_rt) < 1.e-15);
            assert (util::rel_error(f.d2F_dbeta2, -4257540.4021439729_rt) < 1.e-15);
        }
    }

    // test finite-difference approximations to derivatives

    {

        const real_t h = 0.05_rt;

        util::green_println("dF/dη");

        for (auto k : {-0.5_rt, 0.5_rt, 1.5_rt, 2.5_rt}) {
            for (auto eta : {-70.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
                for (auto beta : {1.e-7_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                    FermiIntegral<real_t> f(k, eta, beta);
                    f.evaluate(2);

                    const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                    // check dF/dη
                    auto [diff, err] =
                        fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                        {
                            FermiIntegral<real_t> _f(k, _eta, beta);
                            _f.evaluate(0);
                            return _f.F;
                        }, eta, _h);

                    real_t rel_err = std::abs(f.dF_deta - diff) / std::abs(f.dF_deta);

                    util::threshold_println(rel_err,
                                            "k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, dF/dη= {:15.8g}, error = {:15.8g}",
                                            k, eta, beta, f.dF_deta, rel_err);
                }
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println();
        util::green_println("dF/dβ");

        for (auto k : {-0.5_rt, 0.5_rt, 1.5_rt, 2.5_rt}) {
            for (auto eta : {-70.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
                for (auto beta : {1.e-7_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                    FermiIntegral<real_t> f(k, eta, beta);
                    f.evaluate(2);

                    const real_t _h = (beta == 0) ? 1.e-4_rt : h * std::abs(beta);

                    // check dF/dβ
                    auto [diff, err] =
                        fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                        {
                            FermiIntegral<real_t> _f(k, eta, _beta);
                            _f.evaluate(0);
                            return _f.F;
                        }, beta, _h);

                    real_t rel_err = std::abs(f.dF_dbeta - diff) / std::abs(f.dF_dbeta);

                    util::threshold_println(rel_err,
                                            "k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, dF/dβ = {:15.8g}, error = {:15.8g}",
                                            k, eta, beta, f.dF_dbeta, rel_err);
                }
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println();
        util::green_println("d²F/dη²");

        for (auto k : {-0.5_rt, 0.5_rt, 1.5_rt, 2.5_rt}) {
            for (auto eta : {-70.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
                for (auto beta : {1.e-7_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                    FermiIntegral<real_t> f(k, eta, beta);
                    f.evaluate(2);

                    const real_t _h = (eta == 0) ? h : h * std::abs(eta);

                    // check d²F/dη²

                    // first do a second deriv difference on F
                    auto [diff, err] =
                        fd::adaptive_diff2<real_t>([=] (real_t _eta) -> real_t
                        {
                            FermiIntegral<real_t> _f(k, _eta, beta);
                            _f.evaluate(0);
                            return _f.F;
                        }, eta, _h);

                    real_t rel_err = std::abs(f.d2F_deta2 - diff) / std::abs(f.d2F_deta2);

                    // next do a first order deriv on dF/dη
                    auto [diff2, err2] =
                        fd::adaptive_diff<real_t>([=] (real_t _eta) -> real_t
                        {
                            FermiIntegral<real_t> _f(k, _eta, beta);
                            _f.evaluate(1);
                            return _f.dF_deta;
                        }, eta, _h);

                    real_t rel_err2 = std::abs(f.d2F_deta2 - diff2) / std::abs(f.d2F_deta2);


                    util::threshold_println(rel_err,
                                            "k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, d²F/dη² = {:15.8g}, error (D2F) = {:15.8g}, error (DF') = {:15.8g}",
                                            k, eta, beta, f.dF_dbeta, rel_err, rel_err2);
                }
            }
        }
    }

    {

        const real_t h = 0.05_rt;

        std::println();
        util::green_println("d²F/dβ²");

        for (auto k : {-0.5_rt, 0.5_rt, 1.5_rt, 2.5_rt}) {
            for (auto eta : {-70.0_rt, 0.0_rt, 50.0_rt, 500.0_rt, 10000.0_rt}) {
                for (auto beta : {1.e-7_rt, 1.e-3_rt, 30.0_rt, 100.0_rt}) {

                    FermiIntegral<real_t> f(k, eta, beta);
                    f.evaluate(2);

                    const real_t _h = (beta == 0) ? h : h * std::abs(beta);

                    // check d²F/dβ²

                    // first do a second deriv difference on F
                    auto [diff, err] =
                        fd::adaptive_diff2<real_t>([=] (real_t _beta) -> real_t
                        {
                            FermiIntegral<real_t> _f(k, eta, _beta);
                            _f.evaluate(0);
                            return _f.F;
                        }, beta, _h);

                    real_t rel_err = std::abs(f.d2F_dbeta2 - diff) / std::abs(f.d2F_dbeta2);

                    // next do a second order deriv on dF/dβ
                    auto [diff2, err2] =
                        fd::adaptive_diff<real_t>([=] (real_t _beta) -> real_t
                        {
                            FermiIntegral<real_t> _f(k, eta, _beta);
                            _f.evaluate(1);
                            return _f.dF_dbeta;
                        }, beta, _h);

                    real_t rel_err2 = std::abs(f.d2F_dbeta2 - diff2) / std::abs(f.d2F_dbeta2);

                    util::threshold_println(rel_err,
                                            "k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, d²F/dβ² = {:15.8g}, error (D2F) = {:15.8g}, error (DF') = {:15.8g}",
                                            k, eta, beta, f.dF_dbeta, rel_err, rel_err2);
                }
            }
        }
    }
}

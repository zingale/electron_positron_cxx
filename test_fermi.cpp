#include <cassert>
#include <print>

#include "fermi_integrals.H"
#include "difference_utils.H"


template <typename T>
T rel_error(T a, T b) {
    T err = std::abs(a - b) / std::abs(b);
    return err;
}

int main() {

    using RealT = long double;


    // compare against the Gong et al. code

    {

        // k = 1/2

        RealT k = 0.5L;

        {
            RealT eta = -100L;
            RealT beta = 0.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 3.2968314946796124e-044L) < 1.e-15);
            assert (rel_error(f.dF_deta, 3.2968314946796139e-044L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 1.2363118105048548e-044L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 3.2968314946796129e-044L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 1.2363118105048556e-044L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -7.7269488156553455e-045L) < 1.e-15);
        }

        {
            RealT eta = -50L;
            RealT beta = 1.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 2.2314386397062104e-022L) < 1.e-15);
            assert (rel_error(f.dF_deta, 2.2314386397062108e-022L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 4.4513616188460852e-023L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 2.2314386397062099e-022L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 4.4513616188460870e-023L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -1.0696566775877673e-023L) < 1.e-15);
        }

        {
            RealT eta = 1.0L;
            RealT beta = 10.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 4.3093121530612724L) < 1.e-15);
            assert (rel_error(f.dF_deta, 3.0909917562756424L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 0.19002958762134178L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 1.6796585258434078L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 0.13977966489586524L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -8.4816696742261236e-003L) < 1.e-15);
        }

        {
            RealT eta = 500.0L;
            RealT beta = 100.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 883930.45936891437L) < 1.e-15);
            assert (rel_error(f.dF_deta, 3535.6046159037460L) < 1.e-14);
            assert (rel_error(f.dF_dbeta, 4419.2988177917523L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 7.0710678132825819L) < 1.e-11);
            assert (rel_error(f.d2F_detadbeta, 17.677315986879442L) < 1.e-14);
            assert (rel_error(f.d2F_dbeta2, -22.094727366363866L) < 1.e-15);
        }

    }

    {

        // k = -1/2

        RealT k = -0.5L;

        {
            RealT eta = -100L;
            RealT beta = 100.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 2.7816742731063666e-043L) < 1.e-15);
            assert (rel_error(f.dF_deta, 2.7816742731063678e-043L) <1.e-15);
            assert (rel_error(f.dF_dbeta, 1.2653970717385436e-045L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 2.7816742731063674e-043L) <1.e-15);
            assert (rel_error(f.d2F_detadbeta, 1.2653970717385444e-045L) <1.e-15);
            assert (rel_error(f.d2F_dbeta2, -5.9528644489672306e-048L) <1.e-15);
        }

        {
            RealT eta = -50L;
            RealT beta = 0.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 3.4186200954570750e-022L) < 1.e-15);
            assert (rel_error(f.dF_deta, 3.4186200954570741e-022L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 4.2732751193213414e-023L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 3.4186200954570736e-022L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 4.2732751193213450e-023L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -1.6024781697455028e-023L) < 1.e-15);
        }

        {
            RealT eta = 100.0L;
            RealT beta = 100.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 707.87776608227455L) < 1.e-15);
            assert (rel_error(f.dF_deta, 7.0717751162111329L) < 1.e-14);
            assert (rel_error(f.dF_dbeta, 3.5323860528718423L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, -7.0773544184610238e-006L) < 1.e-9);
            assert (rel_error(f.d2F_detadbeta, 3.5351802891431347e-002L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -1.7633986737239839e-002L) < 1.e-15);
        }

    }

    {
        // k = 3/2

        RealT k = 1.5L;

        {
            RealT eta = -75L;
            RealT beta = 10.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 1.2553979904636453e-032L) < 1.e-15);
            assert (rel_error(f.dF_deta, 1.2553979904636447e-032L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 5.7230338002235925e-034L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 1.2553979904636453e-032L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 5.7230338002235968e-034L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -2.6212332799199379e-035L) < 1.e-15);
        }

        {
            RealT eta = -20L;
            RealT beta = 100.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 2.9294159670904935e-008L) < 1.e-15);
            assert (rel_error(f.dF_deta, 2.9294159663320433e-008L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 1.4502712363335996e-010L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 2.9294159648151430e-008L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 1.4502712359617485e-010L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -7.1804917485691793e-013L) < 1.e-15);
        }

        {
            RealT eta = 40L;
            RealT beta = 10000.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 1517805.2872690351L) < 1.e-15);
            assert (rel_error(f.dF_deta, 113369.99663886119L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 75.889697517711070L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 5656.8613205601887L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 5.6684715477425227L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -3.7944565338813361e-003L) < 1.e-15);
        }
    }

    {

        // k = 5/2

        RealT k = 2.5L;

        {
            RealT eta = -55L;
            RealT beta = 15.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 2.1821367305243183e-023L) < 1.e-15);
            assert (rel_error(f.dF_deta, 2.1821367305243180e-023L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 6.9671437808138560e-025L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 2.1821367305243186e-023L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 6.9671437808138587e-025L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -2.2261771808312711e-026L) < 1.e-15);
        }

        {
            RealT eta = 10.0L;
            RealT beta = 1.e-4L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 1034.9073854351607L) < 1.e-15);
            assert (rel_error(f.dF_deta, 335.76592223107093L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 2231.0524130570820L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 80.071134602962601L) < 1.e-15);
            assert (rel_error(f.d2F_detadbeta, 905.09777777902468L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -5201.1741381184684L) < 1.e-15);
        }

        {
            RealT eta = 100.0L;
            RealT beta = 1.0L;

            FermiIntegral<RealT> f(k, eta, beta);
            f.evaluate(2);

            assert (rel_error(f.F, 17946771.869476113L) < 1.e-15);
            assert (rel_error(f.dF_deta, 714843.05556091876L) < 1.e-15);
            assert (rel_error(f.dF_dbeta, 8740888.8928929009L) < 1.e-15);
            assert (rel_error(f.d2F_deta2, 21361.250145668106L) < 1.e-14);
            assert (rel_error(f.d2F_detadbeta, 350417.80107778130L) < 1.e-15);
            assert (rel_error(f.d2F_dbeta2, -4257540.4021439729L) < 1.e-15);
        }
    }

    // test finite-difference approximations to derivatives

    RealT h = 0.05L;

    std::println("dF/dη");

    for (auto k : {-0.5L, 0.5L, 1.5L, 2.5L}) {
        for (auto eta : {-70.0L, 0.0L, 50.0L, 500.0L, 10000.0L}) {
            for (auto beta : {1.e-7L, 1.e-3L, 30.0L, 100.0L}) {

                FermiIntegral<RealT> f(k, eta, beta);
                f.evaluate(2);

                RealT _h = (eta == 0) ? h : h * std::abs(eta);

                // check dF/dη
                auto [diff, err] =
                    fd::adaptive_diff<RealT>([=] (RealT _eta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, _eta, beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, eta, _h);

                RealT rel_err = std::abs(f.dF_deta - diff) / std::abs(f.dF_deta);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, dF/dη= {:15.8g}, error = {:15.8g}",
                             k, eta, beta, f.dF_deta, rel_err);
            }
        }
    }

    std::println();
    std::println("dF/dβ");

    for (auto k : {-0.5L, 0.5L, 1.5L, 2.5L}) {
        for (auto eta : {-70.0L, 0.0L, 50.0L, 500.0L, 10000.0L}) {
            for (auto beta : {1.e-7L, 1.e-3L, 30.0L, 100.0L}) {

                FermiIntegral<RealT> f(k, eta, beta);
                f.evaluate(2);

                RealT _h = (beta == 0) ? 1.e-4L : h * std::abs(beta);

                // check dF/dβ
                auto [diff, err] =
                    fd::adaptive_diff<RealT>([=] (RealT _beta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, eta, _beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, beta, _h);

                RealT rel_err = std::abs(f.dF_dbeta - diff) / std::abs(f.dF_dbeta);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, dF/dβ = {:15.8g}, error = {:15.8g}",
                             k, eta, beta, f.dF_dbeta, rel_err);
            }
        }
    }

    std::println();
    std::println("d²F/dη²");

    for (auto k : {-0.5L, 0.5L, 1.5L, 2.5L}) {
        for (auto eta : {-70.0L, 0.0L, 50.0L, 500.0L, 10000.0L}) {
            for (auto beta : {1.e-7L, 1.e-3L, 30.0L, 100.0L}) {

                FermiIntegral<RealT> f(k, eta, beta);
                f.evaluate(2);

                RealT _h = (eta == 0) ? h : h * std::abs(eta);

                // check dF/dβ
                auto [diff, err] =
                    fd::adaptive_diff2<RealT>([=] (RealT _eta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, _eta, beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, eta, _h);

                RealT rel_err = std::abs(f.d2F_deta2 - diff) / std::abs(f.d2F_deta2);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, d²F/dη² = {:15.8g}, error = {:15.8g}",
                             k, eta, beta, f.dF_dbeta, rel_err);
            }
        }
    }

    std::println();
    std::println("d²F/dβ²");

    for (auto k : {-0.5L, 0.5L, 1.5L, 2.5L}) {
        for (auto eta : {-70.0L, 0.0L, 50.0L, 500.0L, 10000.0L}) {
            for (auto beta : {1.e-7L, 1.e-3L, 30.0L, 100.0L}) {

                FermiIntegral<RealT> f(k, eta, beta);
                f.evaluate(2);

                RealT _h = (beta == 0) ? h : h * std::abs(beta);

                // check dF/dβ
                auto [diff, err] =
                    fd::adaptive_diff2<RealT>([=] (RealT _beta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, eta, _beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, beta, _h);

                RealT rel_err = std::abs(f.d2F_dbeta2 - diff) / std::abs(f.d2F_dbeta2);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:9.3g}, d²F/dβ² = {:15.8g}, error = {:15.8g}",
                             k, eta, beta, f.dF_dbeta, rel_err);
            }
        }
    }

}

#include <iostream>
#include <print>

#include "fermi_integrals.H"
#include "difference_utils.H"

int main() {

    using RealT = long double;

    RealT h = 1.e-2L;

    std::println("dF/dη");

    for (auto k : {-0.5L, 0.5L, 1.5L, 2.5L}) {
        for (auto eta : {-70.0L, 0.0L, 50.0L, 500.0L, 10000.0L}) {
            for (auto beta : {1.e-7L, 1.e-3L, 30.0L, 100.0L}) {

                FermiIntegral<RealT> f(k, eta, beta);
                f.evaluate(2);

                RealT _h = (eta == 0) ? h : h * std::abs(eta);

                // check dF/dη
                auto [diff, err] =
                    adaptive_diff<RealT>([=] (RealT _eta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, _eta, beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, eta, _h);

                RealT rel_err = std::abs(f.dF_deta - diff) / std::abs(f.dF_deta);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:7.3f}, dF/dη= {:15.8g}, error = {:15.8g}",
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
                    adaptive_diff<RealT>([=] (RealT _beta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, eta, _beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, beta, _h);

                RealT rel_err = std::abs(f.dF_dbeta - diff) / std::abs(f.dF_dbeta);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:7.3f}, dF/dβ = {:15.8g}, error = {:15.8g}",
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
                    adaptive_diff2<RealT>([=] (RealT _eta) -> RealT
                    {
                        FermiIntegral<RealT> _f(k, _eta, beta);
                        _f.evaluate(0);
                        return _f.F;
                    }, eta, _h);

                RealT rel_err = std::abs(f.d2F_deta2 - diff) / std::abs(f.d2F_deta2);

                std::println("k = {:5.2f}, eta = {:9.3f}, beta = {:7.3f}, d²F/dη² = {:15.8g}, error = {:15.8g}",
                             k, eta, beta, f.dF_dbeta, rel_err);
            }
        }
    }

}

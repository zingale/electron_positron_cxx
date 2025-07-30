#include <iostream>
#include <print>

#include "electron_positron.H"
#include "difference_utils.H"

using RealT = long double;


void
test_ne_rho_derivs() {

    ElectronPositronEOS<RealT> eos;
    constexpr RealT Ye{0.5L};

    constexpr RealT eps{0.01};

    for (auto T : {1.e4L, 1.e6L, 1.e9L}) {
        for (auto rho : {1.e-2L, 1.e2L, 1.e5L, 1.e9L}) {
            auto es = eos.pe_state(rho, T, Ye);
            auto drho{eps * rho};
            auto [deriv, err] = fd::adaptive_diff<RealT>([&] (RealT _rho) -> RealT
                {
                    auto es_eps = eos.pe_state(_rho, T, Ye);
                    return es_eps.n_e;
                }, rho, drho);

            RealT rel_err = std::abs(es.dne_drho - deriv) / std::abs(es.dne_drho);
            std::println("ρ = {:10} T = {:10}, ∂n⁻/∂ρ = {:15.8g}, error = {:15.5g}",
                         rho, T, es.dne_drho, rel_err);
        }
    }

}


int main() {

    test_ne_rho_derivs();

}

#include <cassert>
#include <iostream>
#include <print>

#include "real_type.H"
#include "helmholtz.H"


auto main() -> int
{

    std::println("Enter rho, T, Ye (space-separated): ");
    double rho_{};
    double T_{};
    double Ye_{};

    std::cin >> rho_ >> T_ >> Ye_;

    assert (rho_ != 0.0 && T_ != 0.0 && Ye_ != 0.0);

    auto rho = static_cast<real_t>(rho_);
    auto T = static_cast<real_t>(T_);
    auto Ye = static_cast<real_t>(Ye_);

    auto helm = get_helmholtz_terms(rho, T, Ye);

    std::cout << helm << std::endl;

}

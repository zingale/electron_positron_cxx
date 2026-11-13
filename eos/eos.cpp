#include <cassert>
#include <iostream>
#include <print>
#include <string>

#include "real_type.H"
#include "electron_positron.H"


auto main() -> int
{

    std::println("Enter rho, T, Ye (space-separated): ");
    std::string rho_;
    std::string T_;
    std::string Ye_;

    std::cin >> rho_ >> T_ >> Ye_;

    const real_t rho = str_to_real_t(rho_.c_str());
    const real_t T = str_to_real_t(T_.c_str());
    const real_t Ye = str_to_real_t(Ye_.c_str());

    assert (rho != 0.0_rt && T != 0.0_rt && Ye != 0.0_rt);

    ElectronPositronEOS<real_t> eos;
    auto state = eos.pe_state(rho, T, Ye);

    std::cout << state << std::endl;

}

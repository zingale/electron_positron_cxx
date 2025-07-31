#include <cassert>
#include <iostream>
#include <print>

#include "real_type.H"
#include "electron_positron.H"


int main() {

    std::println("Enter rho, T, Ye (space-separated): ");
    double _rho{};
    double _T{};
    double _Ye{};

    std::cin >> _rho >> _T >> _Ye;

    assert (_rho != 0.0 && _T != 0.0 && _Ye != 0.0);

    real_t rho = static_cast<real_t>(_rho);
    real_t T = static_cast<real_t>(_T);
    real_t Ye = static_cast<real_t>(_Ye);

    std::println("{} {} {}", rho, T, Ye);

    ElectronPositronEOS<real_t> eos;
    auto state = eos.pe_state(rho, T, Ye);

    std::cout << state << std::endl;

}

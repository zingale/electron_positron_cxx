#include <iostream>

#include "electron_positron.H"

int main() {

    using RealT = long double;

    ElectronPositronEOS<RealT> eos;

    auto es = eos.pe_state(1.e4L, 1.e7L, 0.5L);

    std::cout << es << std::endl;

}

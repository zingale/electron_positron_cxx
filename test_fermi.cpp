#include "fermi_integrals.H"

#include <iostream>

int main() {

    FermiIntegral<long double> f(0.5, 500, 100);
    f.evaluate(2);

    std::cout << f << std::endl;

}

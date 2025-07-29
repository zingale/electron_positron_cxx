#include <print>
#include "brent.H"


int main() {

    using RealT = long double;

    RealT a = -4.0L;
    RealT b = 4.0L / 3.0L;

    auto r = brent<RealT>([=] (RealT x) -> RealT
                          {return (x + 3.0L) * std::pow(x - 1.0L, 2);},
        a, b);

    std::println("root = {}", r);
}

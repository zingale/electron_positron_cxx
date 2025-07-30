#include <print>

#include "brent.H"
#include "real_type.H"

int main() {

    real_t a = -4.0_rt;
    real_t b = 4.0_rt / 3.0_rt;

    auto r = brent<real_t>([=] (real_t x) -> real_t
                          {return (x + 3.0_rt) * std::pow(x - 1.0_rt, 2);},
        a, b);

    std::println("root = {}", r);
}


#include "brent.H"
#include "real_type.H"
#include "util.H"
#include "mp_math.H"

int main() {

    const real_t a = -4.0_rt;
    const real_t b = 4.0_rt / 3.0_rt;

    auto r = brent<real_t>([=] (real_t x) -> real_t
                          {return (x + 3.0_rt) * mp::pow(x - 1.0_rt, 2);},
        a, b);

    util::println("root = {}", r);
}

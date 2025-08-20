#include <limits>

#include "real_type.H"
#include "util.H"

int main() {

    // this will be interpreted as a long double regardless of real_t,
    // so we might get error in the later bits when building with 128-
    // or 256-bit precision.

    util::println("size of real_t is {} bytes", sizeof(real_t));
    util::println("machine epsilon is {}", std::numeric_limits<real_t>::epsilon());
    util::println("minimum exponent is 10**{}", std::numeric_limits<real_t>::min_exponent10);
    util::println("maximum exponent is 10**{}", std::numeric_limits<real_t>::max_exponent10);
    util::println("");

    real_t x{1.2345_rt};
    util::println("{}", x);

    // this should work well regardless of precision

    real_t y = str_to_real_t("1.2345");
    util::println("{}", y);

    util::green_println("this is a test {:20.10g}", y);

    real_t z = 1.e-16_rt;
    util::threshold_println(y, "this is high {}", y);
    util::threshold_println(z, "this is small {}", z);

}


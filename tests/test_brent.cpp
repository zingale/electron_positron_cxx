#include <limits>
#include <stdexcept>

#include "brent.H"
#include "real_type.H"
#include "util.H"
#include "mp_math.H"

auto main() -> int
{

    const real_t a = -4.0_rt;
    const real_t b = 4.0_rt / 3.0_rt;

    auto r = brent<real_t>([=] (real_t x) -> real_t
                          {return (x + 3.0_rt) * mp::pow(x - 1.0_rt, 2);},
        a, b);

    util::println("root = {}", r);

    // Scaling the residual must not change the root or bracket validation.
    const real_t small_scale = mp::sqrt(std::numeric_limits<real_t>::min()) *
        std::numeric_limits<real_t>::epsilon();
    for (auto scale : {small_scale, -small_scale}) {
        const auto root = brent<real_t>([=] (real_t x) -> real_t
                                       {return scale * (x * x - 0.25_rt);},
                                       0.0_rt, 1.0_rt);
        if (mp::abs(root - 0.5_rt) > 16.0_rt * std::numeric_limits<real_t>::epsilon()) {
            throw std::runtime_error("Brent solve failed for tiny residuals");
        }

        bool rejected = false;
        try {
            brent<real_t>([=] (real_t x) -> real_t
                          {return scale * (x * x + 1.0_rt);},
                          0.0_rt, 1.0_rt);
        } catch (const std::out_of_range&) {
            rejected = true;
        }
        if (!rejected) {
            throw std::runtime_error("Brent accepted a bracket without a root");
        }
    }
}

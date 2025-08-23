#include "real_type.H"
#include "mp_math.H"
#include "difference_utils.H"
#include "util.H"

template <typename T>
auto f1(T x) -> T
{
    return mp::exp(x) / (mp::sin(x) - x * x);
}

template <typename T>
auto df1dx(T x) -> T
{
    return mp::exp(x) * (-x * x + 2.0_rt * x - constants::sqrt2 * mp::cos(x + 0.25_rt * constants::pi)) /
        mp::pow(x * x - mp::sin(x), 2);
}

template <typename T>
auto d2f1dx2(T x) -> T
{
    return (-2.0_rt * mp::pow(2.0_rt * x - mp::cos(x), 2) +
           2.0_rt * (2.0_rt * x - mp::cos(x)) * (x * x - mp::sin(x)) -
           mp::pow(x * x - mp::sin(x), 2) + (x * x - mp::sin(x)) *
            (mp::sin(x) + 2.0_rt)) *
        mp::exp(x) / mp::pow(x * x - mp::sin(x), 3);
}

auto main() -> int 
{

    const real_t x0 = 1.0_rt;
    const real_t h = 1.e-3_rt;

    auto actual_deriv = df1dx<real_t>(x0);

    // test the 4th order difference
    auto deriv4 = fd::fourth_order_diff<real_t>(f1<real_t>, x0, h);

    util::println("4th order diff, error = {}", mp::abs(deriv4 - actual_deriv) / mp::abs(actual_deriv));

    // test the 6th order difference
    auto deriv6 = fd::sixth_order_diff<real_t>(f1<real_t>, x0, h);

    util::println("6th order diff, error = {}", mp::abs(deriv6 - actual_deriv) / mp::abs(actual_deriv));

    // test the 8th order difference
    auto deriv8 = fd::eighth_order_diff<real_t>(f1<real_t>, x0, h);

    util::println("8th order diff, error = {}", mp::abs(deriv8 - actual_deriv) / mp::abs(actual_deriv));

    // test the adaptive method
    auto [deriv_adaptive, err] = fd::adaptive_diff<real_t>(f1<real_t>, x0, h);

    util::println("adaptive diff, error = {}", mp::abs(deriv_adaptive - actual_deriv) / mp::abs(actual_deriv));

    // test the adaptive second-deriv method
    auto [deriv_adaptive2, err2] = fd::adaptive_diff2<real_t>(f1<real_t>, x0, h);

    auto actual_deriv2 = d2f1dx2(x0);

    util::println("adaptive second diff, error = {}", mp::abs(deriv_adaptive2 - actual_deriv2) / mp::abs(actual_deriv2));

}

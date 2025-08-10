#include <cmath>
#include <print>

#include "real_type.H"
#include "difference_utils.H"


template <typename T>
T f1(T x) {
    return std::exp(x) / (std::sin(x) - x * x);
}

template <typename T>
T df1dx(T x) {
    return std::exp(x) * (-x * x + 2.0_rt * x - sqrt2 * std::cos(x + 0.25_rt * pi)) /
        std::pow(x * x - std::sin(x), 2);
}

template <typename T>
T d2f1dx2(T x) {
    return (-2.0_rt * std::pow(2.0_rt * x - std::cos(x), 2) +
           2.0_rt * (2.0_rt * x - std::cos(x)) * (x * x - std::sin(x)) -
           std::pow(x * x - std::sin(x), 2) + (x * x - std::sin(x)) *
            (std::sin(x) + 2.0_rt)) *
        std::exp(x) / std::pow(x * x - std::sin(x), 3);
}

int main() {

    const real_t x0 = 1.0_rt;
    const real_t h = 1.e-3_rt;

    auto actual_deriv = df1dx<real_t>(x0);

    // test the 4th order difference
    auto deriv4 = fd::fourth_order_diff<real_t>(f1<real_t>, x0, h);

    std::println("4th order diff, error = {}", std::abs(deriv4 - actual_deriv) / std::abs(actual_deriv));

    // test the 6th order difference
    auto deriv6 = fd::sixth_order_diff<real_t>(f1<real_t>, x0, h);

    std::println("6th order diff, error = {}", std::abs(deriv6 - actual_deriv) / std::abs(actual_deriv));

    // test the 8th order difference
    auto deriv8 = fd::eighth_order_diff<real_t>(f1<real_t>, x0, h);

    std::println("8th order diff, error = {}", std::abs(deriv8 - actual_deriv) / std::abs(actual_deriv));

    // test the adaptive method
    auto [deriv_adaptive, err] = fd::adaptive_diff<real_t>(f1<real_t>, x0, h);

    std::println("adaptive diff, error = {}", std::abs(deriv_adaptive - actual_deriv) / std::abs(actual_deriv));

    // test the adaptive second-deriv method
    auto [deriv_adaptive2, err2] = fd::adaptive_diff2<real_t>(f1<real_t>, x0, h);

    auto actual_deriv2 = d2f1dx2(x0);

    std::println("adaptive second diff, error = {}", std::abs(deriv_adaptive2 - actual_deriv2) / std::abs(actual_deriv2));

}

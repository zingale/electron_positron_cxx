#include <cmath>
#include <numbers>
#include <print>

#include "difference_utils.H"


template <typename T>
T f1(T x) {
    return std::exp(x) / (std::sin(x) - x * x);
}

template <typename T>
T f1prime(T x) {
    auto pi = std::numbers::pi_v<T>;
    return std::exp(x) * (-x * x + 2.0L * x - std::sqrt(2.0L) * std::cos(x + 0.25L * pi)) /
        std::pow(x * x - std::sin(x), 2);
}

int main() {

    using RealT = long double;

    RealT x0 = 1.0L;
    RealT h = 1.e-3L;

    auto actual_deriv = f1prime<RealT>(x0);

    // test the 4th order difference
    auto deriv4 = fourth_order_diff<RealT>(f1<RealT>, x0, h);

    std::println("4th order diff, error = {}", std::abs(deriv4 - actual_deriv) / std::abs(actual_deriv));

    // test the 6th order difference
    auto deriv6 = sixth_order_diff<RealT>(f1<RealT>, x0, h);

    std::println("6th order diff, error = {}", std::abs(deriv6 - actual_deriv) / std::abs(actual_deriv));

    // test the 8th order difference
    auto deriv8 = eighth_order_diff<RealT>(f1<RealT>, x0, h);

    std::println("8th order diff, error = {}", std::abs(deriv8 - actual_deriv) / std::abs(actual_deriv));

    // test the adaptive method
    auto [deriv_adaptive, err] = adaptive_diff<RealT>(f1<RealT>, x0, h);

    std::println("adaptive diff, error = {}", std::abs(deriv_adaptive - actual_deriv) / std::abs(actual_deriv));
}

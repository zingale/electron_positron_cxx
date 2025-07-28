#include <cmath>
#include <numbers>
#include <print>

#include "difference_utils.H"


template <typename T>
T f1(T x) {
    return std::exp(x) / (std::sin(x) - x * x);
}

template <typename T>
T df1dx(T x) {
    auto pi = std::numbers::pi_v<T>;
    return std::exp(x) * (-x * x + 2.0L * x - std::sqrt(2.0L) * std::cos(x + 0.25L * pi)) /
        std::pow(x * x - std::sin(x), 2);
}

template <typename T>
T d2f1dx2(T x) {
    return (-2.0L * std::pow(2.0L * x - std::cos(x), 2) +
           2.0L * (2.0L * x - std::cos(x)) * (x * x - std::sin(x)) -
           std::pow(x * x - std::sin(x), 2) + (x * x - std::sin(x)) *
            (std::sin(x) + 2.0L)) *
        std::exp(x) / std::pow(x * x - std::sin(x), 3);
}

int main() {

    using RealT = long double;

    RealT x0 = 1.0L;
    RealT h = 1.e-3L;

    auto actual_deriv = df1dx<RealT>(x0);

    // test the 4th order difference
    auto deriv4 = fd::fourth_order_diff<RealT>(f1<RealT>, x0, h);

    std::println("4th order diff, error = {}", std::abs(deriv4 - actual_deriv) / std::abs(actual_deriv));

    // test the 6th order difference
    auto deriv6 = fd::sixth_order_diff<RealT>(f1<RealT>, x0, h);

    std::println("6th order diff, error = {}", std::abs(deriv6 - actual_deriv) / std::abs(actual_deriv));

    // test the 8th order difference
    auto deriv8 = fd::eighth_order_diff<RealT>(f1<RealT>, x0, h);

    std::println("8th order diff, error = {}", std::abs(deriv8 - actual_deriv) / std::abs(actual_deriv));

    // test the adaptive method
    auto [deriv_adaptive, err] = fd::adaptive_diff<RealT>(f1<RealT>, x0, h);

    std::println("adaptive diff, error = {}", std::abs(deriv_adaptive - actual_deriv) / std::abs(actual_deriv));

    // test the adaptive second-deriv method
    auto [deriv_adaptive2, err2] = fd::adaptive_diff2<RealT>(f1<RealT>, x0, h);

    auto actual_deriv2 = d2f1dx2(x0);

    std::println("adaptive second diff, error = {}", std::abs(deriv_adaptive2 - actual_deriv2) / std::abs(actual_deriv2));

}

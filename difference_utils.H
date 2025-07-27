#ifndef DIFFERENCE_UTILS_H
#define DIFFERENCE_UTILS_H

#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>

// Some high-order finite-difference approximations for the EOS.

// the coefficients for these can be found in many sources, including:
// https://en.wikipedia.org/wiki/Finite_difference_coefficient

template <typename T>
T fourth_order_diff(std::function<T(T)> func, T x0, T delta) {
    // Compute a 4th order accurate centered difference approximation
    // of a function, and allow us to specify the component of the
    // object that is returned (if applicable)

    std::vector<T> fvals;

    for (auto i : {-2, -1, 0, 1, 2}) {
        if (i == 0) {
            fvals.push_back(0.0L);
            continue;
        }

        auto _f = func(x0 + static_cast<T>(i) * delta);
        fvals.push_back(_f);
    }

    T deriv = ((fvals[0] - fvals[4]) +
               8.0L * (-fvals[1] + fvals[3]))
        / (12L * delta);

    return deriv;
}

template <typename T>
T sixth_order_diff(std::function<T(T)> func, T x0, T delta) {
    // Compute a 6th order accurate centered difference approximation
    // of a function, and allow us to specify the component of the
    // object that is returned (if applicable)

    std::vector<T> fvals;

    for (auto i : {-3, -2, -1, 0, 1, 2, 3}) {
        if (i == 0) {
            fvals.push_back(0.0L);
            continue;
        }

        auto _f = func(x0 + static_cast<T>(i) * delta);
        fvals.push_back(_f);
    }

    T deriv = ((-fvals[0] + fvals[6]) +
               9.0 * (fvals[1] - fvals[5]) +
               45.0 * (-fvals[2] + fvals[4])) /
        (60 * delta);

    return deriv;
}

template <typename T>
T eighth_order_diff(std::function<T(T)> func, T x0, T delta) {
    // Compute an 8th order accurate centered difference approximation
    // of a function, and allow us to specify the component of the
    // object that is returned (if applicable)

    std::vector<T> fvals;

    for (auto i : {-4, -3, -2, -1, 0, 1, 2, 3, 4}) {
        if (i == 0) {
            fvals.push_back(0.0L);
            continue;
        }

        auto _f = func(x0 + static_cast<T>(i) * delta);
        fvals.push_back(_f);
    }

    T deriv = ((1.0 / 280.0) * (fvals[0] - fvals[8]) +
               (4.0 / 105.0) * (-fvals[1] + fvals[7]) +
               0.2 * (fvals[2] - fvals[6]) +
               0.8 * (-fvals[3] + fvals[5])) / delta;

    return deriv;
}

template <typename T>
std::pair<T, T>
adaptive_diff(std::function<T(T)> func, T x0, T h, int nlevels=10) {
    // Perform an adaptive centered-difference estimate to the first
    // derivative using Richardson-extrapolation / Ridders' method.

    // see C.J.F. Ridders, "Accurate computation of F'x(x) and F''(x)
    // Advances in Engineering Software, 4, 2, 1978

    // we'll use a second-order centered difference, starting with the
    // initial h.  this has the form:
    //
    //  D_h(f) = 0.5 * (f(x0 + h) - f(x0 - h)) / h + c (h**2) + O(h**4)
    //
    // the key point is that C depends only on the function f (or its
    // derivatives) at point x0.  That means we can compute D_{h/2}(f)
    // and then combine the 2 for a better estimate of the derivative:
    //
    // D' = (4 D_{h/2}(f) - D_h(f)) / 3 + O(h**4)

    // In Ridders' method, we build a tableau of the form:
    //
    //
    //   A(1, 1)  A(1, 2)  A(1, 3)  A(1, 4) ...
    //            A(2, 2)  A(2, 3)  A(2, 4) ...
    //                     A(3, 3)  A(3, 3) ...
    //                              A(4, 4) ...
    //
    // where the entries in each row come via extrapolation of the data
    // in the previous row.
    //
    // In this notation, A(1, 1) = D_h(f),
    //                   A(1, 2) = D_{h/2}(f)
    //                   A(2, 2) = D'
    //
    // we then have:
    //
    // A(1, n) = 0.5 * (f(x + h/2**(n-1)) - f(x - h/2**(n-1))) / (h/2**(n-1))
    //
    // A(m, n) = (4**(m-1) A(m-1, n) - A(m-1, n-1)) / (4**(m-1) - 1)

    // Note: we use 0-based indexing below

    constexpr int max_levels = 50;

    assert (nlevels < max_levels);

    std::array<std::array<T, max_levels>, max_levels> A{};

    assert (h > 0L);

    T _h{h};

    T jump{2.0L};

    auto _fp = func(x0 + _h);
    auto _fm = func(x0 - _h);

    A[0][0] = 0.5L * (_fp - _fm) / _h;
    T deriv = A[0][0];

    auto err = std::numeric_limits<T>::max();

    // loop over columns
    for (int n = 1; n < nlevels; ++n) {
        _h /= jump;
        _fp = func(x0 + _h);
        _fm = func(x0 - _h);

        A[0][n] = 0.5L * (_fp - _fm) / _h;

        // for the current column n, fill in all the rows m > 0
        for (int m = 1; m <= n; ++m) {

            A[m][n] = (std::pow(2.0L * jump, m) * A[m-1][n] - A[m-1][n-1]) /
                (std::pow(2.0L * jump, m) - 1);

            // estimate the error so far, using the current row
            // and the estimates in the row above us
            T err_m = std::max(std::abs(A[m][n] - A[m-1][n]),
                               std::abs(A[m][n] - A[m-1][n-1]));

            if (err_m < err) {
                err = err_m;
                deriv = A[m][n];
            }
        }

        // at this point, A(n, n) should be the most accurate estimate
        // for the derivative, but it's possible that things have
        // become worse -- let's check
        if (std::abs(A[n][n] - A[n-1][n-1]) > 2.0 * err) {
            break;
        }

    }

    return {deriv, err};
}

template <typename T>
std::pair<T, T>
adaptive_diff2(std::function<T(T)> func, T x0, T h, int nlevels=10) {
    // This is like adaptive_diff, but for the second derivative.

    constexpr int max_levels = 50;

    assert (nlevels < max_levels);

    std::array<std::array<T, max_levels>, max_levels> A{};

    assert (h > 0L);

    T _h{h};

    T jump{2.0L};

    auto _fp = func(x0 + _h);
    auto _f0 = func(x0);
    auto _fm = func(x0 - _h);

    A[0][0] = (_fp - 2.0L * _f0 + _fm) / (_h * _h);
    T deriv = A[0][0];

    auto err = std::numeric_limits<T>::max();

    // loop over columns
    for (int n = 1; n < nlevels; ++n) {
        _h /= jump;
        _fp = func(x0 + _h);
        _f0 = func(x0);
        _fm = func(x0 - _h);

        A[0][n] = (_fp - 2.0L * _f0 + _fm) / (_h * _h);

        // for the current column n, fill in all the rows m > 0
        for (int m = 1; m <= n; ++m) {

            A[m][n] = (std::pow(2.0L * jump, m) * A[m-1][n] - A[m-1][n-1]) /
                (std::pow(2.0L * jump, m) - 1);

            // estimate the error so far, using the current row
            // and the estimates in the row above us
            T err_m = std::max(std::abs(A[m][n] - A[m-1][n]),
                               std::abs(A[m][n] - A[m-1][n-1]));

            if (err_m < err) {
                err = err_m;
                deriv = A[m][n];
            }
        }

        // at this point, A(n, n) should be the most accurate estimate
        // for the derivative, but it's possible that things have
        // become worse -- let's check
        if (std::abs(A[n][n] - A[n-1][n-1]) > 2.0 * err) {
            break;
        }

    }

    return {deriv, err};
}


#endif

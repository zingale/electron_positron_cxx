#ifndef BRENT_H
#define BRENT_H

#include <cassert>
#include <functional>
#include <cmath>
#include <limits>
#include <print>

// Brent's method for root finding (following the original Algol code
// in Brent 1971)

template <typename T>
T brent(std::function<T(T)> f, T a, T b, T tol=1.e-100) {

    // we will keep track of 3 points, a, b, and c, such that f(b) *
    // f(c) < 0 and |f(b)| <= |f(c)|.

    // a may be the same as c, b is the best approximation to the root,
    // a is the previous value of b, and the root must be between [b,
    // c]

    T fa = f(a);
    T fb = f(b);

    assert (fa * fb < 0);

    T c = a;
    T fc = fa;

    // d will be the step we want to take from b, and will be filled either by
    // inverse quadratic interpolation or bisection

    T d = b - a;

    // e is the old value of d
    T e = d;

    T machine_eps = std::numeric_limits<T>::epsilon();

    int iter{};
    while (true) {

        if (std::abs(fc) < std::abs(fb)) {
            // make |f(b)| the smallest
            a = b;
            fa = fb;

            b = c;
            fb = fc;

            c = a;
            fc = fa;
        }

        // root should be in [b, c].  a represents the previous value
        // of b, and is used to monitor if we are converging.

        T tol_act = 2 * machine_eps * std::abs(b) + tol;

        // m here is (currently) the increment that when added to b
        // would give bisection

        T m = 0.5 * (c - b);

        if (std::abs(m) <= tol_act || fb == 0.0L) {
            // we've converged
            return b;
        }

        T p{};
        T q{};

        if (std::abs(e) >= tol_act && std::abs(fa) > std::abs(fb)) {
            // attempt interpolation
            if (a == c) {
                // secant method
                T s = fb / fa;
                p = 2 * m * s;
                q = 1 - s;
            } else {
                // inverse quadratic interpolation
                T s = fb / fa;
                T r = fb / fc;
                T t = fa / fc;
                p = s * (2 * m * t * (t - r) - (b - a) * (r - 1));
                q = (t - 1) * (r - 1) * (s - 1);
            }

            if (p > 0) {
                q = -q;
            } else {
                p = -p;
            }

            // accept interpolation only if it falls within bounds
            if (2 * p < std::min(3 * m * q - std::abs(tol_act * q), std::abs(e * q))) {
                e = d;
                d = p / q;
            } else {
                // bisection
                d = m;
                e = m;
            }
        } else {
            // Bisection
            d = m;
            e = m;
        }

        a = b;
        fa = fb;
        if (std::abs(d) > tol_act) {
            b += d;
        } else {
            b += (m > 0 ? tol_act : -tol_act);
        }

        fb = f(b);
        if (fb * fc > 0) {
            c = a;
            fc = fa;
            e = d = b - a;
        }

        iter++;
    }
}


#endif

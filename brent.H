#ifndef BRENT_H
#define BRENT_H

#include <cassert>
#include <functional>
#include <cmath>
#include <limits>
#include <print>
#include <stdexcept>

#include "real_type.H"
#include "mp_math.H"

using namespace literals;

// Brent's method for root finding (following the original Algol code
// in Brent 1971)

template <typename T>
inline
T brent(const std::function<T(T)>& f, T a, T b, T atol=1.e-100_rt) {

    // we will keep track of 3 points, a, b, and c, such that f(b) *
    // f(c) < 0 and |f(b)| <= |f(c)|.

    // a may be the same as c, b is the best approximation to the root,
    // a is the previous value of b, and the root must be between [b,
    // c]

    T fa = f(a);
    T fb = f(b);

    if (fa * fb > 0.0_rt) {
        throw std::out_of_range("root is not contained in bounds");
    }

    T c = a;
    T fc = fa;

    // d will be the step we want to take from b, and will be filled either by
    // inverse quadratic interpolation or bisection

    T d = b - a;

    // e is the old value of d
    T e = d;

    while (true) {

        if (mp::abs(fc) < mp::abs(fb)) {
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

        // we create a tolerance that is a combination
        // of a relative tolerance (twice machine epsilon)
        // and an absolute tolerance
        T tol = 2 * std::numeric_limits<T>::epsilon() * mp::abs(b) + atol;

        // m here is (currently) the increment that when added to b
        // would give bisection

        T m = 0.5_rt * (c - b);

        if (mp::abs(m) <= tol || fb == 0.0_rt) {
            // we've converged
            return b;
        }

        T p{};
        T q{};

        if (mp::abs(e) >= tol && mp::abs(fa) > mp::abs(fb)) {
            // attempt interpolation
            if (a == c) {
                // secant method
                T s = fb / fa;
                p = 2.0_rt * m * s;
                q = 1.0_rt - s;
            } else {
                // inverse quadratic interpolation
                T s = fb / fa;
                T r = fb / fc;
                T t = fa / fc;
                p = s * (2.0_rt * m * t * (t - r) - (b - a) * (r - 1.0_rt));
                q = (t - 1.0_rt) * (r - 1.0_rt) * (s - 1.0_rt);
            }

            if (p > 0) {
                q = -q;
            } else {
                p = -p;
            }

            // accept interpolation only if it falls within bounds
            if (2.0_rt * p < std::min(3.0_rt * m * q - mp::abs(tol * q), mp::abs(e * q))) {
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
        if (mp::abs(d) > tol) {
            b += d;
        } else {
            b += (m > 0.0_rt ? tol : -tol);
        }

        fb = f(b);
        if (fb * fc > 0.0_rt) {
            c = a;
            fc = fa;
            e = d = b - a;
        }

    }
}


#endif

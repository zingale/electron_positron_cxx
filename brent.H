#ifndef BRENT_H
#define BRENT_H

#include <cassert>
#include <functional>
#include <cmath>
#include <limits>
#include <print>

// Brent's method for root finding (following the Wikipedia
// description)

template <typename T>
T brent(std::function<T(T)> f, T a, T b, T tol=8*std::numeric_limits<T>::epsilon()) {

    T fa = f(a);
    T fb = f(b);

    assert (fa * fb < 0);

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    T c = a;
    T fc = fa;

    bool mflag = true;

    T s{};
    T d{};

    while (fb != 0.0 and std::abs(b - a) > tol) {
        if (fa != fc && fb != fc) {
            // inverse quadratic interpolation
            s = a * fb * fc / (fa - fb) / (fa - fc) +
                b * fa * fc / (fb - fa) / (fb - fc) +
                c * fa * fb / (fc - fa) / (fc - fb);
        } else {
            // secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        bool check1 = s < 0.25 * (3 * a + b) || s > b;
        bool check2 = mflag and std::abs(s - b) >= 0.5 * std::abs(b - c);
        bool check3 = !mflag and std::abs(s - b) >= 0.5 * std::abs(c - d);
        bool check4 = mflag and std::abs(b - c) < tol;
        bool check5 = !mflag and std::abs(c - d) < tol;

        if (check1 || check2 || check3 || check4 || check5) {
            // do bisection
            s = 0.5 * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }

        T fs = f(s);
        std::println("s = {}, f(s) = {}", s, fs);

        d = c;
        c = b;

        if (fa * fs < 0.0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    return b;

}


#endif

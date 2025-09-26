#ifndef MP_MATH_H
#define MP_MATH_H

#include "real_type.H"
#include <cmath>
#if defined(USE_BOOST256)
#include <boost/multiprecision/cpp_bin_float.hpp>
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
#include <quadmath.h>
#endif

namespace mp {

    // power

    inline auto pow(const real_t& x, const real_t& y) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::pow(x, y);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return powq(x, y);
#else
        return std::pow(x, y);
#endif
    }

    inline auto sqrt(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::sqrt(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return sqrtq(x);
#else
        return std::sqrt(x);
#endif
    }

    // exponentiation

    inline auto exp(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::exp(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return expq(x);
#else
        return std::exp(x);
#endif
    }

    // log

    inline auto log(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::log(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return logq(x);
#else
        return std::log(x);
#endif
    }

    inline auto log10(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::log10(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return log10q(x);
#else
        return std::log10(x);
#endif
    }

    inline auto log1p(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::log1p(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return log1pq(x);
#else
        return std::log1p(x);
#endif
    }

    // trig

    inline auto sin(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::sin(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return sinq(x);
#else
        return std::sin(x);
#endif
    }

    inline auto cos(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::cos(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return cosq(x);
#else
        return std::cos(x);
#endif
    }

    inline auto cosh(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::cosh(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return coshq(x);
#else
        return std::cosh(x);
#endif
    }

    inline auto tanh(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::tanh(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return tanhq(x);
#else
        return std::tanh(x);
#endif
    }

    // basic math

    inline auto trunc(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::trunc(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return truncq(x);
#else
        return std::trunc(x);
#endif
    }

    inline auto abs(const real_t& x) -> real_t
    {
#if defined(USE_BOOST256)
        return boost::multiprecision::abs(x);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return fabsq(x);
#else
        return std::abs(x);
#endif
    }

    // these are approximate functions based on the singularity EOS /
    // Not-Quite-Transcendental Functions and their Applications
    // by Miller et al.

    // approximation to log_2(x)
    inline auto fastlg2(const real_t& x) -> real_t {
        // split into mant + exp
        int n{};
#if defined(USE_BOOST256)
        const real_t m = boost::multiprecision::frexp(x, &n);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        const real_t m = frexpq(x, &n);
#else
        const real_t m = std::frexp(x, &n);
#endif
        return 2.0_rt * (m - 1.0_rt) + static_cast<real_t>(n);
    }

    // approximation to 2**x
    inline auto fastpow2(const real_t& x) -> real_t {
        // create a mantissa and exp compatible with our log
#if defined(USE_BOOST256)
        const int flr = static_cast<int>(boost::multiprecision::floor(x));
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        const int flr = floorq(x);
#else
        const int flr = std::floor(x);
#endif
        const real_t rem = x - flr;
        const real_t mant = 0.5_rt * (rem + 1.0_rt);
        const int exponent = flr + 1;
#if defined(USE_BOOST256)
        return boost::multiprecision::ldexp(mant, exponent);
#elif defined(USE_FLOAT128) && (defined(__clang__) || defined(__aarch64__))
        return ldexpq(mant, exponent);
#else
        return std::ldexp(mant, exponent);
#endif
    }

}

#endif

#ifndef MP_MATH_H
#define MP_MATH_H

#include "real_type.H"
#include <cmath>
#if defined(USE_BOOST256)
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

namespace mp {

    // power

    inline real_t pow(const real_t& x, const real_t& y) {
#if defined(USE_BOOST256)
        return boost::multiprecision::pow(x, y);
#else
        return std::pow(x, y);
#endif
    }

    inline real_t sqrt(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::sqrt(x);
#else
        return std::sqrt(x);
#endif
    }

    // exponentiation

    inline real_t exp(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::exp(x);
#else
        return std::exp(x);
#endif
    }

    // log

    inline real_t log(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::log(x);
#else
        return std::log(x);
#endif
    }

    inline real_t log10(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::log10(x);
#else
        return std::log10(x);
#endif
    }

    inline real_t log1p(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::log1p(x);
#else
        return std::log1p(x);
#endif
    }

    // trig

    inline real_t sin(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::sin(x);
#else
        return std::sin(x);
#endif
    }

    inline real_t cos(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::cos(x);
#else
        return std::cos(x);
#endif
    }

    inline real_t cosh(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::cosh(x);
#else
        return std::cosh(x);
#endif
    }

    inline real_t tanh(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::tanh(x);
#else
        return std::tanh(x);
#endif
    }

    // basic math

    inline real_t trunc(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::trunc(x);
#else
        return std::trunc(x);
#endif
    }

    inline real_t abs(const real_t& x) {
#if defined(USE_BOOST256)
        return boost::multiprecision::abs(x);
#else
        return std::abs(x);
#endif
    }

}

#endif

#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#include <string>

#if defined(USE_FLOAT128)
#include <quadmath.h>
using real_t = __float128;
#elif defined(USE_LONG_DOUBLE)
using real_t = long double;
#else
using real_t = double;
#endif

inline namespace literals {

// literal for floating-point constants
// one is for 1.0_rt, and the other is for 1_rt

    constexpr real_t operator""_rt(long double x) {
        return static_cast<real_t>(x);
    }

    constexpr real_t operator""_rt(unsigned long long x) {
        return static_cast<real_t>(x);
    }
}

// std::numbers does not work with __float128

inline
real_t str_to_real_t(const char* s) {
#if defined(USE_FLOAT128)
    return strtoflt128(s, nullptr); 
#elif defined (USE_LONG_DOUBLE)
    return std::stold(s);
#else
    return std::stod(s);
#endif

}

namespace constants {
static inline const real_t pi = str_to_real_t("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068");

static inline const real_t sqrt2 = str_to_real_t("1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641573");
}

#endif

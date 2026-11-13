#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#include <string>

#ifdef USE_BOOST256
#include <boost/multiprecision/cpp_bin_float.hpp>
  using real_t = boost::multiprecision::cpp_bin_float_oct;

#elifdef USE_FLOAT128
#include <quadmath.h>
using real_t = __float128;

#elifdef USE_LONG_DOUBLE
using real_t = long double;

#else
using real_t = double;
#endif

// std::numbers does not work with __float128

inline auto str_to_real_t(const char* s) -> real_t
{
#ifdef USE_BOOST256
    return real_t(s);
#elifdef USE_FLOAT128
    return strtoflt128(s, nullptr);
#elifdef USE_LONG_DOUBLE
    return std::stold(s);
#else
    return std::stod(s);
#endif

}

inline namespace literals {

    // Parse wide floating-point literals directly in the selected precision.
    // The cooked operator below takes a long double, which would round the
    // literal before conversion to a 128- or 256-bit real_t.
#if defined(USE_BOOST256) || defined(USE_FLOAT128)
    // This numeric literal operator template is standard C++11.  The compiler
    // passes the original spelling as a character pack: for 1.25_rt, Digits
    // contains '1', '.', '2', '5'.  No floating-point conversion happens first.
    template <char... Digits>
    inline auto operator""_rt() -> real_t
    {
        // Turn the character pack into a null-terminated string.  The existing
        // converter uses strtoflt128 for 128 bits or Boost's string constructor
        // for 256 bits, so all supplied digits reach the selected backend.
        static constexpr char text[]{Digits..., '\0'};

        // Each spelling instantiates its own function and cached value, shared
        // by uses of that spelling.  Function-local static initialization is
        // thread-safe and parses only on first use, even inside EOS kernels.
        // Parsing is a runtime operation: 1.25_rt can initialize a const real_t,
        // but cannot be used in a constexpr initializer in these builds.
        static const real_t value = str_to_real_t(text);
        return value;
    }
#else
    constexpr auto operator""_rt(long double x) -> real_t
    {
        return static_cast<real_t>(x);
    }
#endif

    constexpr auto operator""_rt(unsigned long long x) -> real_t
    {
        return static_cast<real_t>(x);
    }
}

namespace constants {
static inline const real_t pi = str_to_real_t("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068");

static inline const real_t sqrt2 = str_to_real_t("1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641573");
}

#endif

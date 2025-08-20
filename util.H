#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <format>
#include <print>
#include <utility>

#include "real_type.H"

namespace util {

    template <typename T>
    T rel_error(T a, T b) {
        T err = std::abs(a - b) / std::abs(b);
        return err;
    }

    // wrappers around println
    // these add coloring and support for boost floats (by casting down)
    // note: this will lose some precision in output

#if defined(USE_BOOST256)
    // Trait to detect our big-float type
    template<typename T>
    struct is_boost_real_t : std::false_type {};

    template<>
    struct is_boost_real_t<real_t> : std::true_type {};

    // we can't print a boost 256-bit float, so cast it to __float128,
    // which does work with std::println()
    template<typename T>
    constexpr auto printable_cast(const T& x) {
        if constexpr (is_boost_real_t<T>::value) {
            return static_cast<long double>(x);
        } else {
            return x;
        }
    }

    // wrappers for println and format
    template<typename... Args>
    inline void println(const std::format_string<decltype(printable_cast(std::declval<Args>()))...> fmt,
                        Args&&... args)
    {
        std::println(fmt, printable_cast(std::forward<Args>(args))...);
    }

    template<typename... Args>
    inline std::string format(const std::format_string<decltype(printable_cast(std::declval<Args>()))...> fmt,
                              Args&&... args)
    {
        return std::format(fmt, printable_cast(std::forward<Args>(args))...);
    }

    template<typename... Args>
    inline void green_println(const std::format_string<decltype(printable_cast(std::declval<Args>()))...> fmt,
                              Args&&... args)
    {
        auto s = std::format(fmt, printable_cast(std::forward<Args>(args))...);
        std::println("\033[92m{}\033[0m", s);
    }


#else
    using std::println;
    using std::format;

    template<typename... Args>
    inline void green_println(const std::format_string<Args...> fmt, Args&&... args) {
        auto s = std::format(fmt, std::forward<Args>(args)...);
        std::println("\033[92m{}\033[0m", s);
    }
#endif



}
#endif

#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <format>
#include <print>
#include <utility>

namespace util {

template <typename T>
T rel_error(T a, T b) {
    T err = std::abs(a - b) / std::abs(b);
    return err;
}

// Helper to convert __float128 -> long double, leave everything else unchanged
template <typename T>
constexpr auto convert_for_print(T&& value) {
    if constexpr (std::is_same_v<std::remove_cvref_t<T>, __float128>) {
        return static_cast<long double>(value);
    } else {
        return std::forward<T>(value);
    }
}

// a println that converts __float128 to long double, since clang's println
// doesn't work with __float128
template <typename... Args>
void safe_println(std::format_string<decltype(convert_for_print(std::declval<Args>()))...> fmt,
                  Args&&... args)
{
    std::print("{}\n", std::format(fmt, convert_for_print(std::forward<Args>(args))...));
}

template <typename... Args>
void green_println(std::format_string<Args...> fmt, Args&&... args) {
    constexpr auto green = "\033[92m";
    constexpr auto reset = "\033[0m";
    std::print("{}{}{}\n", green,
               std::format(fmt, std::forward<Args>(args)...),
               reset);
}

}

#endif

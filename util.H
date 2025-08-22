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

template <typename... Args>
void green_println(std::format_string<Args...> fmt, Args&&... args) {
    constexpr auto green = "\033[92m";
    constexpr auto reset = "\033[0m";
    std::print("{}{}{}\n", green,
               std::format(fmt, std::forward<Args>(args)...),
               reset);
}

// print in red if we go above a threshold
template <typename... Args>
void threshold_println(real_t val, std::format_string<Args...> fmt, Args&&... args) {
    constexpr real_t threshold=1.e-10_rt;

    constexpr auto red = "\033[91m";
    constexpr auto reset = "\033[0m";
    if (std::abs(val) > threshold) {
        std::print("{}{}{}\n", red,
                   std::format(fmt, std::forward<Args>(args)...),
                   reset);
    } else {
        std::print("{}\n", std::format(fmt, std::forward<Args>(args)...));
    }
}

}
#endif

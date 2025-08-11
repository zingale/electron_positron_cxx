#ifndef UTIL_H
#define UTIL_H

#include <cmath>

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

}
#endif

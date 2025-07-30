#ifndef UTIL_H
#define UTIL_H

#include <cmath>

template <typename T>
T rel_error(T a, T b) {
    T err = std::abs(a - b) / std::abs(b);
    return err;
}

#endif

Precision
---------

Presently this uses 80-bit precision (via `long double`).

To get to 128-bit, I can try:

GCC/Clang `__float128` from `quadmath`:

```
#include <iostream>
#include <limits>
#include <quadmath.h>
#include <print>

int main() {
    __float128 a = 1.0Q;
    __float128 b = 3.0Q;

    __float128 result = a / b;

    // Convert to string for output
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.36Qg", result);

    std::cout << static_cast<double>(std::numeric_limits<__float128>::epsilon()) << std::endl;
    std::println("epsilon = {}", std::numeric_limits<__float128>::epsilon());

    std::cout << "a / b = " << buf << std::endl;

    return 0;
}
```

note that you can't sent a quad to `std::cout`, so you need to cast directly, but `std::print()` works.

This can be compiled via:

```
g++ -fext-numeric-literals -std=c++23 -o quad quad.cpp -lquadmath
```

another option is Boost.Multiprecision


#include <print>

#include "degeneracy_parameter_bounds.H"

int main() {

    auto [eta_min, eta_max] = bounds::get_eta_bounds<long double>(2.e4L, 2.e6L);

    std::println("eta bounds = [{}, {}]", eta_min, eta_max);

}


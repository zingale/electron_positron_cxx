#include <print>

#include "degeneracy_parameter_bounds.H"
#include "real_type.H"

using namespace literals;

int main() {

    auto [eta_min, eta_max] = bounds::get_eta_bounds<real_t>(2.e4_rt, 2.e6_rt);

    std::println("eta bounds = [{}, {}]", eta_min, eta_max);

}


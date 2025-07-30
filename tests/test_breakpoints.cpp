#include "fermi_integrals.H"

#include <print>

#include "real_type.H"

using namespace literals;

int main() {

    BreakPoints<real_t> b(0);

    auto [s1, s2, s3] = b.get_points(100.0);

    std::print("s1 = {}, s2 = {}, s3 = {}", s1, s2, s3);

}

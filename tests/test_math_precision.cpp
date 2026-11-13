#include <limits>
#include <stdexcept>

#include "fundamental_constants.H"
#include "mp_math.H"

// These checks use the selected type's epsilon so a narrowed intermediate
// cannot pass merely because it is accurate to double or long double precision.
auto main() -> int
{
    const auto check = [](real_t actual, real_t expected) {
        const real_t tolerance = 16_rt * std::numeric_limits<real_t>::epsilon() *
            mp::abs(expected);
        if (!(mp::abs(actual - expected) <= tolerance)) {
            throw std::runtime_error("math did not retain the selected precision");
        }
    };

    if (0.1_rt != str_to_real_t("0.1") ||
        1.234567890123456789012345678901234567890123456789012345678901234567890123456789_rt !=
        str_to_real_t("1.234567890123456789012345678901234567890123456789012345678901234567890123456789") ||
        1.e-100_rt != str_to_real_t("1.e-100")) {
        throw std::runtime_error("literal rounded before conversion to real_t");
    }
    check(C::k, str_to_real_t("1.3806490000000002e-16"));
    check(C::h, str_to_real_t("6.62607015e-27"));

    const real_t e = str_to_real_t("2.718281828459045235360287471352662497757247093699959574966967627724076630353547595");
    const real_t ln2 = str_to_real_t("0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156");
    const real_t log10_2 = str_to_real_t("0.3010299956639811952137388947244930267681898814621085413104274611271081892744245095");
    check(mp::pow(2_rt, 0.5_rt), constants::sqrt2);
    check(mp::sqrt(2_rt), constants::sqrt2);
    check(mp::exp(1_rt), e);
    check(mp::log(2_rt), ln2);
    check(mp::log10(2_rt), log10_2);
    check(mp::log1p(1_rt), ln2);
    check(mp::sin(constants::pi / 6_rt), 0.5_rt);
    check(mp::cos(constants::pi / 3_rt), 0.5_rt);
    check(mp::cosh(1_rt), (e + 1_rt / e) / 2_rt);
    check(mp::tanh(1_rt), (e * e - 1_rt) / (e * e + 1_rt));

    // Exercise tiny increments as well as ordinary transcendental results.
    const real_t delta = 8_rt * std::numeric_limits<real_t>::epsilon();
    check(mp::log1p(delta), delta);
    check(mp::abs(-delta), delta);
    check(mp::trunc(-1.75_rt), -1_rt);
    check(mp::fastpow2(mp::fastlg2(1.25_rt)), 1.25_rt);
}

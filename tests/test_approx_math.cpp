#include <limits>

#include "real_type.H"
#include "util.H"
#include "mp_math.H"

auto main() -> int
{


    const real_t x = str_to_real_t("2.45e10");
    const real_t y = str_to_real_t("5.42e-8");
    const real_t z = str_to_real_t("1.23456789");

    auto fastlg2x = mp::fastlg2(x);
    util::println("x = {}, fastlg2(x) = {}, fastexp(fastlg2(x)) = {}", x, fastlg2x, mp::fastpow2(fastlg2x));

    auto fastlg2y = mp::fastlg2(y);
    util::println("y = {}, fastlg2(y) = {}, fastexp(fastlg2(y)) = {}", y, fastlg2y, mp::fastpow2(fastlg2y));

    auto fastlg2z = mp::fastlg2(z);
    util::println("z = {}, fastlg2(z) = {}, fastexp(fastlg2(z)) = {}", z, fastlg2z, mp::fastpow2(fastlg2z));


    const real_t Tmin = 1.e4_rt;
    const real_t Tmax = 1.e11_rt;
    const int nT = 15;

    const real_t dlg2T = (mp::fastlg2(Tmax) - mp::fastlg2(Tmin)) / (static_cast<real_t>(nT - 1));

    for (int i = 0; i < nT; ++i) {
        util::println("{}", mp::fastpow2(mp::fastlg2(Tmin) + dlg2T * i));
    }

}


#include <limits>
#include <stdexcept>

#include "electron_positron.H"
#include "difference_utils.H"

auto main() -> int
{
    ElectronPositronEOS<real_t> eos;

    // CI condition 1: thermal electron derivatives are as small as the
    // positron contribution.  Reference values come from a 256-bit calculation
    // of the positron derivatives, combined using charge neutrality.
    const auto cold = eos.pe_state(1.e4_rt, 1.e7_rt, 0.5_rt);
    const auto check_relative = [](real_t actual, real_t expected, real_t tolerance,
                                   const char* message) {
        if (!(mp::abs(actual - expected) <= tolerance * mp::abs(expected))) {
            throw std::runtime_error(message);
        }
    };
    if constexpr (std::numeric_limits<real_t>::min_exponent10 < -505) {
        check_relative(cold.d2ne_drhodT,
                       -1.6063839e-500_rt, 1.e-7_rt, "incorrect electron rho-T derivative");
        check_relative(cold.d3ne_drho2dT,
                       9.988162326507118093246637986391132e-504_rt, 1.e-10_rt,
                       "incorrect electron rho-rho-T derivative");
        check_relative(cold.d3n_drho2dT,
                       1.9976324653014236186493275972782264e-503_rt, 1.e-10_rt,
                       "incorrect total rho-rho-T derivative");
    } else {
        // These pair contributions underflow in double precision.  The
        // derivatives must vanish rather than retain electron cancellation
        // noise from much larger intermediate values.
        if (cold.dne_dT != 0_rt || cold.d2ne_drhodT != 0_rt ||
            cold.d3ne_drho2dT != 0_rt || cold.d3n_drho2dT != 0_rt) {
            throw std::runtime_error("spurious derivative in an electron-dominated state");
        }
    }

    // With appreciable pairs, independently differentiate dn/dT twice in rho.
    // This checks the mixed derivative against nearby EOS states rather than
    // merely checking that electron and positron fields were copied correctly.
    const real_t rho = 1.e7_rt;
    const real_t temp = 5.e9_rt;
    const auto hot = eos.pe_state(rho, temp, 0.5_rt);
    const auto [diff, error] = fd::adaptive_diff2<real_t>(
        [&] (real_t r) -> real_t {
            ElectronPositronEOS<real_t, 1> first_order_eos;
            return first_order_eos.pe_state(r, temp, 0.5_rt).dn_dT;
        }, rho, 0.05_rt * rho, 8);
    check_relative(hot.d3n_drho2dT, diff, 1.e-6_rt,
                   "mixed number-density derivative disagrees with finite differences");
}

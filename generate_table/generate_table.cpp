#include <format>
#include <fstream>
#include <print>
#include <string>
#include <vector>

#include "real_type.H"
#include "electron_positron.H"
#include "helmholtz.H"
#include "util.H"

// create a table in the format of the Timmes & Swesty (2000) EOS.


#define USE_FAST_MATH 1

constexpr int rho_pts{841};  // density
constexpr int T_pts{321};  // temperature

#ifdef USE_FAST_MATH

const real_t rho_lo{mp::fastlg2(str_to_real_t("1.e-10"))};
const real_t rho_hi{mp::fastlg2(str_to_real_t("1.e11"))};

const real_t T_lo{mp::fastlg2(str_to_real_t("1.e3_rt"))};
const real_t T_hi{mp::fastlg2(str_to_real_t("1.e11_rt"))};

#else

const real_t rho_lo{-10.0e0_rt};
const real_t rho_hi{11.0e0_rt};

const real_t T_lo{3.0e0_rt};
const real_t T_hi{11.0e0_rt};

#endif

const real_t dlogrho = (rho_hi - rho_lo) / (static_cast<real_t>(rho_pts-1));
const real_t dlogT = (T_hi - T_lo) / (static_cast<real_t>(T_pts-1));


auto main() -> int
{

    std::string qnpts;
#if defined(QUAD20)
    qnpts = "20";
#elif defined(QUAD50)
    qnpts = "50";
#elif defined(QUAD100)
    qnpts = "100";
#elif defined(QUAD200)
    qnpts = "200";
#elif defined(QUAD400)
    qnpts = "400";
#elif defined(QUAD800)
    qnpts = "800";
#endif

    std::string precision;
#if defined(USE_BOOST256)
    precision = "256";
#elif defined(USE_FLOAT128)
    precision = "128";
#elif defined(USE_LONG_DOUBLE)
    precision = "80";
#else
    precision = "64";
#endif

    // we create the table using Ye = 1
    const real_t Ye{1.0_rt};

    // we'll compute everything first and then output

    std::vector<Helmholtz<real_t>> helm_v(T_pts * rho_pts);
    std::vector<EOSState<real_t>> eos_v(T_pts * rho_pts);

    // we will have OpenMP schedule things such that each thread is
    // working on similar temperatures, to better load balance

    #pragma omp parallel for schedule(static, 1)
    for (int j = 0; j < T_pts; ++j) {
#ifdef USE_FAST_MATH
        real_t T = mp::fastpow2(T_lo + static_cast<real_t>(j) * dlogT);
#else
        real_t T = std::pow(10.0_rt, T_lo + static_cast<real_t>(j) * dlogT);
#endif
        for (int i = 0; i < rho_pts; ++i) {
#ifdef USE_FAST_MATH
            real_t rho = mp::fastpow2(rho_lo + static_cast<real_t>(i) * dlogrho);
#else
            real_t rho = std::pow(10.0_rt, rho_lo + static_cast<real_t>(i) * dlogrho);
#endif

            // in the table, density varies fastest

            int index = j * rho_pts + i;

            util::println("rho = {}, T = {}", rho, T);

            auto [helm, eos] = get_helmholtz_terms<real_t>(rho, T, Ye);

            helm_v[index] = helm;
            eos_v[index] = eos;

        }

    }

#ifdef USE_FAST_MATH
    std::ofstream of(std::format("helm_table_fastmath_p{}_q{}.dat", precision, qnpts));
#else
    std::ofstream of(std::format("helm_table_p{}_q{}.dat", precision, qnpts));
#endif

    // there are 4 tables.  First is the free energy table

    for (auto & h : helm_v) {
        of << util::format("{:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}\n",
                           h.F, h.dF_drho, h.dF_dT,
                           h.d2F_drho2, h.d2F_dT2, h.d2F_drhodT,
                           h.d3F_drho2dT, h.d3F_drhodT2, h.d4F_drho2dT2);
    }

    // next is the pressure derivative table

    for (auto & e : eos_v) {
        of << util::format("{:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}\n",
                           e.dp_drho, e.d2p_drho2, e.d2p_drhodT, e.d3p_drho2dT);
    }

    // next is degeneracy parameter

    for (auto & e : eos_v) {
        of << util::format("{:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}\n",
                           e.eta, e.deta_drho, e.deta_dT, e.d2eta_drhodT);
    }

    // finally number density

    for (auto & e : eos_v) {
        of << util::format("{:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}\n",
                           e.n, e.dn_drho, e.dn_dT, e.d2n_drhodT);
    }

    // now some metadata
    of << "# generated with electron_positron_cxx\n";
#ifdef USE_FAST_MATH
    of << "# fast-math approximate log2 and pow2 used\n";
#endif
    of << util::format("# rho_lo = {}, rho_hi = {}, rho_pts = {}\n",
                       rho_lo, rho_hi, rho_pts);
    of << util::format("# T_lo = {}, T_hi = {}, T_pts = {}\n",
                       T_lo, T_hi, T_pts);
    of << util::format("# precision = {} bits\n", precision);
    of << util::format("# number of quadrature points = {}\n", qnpts);

}

#include <format>
#include <fstream>
#include <string>
#include <vector>

#include "real_type.H"
#include "electron_positron.H"
#include "helmholtz.H"
#include "util.H"

// create a table in the format of the Timmes & Swesty (2000) EOS.

constexpr int rho_pts{421};  // density
constexpr int T_pts{161};  // temperature

const real_t rho_lo{-10.0e0_rt};
const real_t rho_hi{11.0e0_rt};
const real_t dlogrho = (rho_hi - rho_lo) / (static_cast<real_t>(rho_pts-1));

const real_t T_lo{3.0e0_rt};
const real_t T_hi{11.0e0_rt};
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

    // we need one third derivative to compute via differencing
    std::vector<real_t> d3p_drho2dT_v(T_pts * rho_pts);

    // we will have OpenMP schedule things such that each thread is
    // working on similar temperatures, to better load balance

    #pragma omp parallel for schedule(static, 1)
    for (int j = 0; j < T_pts; ++j) {
        real_t T = std::pow(10.0_rt, T_lo + static_cast<real_t>(j) * dlogT);

        for (int i = 0; i < rho_pts; ++i) {
            real_t rho = std::pow(10.0_rt, rho_lo + static_cast<real_t>(i) * dlogrho);

            // in the table, density varies fastest

            int index = j * rho_pts + i;

            util::println("rho = {}, T = {}", rho, T);

            auto [helm, eos] = get_helmholtz_terms<real_t>(rho, T, Ye);

            helm_v[index] = helm;
            eos_v[index] = eos;

            // we need compute a derivative manually
            ElectronPositronEOS<real_t> ep_eos;

            // we'll do ∂/∂ρ(∂²p/∂ρ∂T)
            const real_t eps{0.01_rt};
            auto drho{eps * rho};
            auto [deriv, _err] = fd::adaptive_diff<real_t>([&] (real_t _rho) -> real_t
                {
                    auto es_eps = ep_eos.pe_state(_rho, T, Ye);
                    return es_eps.d2p_drhodT;
                }, rho, drho);

            d3p_drho2dT_v[index] = deriv;

        }

    }

    std::ofstream of(std::format("helm_table_p{}_q{}.dat", precision, qnpts));

    // there are 4 tables.  First is the free energy table

    for (auto & h : helm_v) {
        of << util::format("{:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}\n",
                           h.F, h.dF_drho, h.dF_dT,
                           h.d2F_drho2, h.d2F_dT2, h.d2F_drhodT,
                           h.d3F_drho2dT, h.d3F_drhodT2, h.d4F_drho2dT2);
    }

    // next is the pressure derivative table

    for (auto [e, d3p_drho2dT] : std::views::zip(eos_v, d3p_drho2dT_v)) {
        of << util::format("{:32.24g}  {:32.24g}  {:32.24g}  {:32.24g}\n",
                           e.dp_drho, e.d2p_drho2, e.d2p_drhodT, d3p_drho2dT);
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

}

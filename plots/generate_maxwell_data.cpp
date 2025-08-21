#include <print>
#include <vector>
#include <fstream>

#include "real_type.H"
#include "mp_math.H"
#include "electron_positron.H"
#include "maxwell_relations.H"
#include "util.H"

int main() {

    const real_t Ye{0.5_rt};

    const std::vector<real_t> Ts{1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11};

    std::vector<real_t> rhos;
    real_t rho_min = 1.e-4_rt;
    real_t rho_max = 1.e10_rt;
    int npts = 71;

    real_t dlogrho = (mp::log10(rho_max) - mp::log10(rho_min)) / static_cast<real_t>(npts-1);

    for (int i = 0; i < npts; ++i) {
        rhos.push_back(mp::pow(10.0_rt, mp::log10(rho_min) + i * dlogrho));
    }

    std::string qnpts;
#if defined(QUAD50)
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

    std::ofstream of(std::format("maxwell_p{}_npts{}.txt", precision, qnpts));

    for (auto T : Ts) {
        for (auto rho : rhos) {

            auto [scale1, error1] = maxwell_1<real_t>(rho, T, Ye);
            auto [scale2, error2] = maxwell_2<real_t>(rho, T, Ye);
            auto [scale3, error3] = maxwell_3<real_t>(rho, T, Ye);

            of << util::format("{:8.3g} {:8.3g} {:15.10g} {:15.10g} {:15.10g}\n",
                               T, rho, error1, error2, error3);
        }
    }

}

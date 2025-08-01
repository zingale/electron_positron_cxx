#ifndef ELECTRON_POSITRON_H
#define ELECTRON_POSITRON_H

#include <iostream>
#include <format>
#include <stdexcept>

#include "real_type.H"
#include "brent.H"
#include "degeneracy_parameter_bounds.H"
#include "fermi_integrals.H"
#include "fundamental_constants.H"

// our Fermi integrals will use a dimensionless temperature
constexpr real_t inv_compton_wavelength = C::m_e * C::c_light / C::h;
constexpr real_t rest_mass = C::m_e * C::c_light * C::c_light;

constexpr real_t coeff = 8 * pi * sqrt2 *
    inv_compton_wavelength * inv_compton_wavelength * inv_compton_wavelength;


template <typename T>
class EOSState {

public:
    T rho{};
    T temp{};
    T Y_e{};

    T n_e{};
    T n_pos{};

    T dne_drho{};
    T dne_dT{};

    T dnp_drho{};
    T dnp_dT{};

    T p_e{};
    T p_pos{};

    T dpe_drho{};
    T dpe_dT{};

    T dpp_drho{};
    T dpp_dT{};

    T e_e{};
    T e_pos{};

    T dee_drho{};
    T dee_dT{};

    T dep_drho{};
    T dep_dT{};

    T eta{};

    EOSState() = default;

};

template<typename T>
inline
std::ostream& operator<< (std::ostream& os, EOSState<T>& es) {

    std::string out{"\n"};

    out += std::format(" ρ = {:18.10g};       T = {:18.10g};      Yₑ = {:18.10g};\n", es.rho, es.temp, es.Y_e);
    out += std::format(" η = {:18.10g}\n", es.eta);
    out += std::format("n⁻ = {:18.10g};  ∂n⁻/∂ρ = {:18.10g};  ∂n⁻/∂T = {:18.10g};\n", es.n_e, es.dne_drho, es.dne_dT);
    out += std::format("n⁺ = {:18.10g};  ∂n⁺/∂ρ = {:18.10g};  ∂n⁺/∂T = {:18.10g};\n", es.n_pos, es.dnp_drho, es.dnp_dT);
    out += std::format("p⁻ = {:18.10g};  ∂p⁻/∂ρ = {:18.10g};  ∂p⁻/∂T = {:18.10g};\n", es.p_e, es.dpe_drho, es.dpe_dT);
    out += std::format("p⁺ = {:18.10g};  ∂p⁺/∂ρ = {:18.10g};  ∂p⁺/∂T = {:18.10g};\n", es.p_pos, es.dpp_drho, es.dpp_dT);
    out += std::format("e⁻ = {:18.10g};  ∂e⁻/∂ρ = {:18.10g};  ∂e⁻/∂T = {:18.10g};\n", es.e_e, es.dee_drho, es.dee_dT);
    out += std::format("e⁺ = {:18.10g};  ∂e⁺/∂ρ = {:18.10g};  ∂e⁺/∂T = {:18.10g};\n", es.e_pos, es.dep_drho, es.dep_dT);

    os << out;
    return os;
}

// for the root finding eta

template <typename T>
T n_e_constraint(T eta, T beta) {
    FermiIntegral<T> f12(0.5_rt, eta, beta);
    f12.evaluate(0);

    FermiIntegral<T> f32(1.5_rt, eta, beta);
    f32.evaluate(0);

    T n_e = coeff * beta * std::sqrt(beta) * (f12.F + beta * f32.F);
    return n_e;
}


template <typename T>
T n_p_constraint(T eta, T beta) {
    T eta_pos = -eta - 2.0 / beta;

    FermiIntegral<T> f12_pos(0.5_rt, eta_pos, beta);
    f12_pos.evaluate(0);

    FermiIntegral<T> f32_pos(1.5_rt, eta_pos, beta);
    f32_pos.evaluate(0);

    T n_pos = coeff * beta * std::sqrt(beta) * (f12_pos.F + beta * f32_pos.F);
    return n_pos;
}


template <typename T>
class ElectronPositronEOS {

    bool include_positrons{};

public:

    explicit ElectronPositronEOS(bool _include_positrons=true) :
        include_positrons{_include_positrons} {}

    EOSState<T> pe_state(T rho, T temp,  T Ye) {

        EOSState<T> es;

        es.rho = rho;
        es.temp = temp;
        es.Y_e = Ye;

        T n_e_net = Ye * C::N_A * rho;

        T beta = C::k * temp / rest_mass;

        T beta12 = std::sqrt(beta);
        T beta32 = beta * beta12;
        T beta52 = beta * beta32;

        // compute the degeneracy parameter
        auto [eta_min, eta_max] = bounds::get_eta_bounds(rho * Ye, temp, include_positrons);

        T eta{};
        bool with_positrons = include_positrons; // don't capture "this"

        try {
            eta = brent<T>([=] (T _eta) -> T
                {
                    auto n_e = n_e_constraint(_eta, beta);
                    T n_pos = (with_positrons) ? n_p_constraint(_eta, beta) : 0.0_rt;
                    return n_e_net - (n_e - n_pos);
                }, eta_min, eta_max);
        } catch (const std::out_of_range& e) {
            // try with larger bounds
            try {
                eta = brent<T>([=] (T _eta) -> T
                    {
                        auto n_e = n_e_constraint(_eta, beta);
                        T n_pos = (with_positrons) ? n_p_constraint(_eta, beta) : 0.0_rt;
                        return n_e_net - (n_e - n_pos);
                    }, -100.0_rt, 1.e7_rt);
            } catch (...) {
                std::cerr << "unable to solve for eta" << std::endl;
                abort();
            }
        } catch (...) {
            std::cerr << "unknown error doing eta solve" << std::endl;
            abort();
        }

        es.eta = eta;

        // for positrons
        T eta_pos = -eta - 2.0_rt / beta;

        // compute the number density, pressure and energy
        constexpr T pcoeff = coeff * (2.0 / 3.0) * rest_mass;
        constexpr T ecoeff = coeff * rest_mass;

        FermiIntegral<T> f12(0.5_rt, eta, beta);
        f12.evaluate(1);

        FermiIntegral<T> f32(1.5_rt, eta, beta);
        f32.evaluate(1);

        FermiIntegral<T> f52(2.5_rt, eta, beta);
        f52.evaluate(1);

        es.n_e = coeff * beta32 * (f12.F + beta * f32.F);
        es.p_e = pcoeff * beta52 * (f32.F + 0.5 * beta * f52.F);

        // this is (energy / volume, what we usually write as rho e)
        T E_e = ecoeff * beta52 * (f32.F + beta * f52.F);

        FermiIntegral<T> f12_pos(0.5_rt, eta_pos, beta);
        FermiIntegral<T> f32_pos(1.5_rt, eta_pos, beta);
        FermiIntegral<T> f52_pos(2.5_rt, eta_pos, beta);

        T E_pos{};
        if (include_positrons) {
            f12_pos.evaluate(1);
            f32_pos.evaluate(1);
            f52_pos.evaluate(1);

            es.n_pos = coeff * beta32 * (f12_pos.F + beta * f32_pos.F);
            es.p_pos = pcoeff * beta52 * (f32_pos.F + 0.5 * beta * f52_pos.F);
            E_pos = ecoeff * beta52 * (f32_pos.F + beta * f52_pos.F) + 2 * rest_mass * es.n_pos;
        }

        // compute the derivatives of eta and beta with respect to
        // density and temperature
        T dne_deta = coeff * beta32 * (f12.dF_deta + beta * f32.dF_deta);
        T dne_dbeta = 0.5 * coeff * beta12 * (3.0 * f12.F + 5.0 * beta * f32.F +
                                              2 * beta * (f12.dF_dbeta + beta * f32.dF_dbeta));

        T dnp_deta{};
        T dnp_dbeta{};
        if (include_positrons) {
            dnp_deta = coeff * beta32 * (-f12_pos.dF_deta - beta * f32_pos.dF_deta);
            dnp_dbeta = 0.5 * coeff / beta12 * (beta * (3.0 * f12_pos.F + 5.0 * beta * f32_pos.F) +
                                                2.0 * beta * beta * (f12_pos.dF_dbeta +
                                                                     beta * f32_pos.dF_dbeta) +
                                                4.0 * (f12_pos.dF_deta + beta * f32_pos.dF_deta));
        }

        //dbeta_drho = 0.0
        T dbeta_dT = C::k / rest_mass;

        T deta_drho = C::N_A * Ye / (dne_deta - dnp_deta);
        T deta_dT = -dbeta_dT * (dne_dbeta - dnp_dbeta) / (dne_deta - dnp_deta);

        // Compute partials of number density with density and temperature
        // For debugging
        es.dne_drho = dne_deta * deta_drho;
        es.dne_dT = dne_deta * deta_dT + dne_dbeta * dbeta_dT;

        if (include_positrons) {
            es.dnp_drho = dnp_deta * deta_drho;
            es.dnp_dT = dnp_deta * deta_dT + dnp_dbeta * dbeta_dT;
        }

        // Compute partials of pressure with density and temperature via the chain rule
        T dpe_deta = pcoeff * beta52 * (f32.dF_deta + 0.5 * beta * f52.dF_deta);
        T dpe_dbeta = 0.25 * pcoeff * beta32 * (10.0 * f32.F + 7.0 * beta * f52.F +
                                                4.0 * beta * (f32.dF_dbeta + 0.5 * beta * f52.dF_dbeta));

        es.dpe_drho = dpe_deta * deta_drho;
        es.dpe_dT = dpe_deta * deta_dT + dpe_dbeta * dbeta_dT;

        if (include_positrons) {
            T dpp_deta = -pcoeff * beta52 * (f32_pos.dF_deta + 0.5 * beta * f52_pos.dF_deta);
            T dpp_dbeta = pcoeff * beta12 * (beta * (2.5 * f32_pos.F + 1.75 * beta * f52_pos.F) +
                                             beta * beta * (f32_pos.dF_dbeta + 0.5 * beta * f52_pos.dF_dbeta) +
                                             2.0 * (f32_pos.dF_deta + 0.5 * beta * f52_pos.dF_deta));

            es.dpp_drho = dpp_deta * deta_drho;
            es.dpp_dT = dpp_deta * deta_dT + dpp_dbeta * dbeta_dT;
        }

        // Compute partials of energy with density and temperature
        T dEe_deta = ecoeff * beta52 * (f32.dF_deta + beta * f52.dF_deta);
        T dEe_dbeta = 0.5 * ecoeff * beta32 * (5.0 * f32.F + 7.0 * beta * f52.F +
                                               2.0 * beta * (f32.dF_dbeta + beta * f52.dF_dbeta));

        T dEe_drho = dEe_deta * deta_drho;   // dbeta_drho = 0
        T dEe_dT = dEe_deta * deta_dT + dEe_dbeta * dbeta_dT;

        T dEp_drho{};
        T dEp_dT{};
        if (include_positrons) {
            // these don't include the 2 m_e c**2 n_p part
            T dEp_deta = -ecoeff * beta52 * (f32_pos.dF_deta + beta * f52_pos.dF_deta);
            T dEp_dbeta = ecoeff * beta12 * (beta * (2.5 * f32_pos.F + 3.5 * beta * f52_pos.F) +
                                             beta * beta * (f32_pos.dF_dbeta + beta * f52_pos.dF_dbeta) +
                                             2.0 * (f32_pos.dF_deta + beta * f52_pos.dF_deta));

            dEp_drho = dEp_deta * deta_drho + 2.0 * rest_mass * es.dnp_drho;
            dEp_dT = dEp_deta * deta_dT + dEp_dbeta * dbeta_dT + 2.0 * rest_mass * es.dnp_dT;
        }

        // correct energy to be specific energy
        es.e_e = E_e / rho;
        es.dee_drho = (dEe_drho - E_e / rho) / rho;
        es.dee_dT = dEe_dT / rho;

        if (include_positrons) {
            es.e_pos = E_pos / rho;
            es.dep_drho = (dEp_drho - E_pos / rho) / rho;
            es.dep_dT = dEp_dT / rho;
        }

        return es;
    }

};

#endif

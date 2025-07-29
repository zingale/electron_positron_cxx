#ifndef ELECTRON_POSITRON_H
#define ELECTRON_POSITRON_H

#include <iostream>
#include <format>

#include "brent.H"
#include "degeneracy_parameter_bounds.H"
#include "fermi_integrals.H"
#include "fundamental_constants.H"

template <typename T>
class EOSState {

public:
    T n_e{};
    T n_pos{};

    T p_e{};
    T p_pos{};

    T e_e{};
    T e_pos{};

    T eta{};

    EOSState() {}

};

template<typename T>
inline
std::ostream& operator<< (std::ostream& os, EOSState<T>& es) {

    std::string out{};
    out += std::format("η  = {:20.10g}\n", es.eta);
    out += std::format("n⁻ = {:20.10g}\n", es.n_e);
    out += std::format("n⁺ = {:20.10g}\n", es.n_pos);
    out += std::format("p⁻ = {:20.10g}\n", es.p_e);
    out += std::format("p⁺ = {:20.10g}\n", es.p_pos);
    out += std::format("e⁻ = {:20.10g}\n", es.e_e);
    out += std::format("e⁺ = {:20.10g}\n", es.e_pos);

    os << out;
    return os;
}


template <typename T>
class ElectronPositronEOS {

    bool include_positrons{};

public:

    explicit ElectronPositronEOS(bool _include_positrons=true) :
        include_positrons{_include_positrons} {}

    EOSState<T> pe_state(T rho, T temp,  T Ye) {

        EOSState<T> es;

        T n_e_net = Ye * C::N_A * rho;

        // our Fermi integrals will use a dimensionless temperature
        constexpr T inv_compton_wavelength = C::m_e * C::c_light / C::h;
        constexpr T rest_mass = C::m_e * C::c_light * C::c_light;

        T beta = C::k * temp / rest_mass;

        T beta12 = std::sqrt(beta);
        T beta32 = beta * beta12;
        T beta52 = beta * beta32;

        constexpr T coeff = 8 * std::numbers::pi_v<T> * std::numbers::sqrt2_v<T> *
            inv_compton_wavelength * inv_compton_wavelength * inv_compton_wavelength;

        // compute the degeneracy parameter
        auto [eta_min, eta_max] = bounds::get_eta_bounds(rho * Ye, temp, include_positrons);

        T eta{};
        bool with_positrons = include_positrons; // don't capture "this"

        eta = brent<T>([=] (T _eta)
                       {
                           FermiIntegral<T> f12(0.5L, _eta, beta);
                           f12.evaluate(0);

                           FermiIntegral<T> f32(1.5L, _eta, beta);
                           f32.evaluate(0);

                           T n_e = coeff * beta32 * (f12.F + beta * f32.F);

                           T n_pos{};
                           if (with_positrons) {
                               T eta_pos = -_eta - 2.0 / beta;

                               FermiIntegral<T> f12_pos(0.5L, eta_pos, beta);
                               f12_pos.evaluate(0);

                               FermiIntegral<T> f32_pos(1.5L, eta_pos, beta);
                               f32_pos.evaluate(0);

                               n_pos = coeff * beta32 * (f12_pos.F + beta * f32_pos.F);
                           }
                           return n_e_net - (n_e - n_pos);
                       },
                       eta_min, eta_max);
        es.eta = eta;

        return es;
    }

};

#endif

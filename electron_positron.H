#ifndef ELECTRON_POSITRON_H
#define ELECTRON_POSITRON_H

#include <iostream>
#include <format>

#include "fermi_integrals.H"
#include "fundamental_constants.H"

template <typename T>
class EOSState {

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

        constexpr T coeff = 8 * std::numbers::pi_v<T> * std::numbers::sqrt2_v<T> *
            inv_compton_wavelength * inv_compton_wavelength * inv_compton_wavelength;

        return es;
    }

};

#endif

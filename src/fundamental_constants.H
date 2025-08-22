#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

#include "real_type.H"

using namespace literals;

namespace C
{
    // speed of light in vacuum
    const real_t c_light = 29979245800.0_rt;  // cm/s

    // boltzmann's constant
    const real_t k = 1.3806490000000002e-16_rt;  // erg/K

    // planck's constant over 2pi
    const real_t hbar = 1.0545718176461563e-27_rt;  // erg

    // planck's constant
    const real_t h = 6.62607015e-27_rt;  // erg s

    // avogradro's Number
    const real_t N_A = 6.02214076e+23_rt;  // mol^-1

    // mass of proton
    const real_t m_p = 1.67262192595e-24_rt;  // g

    // mass of neutron
    const real_t m_n = 1.6749275005600003e-24_rt;  // g

    // mass of electron
    const real_t m_e = 9.1093837139e-28_rt;  // g

    // atomic mass unit
    const real_t m_u = 1.66053906892e-24_rt; // g

    const real_t inv_compton_wavelength = m_e * c_light / h;
    const real_t rest_mass = m_e * c_light * c_light;

    const real_t dbeta_dT = k / rest_mass;
}

#endif

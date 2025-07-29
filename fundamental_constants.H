#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

namespace C
{
    // speed of light in vacuum
    constexpr long double c_light = 29979245800.0L;  // cm/s

    // boltzmann's constant
    constexpr long double k = 1.3806490000000002e-16L;  // erg/K

    // planck's constant over 2pi
    constexpr long double hbar = 1.0545718176461563e-27L;  // erg

    // planck's constant
    constexpr long double h = 6.62607015e-27L;  // erg s

    // avogradro's Number
    constexpr long double N_A = 6.02214076e+23L;  // mol^-1

    // mass of proton
    constexpr long double m_p = 1.67262192595e-24L;  // g

    // mass of neutron
    constexpr long double m_n = 1.6749275005600003e-24L;  // g

    // mass of electron
    constexpr long double m_e = 9.1093837139e-28L;  // g

    // atomic mass unit
    constexpr long double m_u = 1.66053906892e-24L; // g

}

#endif

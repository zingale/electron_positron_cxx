#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include "real_type.H"
#include "difference_utils.H"
#include "eos_types.H"
#include "electron_positron.H"

template <typename T>
inline auto get_helmholtz_terms(T rho, T temp, T Ye) -> Helmholtz<T>
{

    Helmholtz<T> helm;

    // get the thermodynamics for this state
    ElectronPositronEOS<T> eos;
    auto state = eos.pe_state(rho, temp, Ye);

    // fill in the derivatives we know
    helm.F = state.e - temp * state.s;


    helm.dF_drho = state.p / rho / rho;

    helm.d2F_dT2 = -state.ds_


    // for the 4th derivative term, we need to do a differencing

    return helm;

}

#endif

#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include "real_type.H"
#include "difference_utils.H"
#include "eos_types.H"
#include "electron_positron.H"

template <typename T>
inline auto get_helmholtz_terms(T rho, T temp, T Ye) 
    -> std::pair<Helmholtz<T>, EOSState<T>>
{

    Helmholtz<T> helm;
    helm.rho = rho;
    helm.temp = temp;
    helm.Y_e = Ye;

    // get the thermodynamics for this state
    ElectronPositronEOS<T> eos;
    auto state = eos.pe_state(rho, temp, Ye);

    // fill in the derivatives we know
    helm.F = state.e - temp * state.s;

    const T rho_inv = 1.0_rt / rho;

    helm.dF_drho = state.p * rho_inv * rho_inv;
    helm.dF_dT = -state.s;

    helm.d2F_drho2 = (state.dp_drho - 2.0_rt * rho_inv * state.p) * rho_inv * rho_inv;
    helm.d2F_dT2 = -state.ds_dT;
    helm.d2F_drhodT = state.dp_dT * rho_inv * rho_inv;

    helm.d3F_drho2dT = (state.d2p_drhodT - 2.0_rt * rho_inv * state.dp_dT) * rho_inv * rho_inv;
    helm.d3F_drhodT2 = -state.d2s_drhodT;

    helm.d4F_drho2dT2 = -state.d3s_drho2dT;

    return {helm, state};

}

#endif

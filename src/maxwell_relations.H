#ifndef MAXWELL_RELATIONS_H
#define MAXWELL_RELATIONS_H

#include <utility>

#include "real_type.H"
#include "electron_positron.H"

template <typename T>
std::pair<T, T>
inline
maxwell_1(T rho, T temp, T Ye) {

    ElectronPositronEOS<T> eos;
    auto es = eos.pe_state(rho, temp, Ye);

    const real_t ptot = es.p_e + es.p_pos;
    const real_t de_drho = es.dee_drho + es.dep_drho;
    const real_t dp_dT = es.dpe_dT + es.dpp_dT;

    const real_t term = mp::abs((ptot - (rho * rho * de_drho + temp * dp_dT)) / ptot);

    return {ptot, term};

}

template <typename T>
std::pair<T, T>
inline
maxwell_2(T rho, T temp, T Ye) {

    ElectronPositronEOS<T> eos;
    auto es = eos.pe_state(rho, temp, Ye);

    const real_t de_dT = es.dee_dT + es.dep_dT;
    const real_t ds_dT = es.dse_dT + es.dsp_dT;

    const real_t term = mp::abs((de_dT - temp * ds_dT) / de_dT);

    return {de_dT, term};


}

template <typename T>
std::pair<T, T>
inline
maxwell_3(T rho, T temp, T Ye) {

    ElectronPositronEOS<T> eos;
    auto es = eos.pe_state(rho, temp, Ye);

    const real_t ds_drho = es.dse_drho + es.dsp_drho;
    const real_t dp_dT = es.dpe_dT + es.dpp_dT;

    const real_t term = mp::abs((ds_drho + 1.0_rt / (rho * rho) * dp_dT) / ds_drho);

    return {ds_drho, term};
}

#endif

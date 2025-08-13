#ifndef ELECTRON_POSITRON_H
#define ELECTRON_POSITRON_H

#include <iostream>
#include <format>
#include <stdexcept>

#include "eos_types.H"
#include "real_type.H"
#include "brent.H"
#include "degeneracy_parameter_bounds.H"
#include "fermi_integrals.H"
#include "fundamental_constants.H"

// our Fermi integrals will use a dimensionless temperature

constexpr real_t coeff = 8.0_rt * pi * sqrt2 *
    C::inv_compton_wavelength * C::inv_compton_wavelength * C::inv_compton_wavelength;


// for the root finding eta

template <typename T>
inline
T n_e_constraint(T eta, T beta) {
    FermiIntegral<T> f12(0.5_rt, eta, beta);
    f12.evaluate(0);

    FermiIntegral<T> f32(1.5_rt, eta, beta);
    f32.evaluate(0);

    T n_e = coeff * beta * std::sqrt(beta) * (f12.F + beta * f32.F);
    return n_e;
}


template <typename T>
inline
T n_p_constraint(T eta, T beta) {
    T eta_tilde = -eta - 2.0_rt / beta;

    FermiIntegral<T> f12_pos(0.5_rt, eta_tilde, beta);
    f12_pos.evaluate(0);

    FermiIntegral<T> f32_pos(1.5_rt, eta_tilde, beta);
    f32_pos.evaluate(0);

    T n_pos = coeff * beta * std::sqrt(beta) * (f12_pos.F + beta * f32_pos.F);
    return n_pos;
}

template <typename T>
inline
T get_eta(T rhoYe, T temp) {
    // compute the degeneracy parameter

    T n_e_net = rhoYe * C::N_A;
    T beta = C::dbeta_dT * temp;

    auto [eta_min, eta_max] = bounds::get_eta_bounds(rhoYe, temp);

    T eta{};

    try {
        eta = brent<T>([=] (T _eta) -> T
            {
                auto n_e = n_e_constraint(_eta, beta);
                auto n_pos = n_p_constraint(_eta, beta);
                return n_e_net - (n_e - n_pos);
            }, eta_min, eta_max);
    } catch (const std::out_of_range& e) {
        // try with larger bounds
        try {
            eta = brent<T>([=] (T _eta) -> T
                {
                    auto n_e = n_e_constraint(_eta, beta);
                    auto n_pos = n_p_constraint(_eta, beta);
                    return n_e_net - (n_e - n_pos);
                }, -100.0_rt, 1.e7_rt);
        } catch (...) {
            std::cerr << "unable to solve for eta" << std::endl;
            abort();
        }
    }

    return eta;
}


template <typename T>
inline
std::pair<EtaBetaDerivs<T>, EtaBetaDerivs<T>>
get_n_derivs(const T beta,
             const FermiIntegral<T>& fx, const FermiIntegral<T>& fy,
             const FermiIntegral<T>& fx_pos, const FermiIntegral<T>& fy_pos) {
    // convenience wrapper for testing number density derivatives with respect
    // to η and β

    return get_eta_beta_derivs<T>(coeff, 1.0_rt, 3.0_rt, beta, fx, fy, fx_pos, fy_pos);
}


template <typename T>
inline
std::pair<EtaBetaDerivs<T>, EtaBetaDerivs<T>>
get_eta_beta_derivs(const T A, const T c, const T p, const T beta,
                    const FermiIntegral<T>& fx, const FermiIntegral<T>& fy,
                    const FermiIntegral<T>& fx_pos, const FermiIntegral<T>& fy_pos) {
    // we can write number density, pressure, and energy density for electrons as
    //
    // q = A β**(p/2) (F_x(η, β) + c β F_y(η, β))
    //
    // and for positrons
    //
    // w = A β**(p/2) (F_x(-η - 2/β, β) + c β F_y(-η - 2/β, β))

    // Here we compute the η and β first- and second-derivatives of
    // both quantities.  These are then used via the chain rule to
    // find the thermodynamic derivatives.

    EtaBetaDerivs<T> dq;
    EtaBetaDerivs<T> dw;

    T betap2 = std::pow(beta, p / 2.0_rt);
    T beta2 = beta * beta;
    T beta3 = beta * beta2;
    T beta4 = beta * beta3;
    T betai = 1.0_rt / beta;

    // first derivatives

    // ∂q/∂η
    dq.deta = A * betap2 * (c * beta * fy.dF_deta + fx.dF_deta);

    // ∂q/∂β
    dq.dbeta = A * betap2 * (0.5_rt * betai * (p * fx.F + c * (p + 2.0_rt) * beta * fy.F) +
                             (fx.dF_dbeta + c * beta * fy.dF_dbeta));

    // ∂w/∂η
    dw.deta = -A * betap2 * (c * beta * fy_pos.dF_deta + fx_pos.dF_deta);

    // ∂w/∂β
    dw.dbeta = 0.5_rt * A / beta2 * (betap2 * beta * (c * (p + 2.0_rt) * beta * fy_pos.F + p * fx_pos.F) +
                                     2.0_rt * betap2 * (beta2 * (fx_pos.dF_dbeta + c * beta * fy_pos.dF_dbeta) +
                                                        2.0_rt * (fx_pos.dF_deta + c * beta * fy_pos.dF_deta)));

    // second derivatives

    // ∂²q/∂η²
    dq.deta2 = A * betap2 * (c * beta * fy.d2F_deta2 + fx.d2F_deta2);

    // ∂²q/∂η∂β
    dq.detadbeta = A * betap2 * (c * beta * fy.d2F_detadbeta + c * fy.dF_deta +
                                 0.5_rt / beta * p * (c * beta * fy.dF_deta + fx.dF_deta) +
                                 fx.d2F_detadbeta);

    // ∂²q/∂β²
    dq.dbeta2 = A * betap2 * (c * beta * fy.d2F_dbeta2 + 2.0_rt * c * fy.dF_dbeta +
                              0.25_rt * p * (p - 2.0_rt) / beta2 * (c * beta * fy.F + fx.F) +
                              p / beta * (c * fy.F + c * beta * fy.dF_dbeta + fx.dF_dbeta) +
                              fx.d2F_dbeta2);

    // ∂²w/∂η²
    dw.deta2 = A * betap2 * (c * beta * fy_pos.d2F_deta2 + fx_pos.d2F_deta2);

    // ∂²w/∂η∂β
    dw.detadbeta = -A * betap2 * (2.0_rt * betai * betai * (fx_pos.d2F_deta2 + c * beta * fy_pos.d2F_deta2) +
                                  0.5_rt * betai * (c * (p + 2.0_rt) * beta * fy_pos.dF_deta + p * fx_pos.dF_deta) +
                                  (fx_pos.d2F_detadbeta + c * beta * fy_pos.d2F_detadbeta));

    // ∂²w/∂β²
    dw.dbeta2 = 0.25_rt * A * betap2 / beta4 * (p * (p - 2.0_rt) * beta2 * fx_pos.F + c * p * (p + 2.0_rt) * beta3 * fy_pos.F +
                                                8.0_rt * beta * ((p - 2.0_rt) * fx_pos.dF_deta + c * p * beta * fy_pos.dF_deta) +
                                                4.0_rt * beta4 * (fx_pos.d2F_dbeta2 + c * beta * fy_pos.d2F_dbeta2) +
                                                4.0_rt * beta3 * (p * fx_pos.dF_dbeta + c * (p + 2.0_rt) * beta * fy_pos.dF_dbeta) +
                                                16.0_rt * beta2 * (fx_pos.d2F_detadbeta + c * beta * fy_pos.d2F_detadbeta) +
                                                16.0_rt * (fx_pos.d2F_deta2 + c * beta * fy_pos.d2F_deta2));

    return {dq, dw};
}


template <typename T>
inline
RhoTDerivs<T>
get_eta_thermo_derivs(const T Ye,
                      const EtaBetaDerivs<T>& dn_e, const EtaBetaDerivs<T>& dn_pos) {
    // Compute the first- and second-derivatives of η with respect to
    // ρ and T.  This requires the derivatives of n⁻ and n⁺ with
    // respect to η and β.

    RhoTDerivs<T> deta;

    // ∂η/∂ρ and ∂η/∂T
    deta.drho = C::N_A * Ye / (dn_e.deta - dn_pos.deta);
    deta.dT = -C::dbeta_dT * (dn_e.dbeta - dn_pos.dbeta) / (dn_e.deta - dn_pos.deta);

    // ∂²η/∂T²
    deta.dT2 = (std::pow(C::dbeta_dT, 2.0_rt) * (dn_pos.dbeta2 - dn_e.dbeta2) +
                2.0_rt * C::dbeta_dT * deta.dT * (dn_pos.detadbeta - dn_e.detadbeta) +
                std::pow(deta.dT, 2.0_rt) * (dn_pos.deta2 - dn_e.deta2)) /
        (dn_e.deta - dn_pos.deta);

    // ∂²η/∂ρ∂T
    deta.drhodT = (C::dbeta_dT * (dn_pos.detadbeta - dn_e.detadbeta) +
                   deta.dT * (dn_pos.deta2 - dn_e.deta2)) /
        (dn_e.deta - dn_pos.deta);

    // ∂²η/∂ρ²
    deta.drho2 = deta.drho * deta.drho * (dn_pos.deta2 - dn_e.deta2) /
        (dn_e.deta - dn_pos.deta);

    return deta;
}


template <typename T>
inline
RhoTDerivs<T>
chain_rule(const EtaBetaDerivs<T>& df, const RhoTDerivs<T>& deta) {
    // Compute the thermodynamic derivatives of f via the chain rule
    // and the derivatives with respect to η and β

    RhoTDerivs<T> df_th;  // thermodynamic derivs

    // ∂f/∂ρ
    df_th.drho = deta.drho * df.deta;

    // ∂f/∂T
    df_th.dT = C::dbeta_dT * df.dbeta + deta.dT * df.deta;

    // ∂²f/∂ρ²
    df_th.drho2 = deta.drho * deta.drho * df.deta2 +
                  deta.drho2 * df.deta;

    // ∂²f/∂T²
    df_th.dT2 = C::dbeta_dT * C::dbeta_dT * df.dbeta2 +
                2.0_rt * C::dbeta_dT * deta.dT * df.detadbeta +
                deta.dT * deta.dT * df.deta2 +
                deta.dT2 * df.deta;

    // ∂²f/∂ρ∂T
    df_th.drhodT = C::dbeta_dT * deta.drho * df.detadbeta +
                   deta.dT * deta.drho * df.deta2 +
                   deta.drhodT * df.deta;

    return df_th;

}

template <typename T>
class ElectronPositronEOS {

public:

    ElectronPositronEOS() = default;

    EOSState<T> pe_state(T rho, T temp, T Ye) {

        EOSState<T> es;

        es.rho = rho;
        es.temp = temp;
        es.Y_e = Ye;

        T rho_inv{1.0_rt / rho};
        T temp_inv{1.0_rt / temp};

        T beta = temp * C::dbeta_dT;

        T beta12 = std::sqrt(beta);
        T beta32 = beta * beta12;
        T beta52 = beta * beta32;

        // get the degeneracy parameter for these thermodynamic
        // conditions
        auto eta = get_eta(rho * Ye, temp);

        es.eta = eta;
        es.beta = beta;

        // for positrons
        T eta_tilde = -eta - 2.0_rt / beta;
        T eta_pos = -eta;  // see C&G 24.334

        // compute the number density, pressure and energy
        constexpr T pcoeff = coeff * (2.0_rt / 3.0_rt) * C::rest_mass;
        constexpr T ecoeff = coeff * C::rest_mass;

        FermiIntegral<T> f12(0.5_rt, eta, beta);
        f12.evaluate(2);

        FermiIntegral<T> f32(1.5_rt, eta, beta);
        f32.evaluate(2);

        FermiIntegral<T> f52(2.5_rt, eta, beta);
        f52.evaluate(2);

        es.n_e = coeff * beta32 * (f12.F + beta * f32.F);
        es.p_e = pcoeff * beta52 * (f32.F + 0.5_rt * beta * f52.F);

        // this is (energy / volume, what we usually write as rho e)
        T E_e = ecoeff * beta52 * (f32.F + beta * f52.F);

        FermiIntegral<T> f12_pos(0.5_rt, eta_tilde, beta);
        f12_pos.evaluate(2);

        FermiIntegral<T> f32_pos(1.5_rt, eta_tilde, beta);
        f32_pos.evaluate(2);

        FermiIntegral<T> f52_pos(2.5_rt, eta_tilde, beta);
        f52_pos.evaluate(2);

        es.n_pos = coeff * beta32 * (f12_pos.F + beta * f32_pos.F);
        es.p_pos = pcoeff * beta52 * (f32_pos.F + 0.5_rt * beta * f52_pos.F);
        T E_pos = ecoeff * beta52 * (f32_pos.F + beta * f52_pos.F) + 2.0_rt * C::rest_mass * es.n_pos;

        // compute the derivatives of η and β with respect to
        // density and temperature

        // first we get the derivatives of number density with respect
        // to η and β
        const auto [dn_e, dn_pos] = get_eta_beta_derivs<T>(coeff, 1.0_rt, 3.0_rt,
                                                           beta, f12, f32, f12_pos, f32_pos);

        // note that ∂β/∂ρ = 0

        // now get the derivatives of η with respect to ρ and T
        auto deta = get_eta_thermo_derivs(Ye, dn_e, dn_pos);

        // store the derivatives
        es.deta_drho = deta.drho;
        es.deta_dT = deta.dT;

        es.d2eta_drho2 = deta.drho2;
        es.d2eta_drhodT = deta.drhodT;
        es.d2eta_dT2 = deta.dT2;

        // now we can compute the first and second derivatives of all
        // the thermodynamic quantities

        auto dne_thermo = chain_rule(dn_e, deta);
        es.dne_drho = dne_thermo.drho;
        es.dne_dT = dne_thermo.dT;
        es.d2ne_drho2 = dne_thermo.drho2;
        es.d2ne_dT2 = dne_thermo.dT2;
        es.d2ne_drhodT = dne_thermo.drhodT;

        auto dnp_thermo = chain_rule(dn_pos, deta);
        es.dnp_drho = dnp_thermo.drho;
        es.dnp_dT = dnp_thermo.dT;
        es.d2np_drho2 = dnp_thermo.drho2;
        es.d2np_dT2 = dnp_thermo.dT2;
        es.d2np_drhodT = dnp_thermo.drhodT;


        // Compute partials of pressure with density and temperature
        // via the chain rule.

        // First we get the derivatives of p with respect to η and β
        const auto [dp_e, dp_pos] = get_eta_beta_derivs<T>(pcoeff, 0.5_rt, 5.0_rt,
                                                           beta, f32, f52, f32_pos, f52_pos);

        auto dpe_thermo = chain_rule(dp_e, deta);
        es.dpe_drho = dpe_thermo.drho;
        es.dpe_dT = dpe_thermo.dT;
        es.d2pe_drho2 = dpe_thermo.drho2;
        es.d2pe_dT2 = dpe_thermo.dT2;
        es.d2pe_drhodT = dpe_thermo.drhodT;

        auto dpp_thermo = chain_rule(dp_pos, deta);
        es.dpp_drho = dpp_thermo.drho;
        es.dpp_dT = dpp_thermo.dT;
        es.d2pp_drho2 = dpp_thermo.drho2;
        es.d2pp_dT2 = dpp_thermo.dT2;
        es.d2pp_drhodT = dpp_thermo.drhodT;

        // Compute partials of energy density with density and
        // temperature via the chain rule.

        // First we get the derivatives of p with respect to η and β
        const auto [dE_e, dE_pos] = get_eta_beta_derivs<T>(ecoeff, 1.0_rt, 5.0_rt,
                                                           beta, f32, f52, f32_pos, f52_pos);

        auto dEe_thermo = chain_rule(dE_e, deta);
        auto dEp_thermo = chain_rule(dE_pos, deta);

        // add in rest mass term for positrons

        dEp_thermo.drho += 2.0_rt * C::rest_mass * es.dnp_drho;
        dEp_thermo.dT += 2.0_rt * C::rest_mass * es.dnp_dT;
        dEp_thermo.drho2 += 2.0_rt * C::rest_mass * es.d2np_drho2;
        dEp_thermo.dT2 += 2.0_rt * C::rest_mass * es.d2np_dT2;
        dEp_thermo.drhodT += 2.0_rt * C::rest_mass * es.d2np_drhodT;

        // correct energy to be specific energy
        es.e_e = E_e * rho_inv;
        es.dee_drho = (dEe_thermo.drho - E_e * rho_inv) * rho_inv;
        es.dee_dT = dEe_thermo.dT * rho_inv;
        es.d2ee_drho2 = (dEe_thermo.drho2 - 2.0_rt * dEe_thermo.drho * rho_inv -
                         2.0_rt * E_e * rho_inv * rho_inv) * rho_inv;
        es.d2ee_dT2 = dEe_thermo.dT2 * rho_inv;
        es.d2ee_drhodT = (dEe_thermo.drhodT - dEe_thermo.dT * rho_inv) * rho_inv;

        es.e_pos = E_pos * rho_inv;
        es.dep_drho = (dEp_thermo.drho - E_pos * rho_inv) * rho_inv;
        es.dep_dT = dEp_thermo.dT * rho_inv;
        es.d2ep_drho2 = (dEp_thermo.drho2 - 2.0_rt * dEp_thermo.drho * rho_inv -
                         2.0_rt * E_pos * rho_inv * rho_inv) * rho_inv;
        es.d2ep_dT2 = dEp_thermo.dT2 * rho_inv;
        es.d2ep_drhodT = (dEp_thermo.drhodT - dEp_thermo.dT * rho_inv) * rho_inv;

        // compute entropy -- this follows from Arnett's text, Eq. B.6
        // also C&G 24.368d
        T h_e = es.e_e + es.p_e * rho_inv;
        T dhe_dT = es.dee_dT + es.dpe_dT * rho_inv;
        T dhe_drho = es.dee_drho + es.dpe_drho * rho_inv - es.p_e * rho_inv * rho_inv;

        T h_pos = es.e_pos + es.p_pos * rho_inv;
        T dhp_dT = es.dep_dT + es.dpp_dT * rho_inv;
        T dhp_drho = es.dep_drho + es.dpp_drho * rho_inv - es.p_pos * rho_inv * rho_inv;

        es.s_e = h_e * temp_inv - es.n_e * C::k * es.eta * rho_inv;

        es.dse_dT = (dhe_dT - h_e * temp_inv) * temp_inv -
            C::k * (es.eta * es.dne_dT + es.n_e * deta.dT) * rho_inv;

        es.dse_drho = dhe_drho * temp_inv -
            C::k * (es.n_e * (deta.drho - es.eta * rho_inv) + es.eta * es.dne_drho ) * rho_inv;

        es.s_pos = h_pos * temp_inv - es.n_pos * C::k * eta_pos * rho_inv;

        es.dsp_dT = (dhp_dT - h_pos * temp_inv) * temp_inv -
            C::k * (eta_pos * es.dnp_dT - es.n_pos * deta.dT) * rho_inv;

        es.dsp_drho = dhp_drho * temp_inv -
            C::k * (es.n_pos * (-deta.drho - eta_pos * rho_inv) + eta_pos * es.dnp_drho) * rho_inv;

        return es;
    }

};

#endif

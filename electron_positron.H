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

// convenience wrapper for testing
template <typename T>
inline
std::pair<EtaBetaDerivs<T>, EtaBetaDerivs<T>>
get_n_derivs(const T beta,
             const FermiIntegral<T>& fx, const FermiIntegral<T>& fy,
             const FermiIntegral<T>& fx_pos, const FermiIntegral<T>& fy_pos) {

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
    // q = A β**(p/2) (F_x(η, β) + c beta F_y(η, β))
    //
    // and for positrons
    //
    // w = A β**(p/2) (F_x(-η - 2/β, β) + c beta F_y(-η - 2/β, β))

    // Here we compute the η and β first- and second-derivatives of
    // both quantities.  These are then used via the chain rule to
    // find the thermodynamic derivatives.

    EtaBetaDerivs<T> dq;
    EtaBetaDerivs<T> dw;

    T betap2 = std::pow(beta, p / 2.0_rt);
    T beta2 = beta * beta;
    T beta3 = beta * beta2;
    T beta4 = beta * beta3;
    T beta5 = beta * beta4;
    T betai = 1.0_rt / beta;

    // first derivatives

    // ∂q/∂η
    dq.deta = A * (c * beta * fy.dF_deta + fx.dF_deta) * betap2;

    // ∂q/∂β
    dq.dbeta = A * betap2 * (0.5_rt * p * (c * beta * fy.F + fx.F) * betai +
                             (c * fy.F + c * beta * fy.dF_dbeta + fx.dF_dbeta));

    // ∂w/∂η
    dw.deta = -A * (c * beta * fy_pos.dF_deta + fx_pos.dF_deta) * betap2;

    // ∂w/∂β
    dw.dbeta = 0.5 * A / beta4 * (p * (c * beta * fy_pos.F + fx_pos.F) * betap2 * beta3 +
                                  2.0_rt * betap2  * beta2 * (c * (beta2 * fy_pos.dF_dbeta + 2.0_rt * fy_pos.dF_deta) * beta +
                                                              (c * fy_pos.F + fx_pos.dF_dbeta) * beta2 +
                                                              2.0_rt * fx_pos.dF_deta));

    // second derivatives

    // ∂²n⁻/∂η²
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
    dw.detadbeta = -A * betap2 * (c * beta * (fy_pos.d2F_detadbeta + 2.0 * betai * betai * fy_pos.d2F_deta2) +
                                  c * fy_pos.dF_deta +
                                  0.5_rt * p * betai * (c * beta * fy_pos.dF_deta + fx_pos.dF_deta) +
                                  fx_pos.d2F_detadbeta +
                                  2.0 * betai * betai * fx_pos.d2F_deta2);

    // ∂²w/∂β²
    dw.dbeta2 = 0.25_rt * A * betap2 / beta4 * (c * p * (p + 2.0_rt) * beta3 * fy_pos.F +
                                                4.0_rt * c * p * beta4 * fy_pos.dF_dbeta + 8.0_rt * c * p * beta2 * fy_pos.dF_deta +
                                                4.0_rt * c * beta5 * fy_pos.d2F_dbeta2 + 8.0 * c * beta4 * fy_pos.dF_dbeta +
                                                16.0_rt * c * beta3 * fy_pos.d2F_detadbeta + 16 * c * beta * fy_pos.d2F_deta2 +
                                                p * (p - 2.0_rt) * beta2 * fx_pos.F +
                                                4.0 * p * beta3 * fx_pos.dF_dbeta + 8.0_rt * p * beta * fx_pos.dF_deta +
                                                4.0 * beta4 * fx_pos.d2F_dbeta2 + 16.0_rt * beta2 * fx_pos.d2F_detadbeta +
                                                -16.0_rt * beta * fx_pos.dF_deta + 16.0_rt * fx_pos.d2F_deta2);

    return {dq, dw};
}


template <typename T>
inline
RhoTDerivs<T>
get_eta_thermo_derivs(const T Ye,
                      const EtaBetaDerivs<T>& dn_e, const EtaBetaDerivs<T>& dn_pos) {

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
    // compute the thermodynamic derivatives of f via the chain rule
    // with η and β

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
                2.0_rt * C::dbeta_dT * deta.dT + df.detadbeta +
                deta.dT * deta.dT + df.deta2 +
                deta.dT2 * df.deta;

    // ∂²f/∂ρ∂T
    df_th.drhodT = C::dbeta_dT * deta.drho * df.detadbeta +
                   deta.dT * deta.drho * df.deta2 +
                   deta.drhodT * df.deta;

    return df_th;

}

template <typename T>
inline
RhoTDerivs<T>
get_eta_thermo_derivs(const T Ye,
                      const EtaBetaDerivs<T>& dn_e, const EtaBetaDerivs<T>& dn_pos) {

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
    // compute the thermodynamic derivatives of f via the chain rule
    // with η and β

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
                2.0_rt * C::dbeta_dT * deta.dT + df.detadbeta +
                deta.dT * deta.dT + df.deta2 +
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

        T n_e_net = Ye * C::N_A * rho;

        T beta = temp * C::dbeta_dT;

        T beta12 = std::sqrt(beta);
        T beta32 = beta * beta12;
        T beta52 = beta * beta32;

        // compute the degeneracy parameter
        auto [eta_min, eta_max] = bounds::get_eta_bounds(rho * Ye, temp);

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

        // compute the derivatives of eta and beta with respect to
        // density and temperature

        const auto [dn_e, dn_pos] = get_eta_beta_derivs<T>(coeff, 1.0_rt, 3.0_rt,
                                                           beta, f12, f32, f12_pos, f32_pos);

        // note that ∂β/∂ρ = 0

        auto deta = get_eta_thermo_derivs(Ye, dn_e, dn_pos);

        es.deta_drho = deta.drho;
        es.deta_dT = deta.dT;

        es.d2eta_drho2 = deta.drho2;
        es.d2eta_drhodT = deta.drhodT;
        es.d2eta_dT2 = deta.dT2;

        // Compute partials of number density with density and temperature
        es.dne_drho = dn_e.deta * deta.drho;
        es.dne_dT = dn_e.deta * deta.dT + dn_e.dbeta * C::dbeta_dT;

        es.dnp_drho = dn_pos.deta * deta.drho;
        es.dnp_dT = dn_pos.deta * deta.dT + dn_pos.dbeta * C::dbeta_dT;

        // Compute partials of pressure with density and temperature via the chain rule
        T dpe_deta = pcoeff * beta52 * (f32.dF_deta + 0.5_rt * beta * f52.dF_deta);
        T dpe_dbeta = 0.25_rt * pcoeff * beta32 * (10.0_rt * f32.F + 7.0_rt * beta * f52.F +
                                                   4.0_rt * beta * (f32.dF_dbeta + 0.5_rt * beta * f52.dF_dbeta));

        es.dpe_drho = dpe_deta * deta.drho;
        es.dpe_dT = dpe_deta * deta.dT + dpe_dbeta * C::dbeta_dT;

        T dpp_deta = -pcoeff * beta52 * (f32_pos.dF_deta + 0.5_rt * beta * f52_pos.dF_deta);
        T dpp_dbeta = pcoeff * beta12 * (beta * ((2.5_rt * f32_pos.F + 1.75_rt * beta * f52_pos.F) +
                                                 beta * (f32_pos.dF_dbeta + 0.5_rt * beta * f52_pos.dF_dbeta)) +
                                         2.0_rt * (f32_pos.dF_deta + 0.5_rt * beta * f52_pos.dF_deta));

        es.dpp_drho = dpp_deta * deta.drho;
        es.dpp_dT = dpp_deta * deta.dT + dpp_dbeta * C::dbeta_dT;

        // Compute partials of energy with density and temperature
        T dEe_deta = ecoeff * beta52 * (f32.dF_deta + beta * f52.dF_deta);
        T dEe_dbeta = 0.5 * ecoeff * beta32 * (5.0_rt * f32.F + 7.0_rt * beta * f52.F +
                                               2.0_rt * beta * (f32.dF_dbeta + beta * f52.dF_dbeta));

        T dEe_drho = dEe_deta * deta.drho;   // dbeta_drho = 0
        T dEe_dT = dEe_deta * deta.dT + dEe_dbeta * C::dbeta_dT;

        // these don't include the 2 m_e c**2 n_p part
        T dEp_deta = -ecoeff * beta52 * (f32_pos.dF_deta + beta * f52_pos.dF_deta);
        T dEp_dbeta = ecoeff * beta12 * (beta * ((2.5_rt * f32_pos.F + 3.5_rt * beta * f52_pos.F) +
                                                 beta * (f32_pos.dF_dbeta + beta * f52_pos.dF_dbeta)) +
                                         2.0_rt * (f32_pos.dF_deta + beta * f52_pos.dF_deta));

        T dEp_drho = dEp_deta * deta.drho + 2.0_rt * C::rest_mass * es.dnp_drho;
        T dEp_dT = dEp_deta * deta.dT + dEp_dbeta * C::dbeta_dT + 2.0_rt * C::rest_mass * es.dnp_dT;

        // correct energy to be specific energy
        es.e_e = E_e * rho_inv;
        es.dee_drho = (dEe_drho - E_e * rho_inv) * rho_inv;
        es.dee_dT = dEe_dT * rho_inv;

        es.e_pos = E_pos * rho_inv;
        es.dep_drho = (dEp_drho - E_pos * rho_inv) * rho_inv;
        es.dep_dT = dEp_dT * rho_inv;

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

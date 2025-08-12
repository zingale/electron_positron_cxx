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
std::pair<EtaBetaDerivs<T>, EtaBetaDerivs<T>>
get_n_derivs(T beta,
             const FermiIntegral<T>& f12, const FermiIntegral<T>& f32,
             const FermiIntegral<T>& f12_pos, const FermiIntegral<T>& f32_pos) {

    // number density derivatives needed for chain-rule computations
    EtaBetaDerivs<T> dn_e;
    EtaBetaDerivs<T> dn_pos;

    T beta12 = std::sqrt(beta);
    T beta32 = beta * beta12;
    T beta2 = beta * beta;
    T beta52 = beta * beta32;

    // ∂n⁻/∂η
    dn_e.deta = coeff * beta32 * (f12.dF_deta + beta * f32.dF_deta);

    // ∂n⁻/∂β
    dn_e.dbeta = 0.5_rt * coeff * beta12 * (3.0_rt * f12.F + 5.0_rt * beta * f32.F +
                                            2.0_rt * beta * (f12.dF_dbeta + beta * f32.dF_dbeta));

    // ∂n⁺/∂η
    dn_pos.deta = coeff * beta32 * (-f12_pos.dF_deta - beta * f32_pos.dF_deta);

    // ∂n⁺/∂β
    dn_pos.dbeta = 0.5_rt * coeff * (beta12 * (3.0_rt * f12_pos.F + 5.0_rt * beta * f32_pos.F +
                                               2.0_rt * beta * (f12_pos.dF_dbeta + beta * f32_pos.dF_dbeta)) +
                                     4.0_rt / beta12 * (f12_pos.dF_deta + beta * f32_pos.dF_deta));

    // now compute the second derivatives of number density

    // ∂²n⁻/∂η²
    dn_e.deta2 = coeff * beta32 * (f12.d2F_deta2 + beta * f32.d2F_deta2);

    // ∂²n⁻/∂η∂β
    dn_e.detadbeta = 0.5_rt * coeff * beta12 * (2.0_rt * beta * (f12.d2F_detadbeta + beta * f32.d2F_detadbeta) +
                                                3.0_rt * f12.dF_deta + 5.0_rt * beta * f32.dF_deta);

    // ∂²n⁻/∂β²
    dn_e.dbeta2 = 0.25_rt * coeff / beta12 * (3.0_rt * f12.F + 15.0_rt * beta * f32.F +
                                              4.0 * beta * (beta * (f12.d2F_dbeta2 + beta * f32.d2F_dbeta2) +
                                                            (3.0_rt * f12.dF_dbeta + 5.0_rt * beta * f32.dF_dbeta)));

    // ∂²n⁺/∂η²
    dn_pos.deta2 = coeff * beta32 * (f12_pos.d2F_deta2 + beta * f32_pos.d2F_deta2);

    // ∂²n⁺/∂η∂β
    dn_pos.detadbeta = -0.5_rt * coeff / beta12 * (2.0_rt * beta2 * (f12_pos.d2F_detadbeta + beta * f32_pos.d2F_detadbeta) +
                                                   beta * (3.0_rt * f12_pos.dF_deta + 5.0_rt * beta * f32_pos.dF_deta) +
                                                   4.0 * (f12_pos.d2F_deta2 + beta * f32_pos.d2F_deta2));

    // ∂²n⁺/∂β²
    dn_pos.dbeta2 = 0.25_rt * coeff / beta52 * (beta2 * ((3.0_rt * f12_pos.F + 15.0_rt * beta * f32_pos.F) +
                                                         4.0_rt * beta2 * (f12_pos.d2F_dbeta2 + beta * f32_pos.d2F_dbeta2) +
                                                         beta  * (12.0_rt * f12_pos.dF_dbeta + 20.0_rt * beta * f32_pos.dF_dbeta) +
                                                         16.0_rt * (f12_pos.d2F_detadbeta + beta * f32_pos.d2F_detadbeta)) +
                                                8.0_rt * beta * (f12_pos.dF_deta + 3.0_rt * beta * f32_pos.dF_deta) +
                                                16.0_rt * (f12_pos.d2F_deta2 + beta * f32_pos.d2F_deta2));

    return {dn_e, dn_pos};
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

        const auto [dn_e, dn_pos] = get_n_derivs<T>(beta, f12, f32, f12_pos, f32_pos);

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

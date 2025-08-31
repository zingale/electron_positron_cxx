#ifndef ELECTRON_POSITRON_H
#define ELECTRON_POSITRON_H

#include <iostream>
#include <format>
#include <stdexcept>

#include "eos_types.H"
#include "real_type.H"
#include "mp_math.H"
#include "brent.H"
#include "degeneracy_parameter_bounds.H"
#include "fermi_integrals.H"
#include "fundamental_constants.H"

// coefficient for number density expression in terms of Fermi
// integrals.

const real_t coeff = 8.0_rt * constants::pi * constants::sqrt2 *
    C::inv_compton_wavelength * C::inv_compton_wavelength * C::inv_compton_wavelength;


// for the root finding eta, we need to know the number density of
// electrons and positrons.  These will be constrained to the number
// density of electrons from charge neutrality.

template <typename T>
inline auto n_e_constraint(T eta, T beta) -> T
{
    FermiIntegral<T> f12(0.5_rt, eta, beta);
    f12.evaluate(0);

    FermiIntegral<T> f32(1.5_rt, eta, beta);
    f32.evaluate(0);

    T n_e = coeff * beta * mp::sqrt(beta) * (f12.F + beta * f32.F);
    return n_e;
}


template <typename T>
inline auto n_p_constraint(T eta, T beta) -> T
{
    T eta_tilde = -eta - 2.0_rt / beta;

    FermiIntegral<T> f12_pos(0.5_rt, eta_tilde, beta);
    f12_pos.evaluate(0);

    FermiIntegral<T> f32_pos(1.5_rt, eta_tilde, beta);
    f32_pos.evaluate(0);

    T n_pos = coeff * beta * mp::sqrt(beta) * (f12_pos.F + beta * f32_pos.F);
    return n_pos;
}

template <typename T>
inline auto get_eta(T rhoYe, T temp) -> T
{
    // compute the degeneracy parameter via root finding on
    // n_e_tot = n_e + n_pos, where n_e and n_pos are found
    // via integrating the Fermi distribution.

    T n_e_net = rhoYe * C::N_A;
    T beta = C::dbeta_dT * temp;

    // get a reasonable guess for the bounds of η

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
inline auto get_n_derivs(const T beta,
                         const FermiIntegral<T>& fx, const FermiIntegral<T>& fy,
                         const FermiIntegral<T>& fx_pos, const FermiIntegral<T>& fy_pos)
    -> std::pair<EtaBetaDerivs<T>, EtaBetaDerivs<T>>
{
    // convenience wrapper for testing number density derivatives with respect
    // to η and β

    return get_eta_beta_derivs<T>(coeff, 1.0_rt, 3.0_rt, beta, fx, fy, fx_pos, fy_pos);
}


template <typename T>
inline auto get_eta_beta_derivs(const T A, const T c, const T p, const T beta,
                                const FermiIntegral<T>& fx, const FermiIntegral<T>& fy,
                                const FermiIntegral<T>& fx_pos, const FermiIntegral<T>& fy_pos) 
    -> std::pair<EtaBetaDerivs<T>, EtaBetaDerivs<T>>
{
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

    T betap2 = mp::pow(beta, p / 2.0_rt);
    T beta2 = beta * beta;
    T beta3 = beta * beta2;
    T beta4 = beta * beta3;
    T beta5 = beta2 * beta3;
    T beta6 = beta2 * beta4;
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

    // third derivatives

    // ∂³q/∂η³
    dq.deta3 = A * betap2 * (fx.d3F_deta3 + c * beta * fy.d3F_deta3);

    // ∂³q/∂η²∂β
    dq.deta2dbeta = 0.5_rt * A * betap2 * betai * (p * (fx.d2F_deta2 + c * beta * fy.d2F_deta2) +
                                                   2.0 * beta * (fx.d3F_deta2dbeta + c * beta * fy.d3F_deta2dbeta + c * fy.d2F_deta2));

    // ∂³q/∂η∂β²
    dq.detadbeta2 = 0.25_rt * A * betap2 / beta2 *
        (p * (p - 2.0_rt) * (fx.dF_deta + c * beta * fy.dF_deta) +
         4.0_rt * p * beta * (fx.d2F_detadbeta + c * beta * fy.d2F_detadbeta + c * fy.dF_deta) +
         4.0_rt * beta2 * (fx.d3F_detadbeta2 + c * beta * fy.d3F_detadbeta2 + 2.0_rt * c * fy.d2F_detadbeta));

    // ∂³q/∂β³
    dq.dbeta3 = 0.125_rt * A * betap2 / beta3 *
        (6.0_rt * p * (p - 2.0_rt) * beta * (fx.dF_dbeta + c * beta * fy.dF_dbeta + c * fy.F) +
         p * (p * p - 6.0_rt * p + 8.0_rt) * (fx.F + c * beta * fy.F) +
         12.0_rt * p * beta2 * (fx.d2F_dbeta2 + c * beta * fy.d2F_dbeta2 + 2.0_rt * c * fy.dF_dbeta) +
         8.0_rt * beta3 * (fx.d3F_dbeta3 + c * beta * fy.d3F_dbeta3 + 3.0_rt * c * fy.d2F_dbeta2));

    // ∂³w/∂η³
    dw.deta3 = -A * betap2 * (fx_pos.d3F_deta3 + c * beta * fy_pos.d3F_deta3);

    // ∂³w/∂η²∂β
    dw.deta2dbeta = 0.5_rt * A * betap2 / beta2 *
        (beta * (2.0_rt * c * (beta2 * fy_pos.d3F_deta2dbeta + 2.0_rt * fy_pos.d3F_deta3) +
                 p * (fx_pos.d2F_deta2 + c * beta * fy_pos.d2F_deta2)) +
         2.0_rt * beta2 * (fx_pos.d3F_deta2dbeta + c * fy_pos.d2F_deta2) +
         4.0_rt * fx_pos.d3F_deta3);

    // ∂³w/∂η∂β²
    dw.detadbeta2 = 0.25 * A * betap2 / beta4 *
        (-c * p * p * beta3 * fy_pos.dF_deta - 4.0_rt * c * p * beta4 * fy_pos.d2F_detadbeta +
         -2.0_rt * c * p * beta3 * fy_pos.dF_deta - 8.0_rt * c * p * beta2 * fy_pos.d2F_deta2 +
         -4.0_rt * c * beta4 * (beta * fy_pos.d3F_detadbeta2 + 2.0 * fy_pos.d2F_detadbeta) +
         -16.0_rt * c * beta * (beta2 * fy_pos.d3F_deta2dbeta + fy_pos.d3F_deta3) +
         -p * p * beta2 * fx_pos.dF_deta - 4.0 * p * beta3 * fx_pos.d2F_detadbeta +
         2.0_rt * p * beta2 * fx_pos.dF_deta - 8.0_rt * p * beta * fx_pos.d2F_deta2 - 4.0_rt * beta4 * fx_pos.d3F_detadbeta2
         -16.0_rt * (beta2 * fx_pos.d3F_deta2dbeta - beta * fx_pos.d2F_deta2 + fx_pos.d3F_deta3));

    // ∂³w/∂β³
    dw.dbeta3 = 0.125_rt * A * betap2 / beta6 *
        (c * p * p * beta4 * (p * fy_pos.F + 6.0_rt * beta * fy_pos.dF_dbeta) +
         12.0_rt * c * p * p * beta3 * fy_pos.dF_deta - 4.0_rt * c * p * beta4 * fy_pos.F +
         12.0_rt * c * p * beta2 * (beta4 * fy_pos.d2F_dbeta2 + beta3 * fy_pos.dF_dbeta + 4.0_rt * beta2 * fy_pos.d2F_detadbeta - 2.0_rt * beta * fy_pos.dF_deta + 4.0_rt * fy_pos.d2F_deta2) +
         8.0_rt * c * beta * (beta6 * fy_pos.d3F_dbeta3 + 3.0_rt * beta5 * fy_pos.d2F_dbeta2  + 6.0_rt * beta4 * fy_pos.d3F_detadbeta2 +
                              12.0_rt * beta2 * fy_pos.d3F_deta2dbeta - 12.0_rt * beta * fy_pos.d2F_deta2 + 8.0_rt * fy_pos.d3F_deta3) +
         p * beta3 * (p * (p - 6.0_rt) * fx_pos.F + 6.0_rt * p * beta * fx_pos.dF_dbeta) +
         12.0_rt * p * p * beta2 * fx_pos.dF_deta + 8.0_rt * p * beta3 * fx_pos.F + 12.0_rt * p * beta5 * fx_pos.d2F_dbeta2 - 12.0_rt * p * beta4 * fx_pos.dF_dbeta +
         48.0_rt * p * beta3 * fx_pos.d2F_detadbeta - 72.0_rt * p * beta2 * fx_pos.dF_deta + 48.0_rt * p * beta * fx_pos.d2F_deta2 +
         8.0_rt * beta6 * fx_pos.d3F_dbeta3 + 48.0_rt * beta4 * fx_pos.d3F_detadbeta2 +
         96.0_rt * (-beta3 * fx_pos.d2F_detadbeta + beta2 * fx_pos.dF_deta + beta2 * fx_pos.d3F_deta2dbeta - 2.0_rt * beta * fx_pos.d2F_deta2) +
         64.0_rt * fx_pos.d3F_deta3);

    return {dq, dw};
}


template <typename T>
inline auto get_eta_thermo_derivs(const T Ye,
                                  const EtaBetaDerivs<T>& dn_e, const EtaBetaDerivs<T>& dn_pos)
    -> RhoTDerivs<T>
{
    // Compute the first- and second-derivatives of η with respect to
    // ρ and T.  This requires the derivatives of n⁻ and n⁺ with
    // respect to η and β.

    RhoTDerivs<T> deta;

    const real_t dbdT = C::dbeta_dT;
    const real_t dbdT2 = C::dbeta_dT * C::dbeta_dT;
    const real_t dbdT3 = dbdT2 * C::dbeta_dT;

    const real_t denomi = 1.0_rt / (dn_e.deta - dn_pos.deta);

    // first derivs

    // ∂η/∂ρ and ∂η/∂T
    deta.drho = C::N_A * Ye * denomi;
    deta.dT = -dbdT * (dn_e.dbeta - dn_pos.dbeta) * denomi;

    // second derivs

    // ∂²η/∂T²
    deta.dT2 = (dbdT2 * (dn_pos.dbeta2 - dn_e.dbeta2) +
                2.0_rt * dbdT * deta.dT * (dn_pos.detadbeta - dn_e.detadbeta) +
                mp::pow(deta.dT, 2.0_rt) * (dn_pos.deta2 - dn_e.deta2)) * denomi;


    // ∂²η/∂ρ∂T
    deta.drhodT = (dbdT * (dn_pos.detadbeta - dn_e.detadbeta) +
                   deta.dT * (dn_pos.deta2 - dn_e.deta2)) * deta.drho * denomi;

    // ∂²η/∂ρ²
    deta.drho2 = deta.drho * deta.drho * (dn_pos.deta2 - dn_e.deta2) * denomi;

    // third derivs

    // ∂³η/∂T³
    deta.dT3 = (dbdT3 * (dn_pos.dbeta3 - dn_e.dbeta3) +
                3.0_rt * dbdT2 * deta.dT * (dn_pos.detadbeta2 - dn_e.detadbeta2) +
                3.0_rt * C::dbeta_dT * mp::pow(deta.dT, 2.0_rt) * (dn_pos.deta2dbeta - dn_e.deta2dbeta) +
                3.0_rt * C::dbeta_dT * deta.dT2 * (dn_pos.detadbeta - dn_e.detadbeta) +
                mp::pow(deta.dT, 3.0_rt) * (dn_pos.deta3 - dn_e.deta3) +
                3.0_rt * deta.dT * deta.dT2 * (dn_pos.deta2 - dn_e.deta2)) * denomi;

    // ∂³η/∂ρ∂T²
    deta.drhodT2 = (dbdT2 * deta.drho * (dn_pos.detadbeta2 - dn_e.detadbeta2) +
                    2.0_rt * dbdT * deta.drho * deta.dT * (dn_pos.deta2dbeta - dn_e.deta2dbeta) +
                    2.0_rt * dbdT * deta.drhodT * (dn_pos.detadbeta - dn_e.detadbeta) +
                    mp::pow(deta.dT, 2.0_rt) * deta.drho * (dn_pos.deta3 - dn_e.deta3) +
                    2.0_rt * deta.dT * deta.drhodT * (dn_pos.deta2 - dn_e.deta2) +
                    deta.dT2 * deta.drho * (dn_pos.deta2 - dn_e.deta2)) * denomi;

    // ∂³η/∂ρ²∂T
    deta.drho2dT = (dbdT * mp::pow(deta.drho, 2.0_rt) * (dn_pos.deta2dbeta - dn_e.deta2dbeta) +
                    dbdT * deta.drho2 * (dn_pos.detadbeta - dn_e.detadbeta) +
                    deta.dT * mp::pow(deta.drho, 2.0_rt) * (dn_pos.deta3 - dn_e.deta3) +
                    deta.dT * deta.drho2 * (dn_pos.deta2 - dn_e.deta2) +
                    2.0_rt * deta.drho * deta.drhodT * (dn_pos.deta2 - dn_e.deta2)) * denomi;

    // ∂³η/∂ρ³
    deta.drho3 = (mp::pow(deta.drho, 2.0_rt) * (dn_pos.deta3 - dn_e.deta3) +
                  3.0_rt * deta.drho2 * (dn_pos.deta2 - dn_e.deta2)) * deta.drho * denomi;

    return deta;
}


template <typename T>
inline auto chain_rule(const EtaBetaDerivs<T>& df,
                       const RhoTDerivs<T>& deta)
    -> RhoTDerivs<T>
{
    // Compute the thermodynamic derivatives of f via the chain rule
    // and the derivatives with respect to η and β

    RhoTDerivs<T> df_th;  // thermodynamic derivs

    const real_t dbdT = C::dbeta_dT;
    const real_t dbdT2 = C::dbeta_dT * C::dbeta_dT;
    const real_t dbdT3 = dbdT2 * C::dbeta_dT;

    // first derivs

    // ∂f/∂ρ
    df_th.drho = deta.drho * df.deta;

    // ∂f/∂T
    df_th.dT = dbdT * df.dbeta + deta.dT * df.deta;

    // second derivs

    // ∂²f/∂ρ²
    df_th.drho2 = deta.drho * deta.drho * df.deta2 +
                  deta.drho2 * df.deta;

    // ∂²f/∂ρ∂T
    df_th.drhodT = dbdT * deta.drho * df.detadbeta +
                   deta.dT * deta.drho * df.deta2 +
                   deta.drhodT * df.deta;

    // ∂²f/∂T²
    df_th.dT2 = dbdT2 * df.dbeta2 +
                2.0_rt * dbdT * deta.dT * df.detadbeta +
                deta.dT * deta.dT * df.deta2 +
                deta.dT2 * df.deta;

    // ∂³f/∂ρ³
    df_th.drho3 = mp::pow(deta.drho, 3.0_rt) * df.deta3 +
                  3.0_rt * deta.drho * deta.drho2 * df.deta2 +
                  deta.drho3 * df.deta;

    // ∂³f/∂ρ²∂T
    df_th.drho2dT = dbdT * deta.drho * deta.drho * df.deta2dbeta +
                    dbdT * deta.drho2 * df.detadbeta +
                    deta.dT * deta.drho * deta.drho * df.deta3 +
                    deta.dT * deta.drho2 * df.deta2 +
                    2.0_rt * deta.drho * deta.drhodT * df.deta2 +
                    deta.drho2dT * df.deta;

    // ∂³f/∂ρ∂T²
    df_th.drhodT2 = dbdT2 * deta.drho * df.detadbeta2 +
                    2.0_rt * dbdT * deta.dT * deta.drho * df.deta2dbeta +
                    2.0_rt * dbdT * deta.drhodT * df.detadbeta +
                    deta.dT * deta.dT * deta.drho * df.deta3 +
                    2.0_rt * deta.dT * deta.drhodT * df.deta2 +
                    deta.drho * deta.dT2 * df.deta2 +
                    deta.drhodT2 * df.deta;
    // ∂³f/∂T³
    df_th.dT3 = dbdT3 * df.dbeta3 +
                3.0_rt * dbdT2 * deta.dT * df.detadbeta2 +
                3.0_rt * dbdT * deta.dT * deta.dT * df.deta2dbeta +
                3.0_rt * dbdT * deta.dT2 * df.detadbeta +
                deta.dT * deta.dT * deta.dT * df.deta3 +
                3.0_rt * deta.dT * deta.dT2 * df.deta2 +
                deta.dT3 * df.deta;

    return df_th;

}

template <typename T>
class ElectronPositronEOS {

public:

    ElectronPositronEOS() = default;

    auto pe_state(T rho, T temp, T Ye) -> EOSState<T>
    {

        EOSState<T> es;

        es.rho = rho;
        es.temp = temp;
        es.Y_e = Ye;

        T rho_inv{1.0_rt / rho};
        T temp_inv{1.0_rt / temp};

        T beta = temp * C::dbeta_dT;

        T beta12 = mp::sqrt(beta);
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

        T pcoeff = coeff * (2.0_rt / 3.0_rt) * C::rest_mass;
        T ecoeff = coeff * C::rest_mass;

        FermiIntegral<T> f12(0.5_rt, eta, beta);
        f12.evaluate(3);

        FermiIntegral<T> f32(1.5_rt, eta, beta);
        f32.evaluate(3);

        FermiIntegral<T> f52(2.5_rt, eta, beta);
        f52.evaluate(3);

        es.n_e = coeff * beta32 * (f12.F + beta * f32.F);
        es.p_e = pcoeff * beta52 * (f32.F + 0.5_rt * beta * f52.F);

        // E_e is (energy / volume, what we usually write as rho e)
        // We'll correct it to specific energy later.

        T E_e = ecoeff * beta52 * (f32.F + beta * f52.F);

        FermiIntegral<T> f12_pos(0.5_rt, eta_tilde, beta);
        f12_pos.evaluate(3);

        FermiIntegral<T> f32_pos(1.5_rt, eta_tilde, beta);
        f32_pos.evaluate(3);

        FermiIntegral<T> f52_pos(2.5_rt, eta_tilde, beta);
        f52_pos.evaluate(3);

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

        es.d3eta_drho3 = deta.drho3;
        es.d3eta_drho2dT = deta.drho2dT;
        es.d3eta_drhodT2 = deta.drhodT2;
        es.d3eta_dT3 = deta.dT3;

        // now we can compute the first and second derivatives of all
        // the thermodynamic quantities

        auto dne_thermo = chain_rule(dn_e, deta);
        es.dne_drho = dne_thermo.drho;
        es.dne_dT = dne_thermo.dT;
        es.d2ne_drho2 = dne_thermo.drho2;
        es.d2ne_drhodT = dne_thermo.drhodT;
        es.d2ne_dT2 = dne_thermo.dT2;
        es.d3ne_drho3 = dne_thermo.drho3;
        es.d3ne_drho2dT = dne_thermo.drho2dT;
        es.d3ne_drhodT2 = dne_thermo.drhodT2;
        es.d3ne_dT3 = dne_thermo.dT3;

        auto dnp_thermo = chain_rule(dn_pos, deta);
        es.dnp_drho = dnp_thermo.drho;
        es.dnp_dT = dnp_thermo.dT;
        es.d2np_drho2 = dnp_thermo.drho2;
        es.d2np_drhodT = dnp_thermo.drhodT;
        es.d2np_dT2 = dnp_thermo.dT2;
        es.d3np_drho3 = dnp_thermo.drho3;
        es.d3np_drho2dT = dnp_thermo.drho2dT;
        es.d3np_drhodT2 = dnp_thermo.drhodT2;
        es.d3np_dT3 = dnp_thermo.dT3;


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
        es.d2ee_drho2 = (dEe_thermo.drho2 - 2.0_rt * dEe_thermo.drho * rho_inv +
                         2.0_rt * E_e * rho_inv * rho_inv) * rho_inv;
        es.d2ee_dT2 = dEe_thermo.dT2 * rho_inv;
        es.d2ee_drhodT = (dEe_thermo.drhodT - dEe_thermo.dT * rho_inv) * rho_inv;

        es.e_pos = E_pos * rho_inv;
        es.dep_drho = (dEp_thermo.drho - E_pos * rho_inv) * rho_inv;
        es.dep_dT = dEp_thermo.dT * rho_inv;
        es.d2ep_drho2 = (dEp_thermo.drho2 - 2.0_rt * dEp_thermo.drho * rho_inv +
                         2.0_rt * E_pos * rho_inv * rho_inv) * rho_inv;
        es.d2ep_dT2 = dEp_thermo.dT2 * rho_inv;
        es.d2ep_drhodT = (dEp_thermo.drhodT - dEp_thermo.dT * rho_inv) * rho_inv;

        // compute entropy -- this follows from Arnett's text, Eq. B.6
        // also C&G 24.368d

        T h_e = es.e_e + es.p_e * rho_inv;
        T dhe_drho = es.dee_drho + es.dpe_drho * rho_inv - es.p_e * rho_inv * rho_inv;
        T dhe_dT = es.dee_dT + es.dpe_dT * rho_inv;
        T d2he_drho2 = es.d2ee_drho2 + (es.d2pe_drho2 + 2.0_rt * rho_inv * (-es.dpe_drho + rho_inv * es.p_e)) * rho_inv;
        T d2he_dT2 = es.d2ee_dT2 + rho_inv * es.d2pe_dT2;
        T d2he_drhodT = es.d2ee_drhodT + rho_inv * (es.d2pe_drhodT - rho_inv * es.dpe_dT);

        T h_pos = es.e_pos + es.p_pos * rho_inv;
        T dhp_drho = es.dep_drho + es.dpp_drho * rho_inv - es.p_pos * rho_inv * rho_inv;
        T dhp_dT = es.dep_dT + es.dpp_dT * rho_inv;
        T d2hp_drho2 = es.d2ep_drho2 + (es.d2pp_drho2 + 2.0_rt * rho_inv * (-es.dpp_drho + rho_inv * es.p_pos)) * rho_inv;
        T d2hp_dT2 = es.d2ep_dT2 + rho_inv * es.d2pp_dT2;
        T d2hp_drhodT = es.d2ep_drhodT + rho_inv * (es.d2pp_drhodT - rho_inv * es.dpp_dT);

        es.s_e = h_e * temp_inv - es.n_e * C::k * es.eta * rho_inv;

        es.dse_drho = dhe_drho * temp_inv -
            C::k * (es.n_e * (deta.drho - es.eta * rho_inv) + es.eta * es.dne_drho ) * rho_inv;

        es.dse_dT = (dhe_dT - h_e * temp_inv) * temp_inv -
            C::k * (es.eta * es.dne_dT + es.n_e * deta.dT) * rho_inv;

        es.d2se_drho2 = -C::k * rho_inv * (es.eta * es.d2ne_drho2 + es.n_e * es.d2eta_drho2 + 2.0_rt * es.deta_drho * es.dne_drho) +
            2.0_rt * C::k * rho_inv * rho_inv * (es.eta * es.dne_drho + es.n_e * es.deta_drho) +
            -2.0_rt * C::k * mp::pow(rho_inv, 3) * es.eta * es.n_e +
            d2he_drho2 * temp_inv;

        es.d2se_dT2 = -C::k * rho_inv * (es.eta * es.d2ne_dT2 + es.n_e * es.d2eta_dT2 + 2.0_rt * es.deta_dT * es.dne_dT) +
            temp_inv * (d2he_dT2 - 2.0_rt * temp_inv * (dhe_dT - h_e * temp_inv));

        es.d2se_drhodT = -C::k * rho_inv * (es.eta * es.d2ne_drhodT + es.n_e * es.d2eta_drhodT + es.deta_dT * es.dne_drho + es.deta_drho * es.dne_dT) +
            C::k * rho_inv * rho_inv * (es.eta * es.dne_dT + es.n_e * es.deta_dT)
            + temp_inv * (d2he_drhodT - dhe_drho * temp_inv);

        es.s_pos = h_pos * temp_inv - es.n_pos * C::k * eta_pos * rho_inv;

        es.dsp_drho = dhp_drho * temp_inv -
            C::k * (es.n_pos * (-deta.drho - eta_pos * rho_inv) + eta_pos * es.dnp_drho) * rho_inv;

        es.dsp_dT = (dhp_dT - h_pos * temp_inv) * temp_inv -
            C::k * (eta_pos * es.dnp_dT - es.n_pos * deta.dT) * rho_inv;

        es.d2sp_drho2 = C::k * rho_inv * (es.eta * es.d2np_drho2 + es.n_pos * es.d2eta_drho2 + 2.0_rt * es.deta_drho * es.dnp_drho) +
            -2.0_rt * C::k * rho_inv * rho_inv * (es.eta * es.dnp_drho + es.n_pos * es.deta_drho) +
            2.0_rt * C::k * mp::pow(rho_inv, 3) * es.eta * es.n_pos +
            d2hp_drho2 * temp_inv;

        es.d2sp_dT2 = C::k * rho_inv * (es.eta * es.d2np_dT2 + es.n_pos * es.d2eta_dT2 + 2.0_rt * es.deta_dT * es.dnp_dT) +
            temp_inv * (d2hp_dT2 - 2.0_rt * temp_inv * (dhp_dT - h_pos * temp_inv));

        es.d2sp_drhodT = C::k * rho_inv * (es.eta * es.d2np_drhodT + es.n_pos * es.d2eta_drhodT + es.deta_dT * es.dnp_drho + es.deta_drho * es.dnp_dT) +
            -C::k * rho_inv * rho_inv * (es.eta * es.dnp_dT + es.n_pos * es.deta_dT)
            + temp_inv * (d2hp_drhodT - dhp_drho * temp_inv);

        // now fill the total values of each of the thermodynamic quantities

        es.n = es.n_e + es.n_pos;
        es.dn_drho = es.dne_drho + es.dnp_drho;
        es.dn_dT = es.dne_dT + es.dnp_dT;
        es.d2n_drho2 = es.d2ne_drho2 + es.d2np_drho2;
        es.d2n_drhodT = es.d2ne_drhodT + es.d2np_drhodT;
        es.d2n_dT2 = es.d2ne_dT2 + es.d2np_dT2;
        es.d3n_drho3 = es.d3ne_drho3 + es.d3np_drho3;
        es.d3n_drho2dT = es.d3ne_drho2dT + es.d3np_drho2dT;
        es.d3n_drhodT2 = es.d3ne_drhodT2 + es.d3np_drhodT2;
        es.d3n_dT3 = es.d3ne_dT3 + es.d3np_dT3;

        es.p = es.p_e + es.p_pos;
        es.dp_drho = es.dpe_drho + es.dpp_drho;
        es.dp_dT = es.dpe_dT + es.dpp_dT;
        es.d2p_drho2 = es.d2pe_drho2 + es.d2pp_drho2;
        es.d2p_dT2 = es.d2pe_dT2 + es.d2pp_dT2;
        es.d2p_drhodT = es.d2pe_drhodT + es.d2pp_drhodT;

        es.e = es.e_e + es.e_pos;
        es.de_drho = es.dee_drho + es.dep_drho;
        es.de_dT = es.dee_dT + es.dep_dT;
        es.d2e_drho2 = es.d2ee_drho2 + es.d2ep_drho2;
        es.d2e_dT2 = es.d2ee_dT2 + es.d2ep_dT2;
        es.d2e_drhodT = es.d2ee_drhodT + es.d2ep_drhodT;

        es.s = es.s_e + es.s_pos;
        es.ds_drho = es.dse_drho + es.dsp_drho;
        es.ds_dT = es.dse_dT + es.dsp_dT;
        es.d2s_drho2 = es.d2se_drho2 + es.d2sp_drho2;
        es.d2s_dT2 = es.d2se_dT2 + es.d2sp_dT2;
        es.d2s_drhodT = es.d2se_drhodT + es.d2sp_drhodT;

        return es;
    }

};

#endif

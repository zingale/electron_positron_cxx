#include <print>

#include <array>

#include "real_type.H"
#include "electron_positron.H"
#include "eos_types.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"
#include "derivative_comparison.H"


auto main() -> int
{

    // test ∂p/∂ρ
    test_rho_deriv<1>(&EOSState<real_t>::dpe_drho, &EOSState<real_t>::p_e, "p⁻");
    test_rho_deriv<1>(&EOSState<real_t>::dpp_drho, &EOSState<real_t>::p_pos, "p⁺");
    test_rho_deriv<1>(&EOSState<real_t>::dp_drho, &EOSState<real_t>::p, "p");

    // test ∂p/∂T
    test_T_deriv<1>(&EOSState<real_t>::dpe_dT, &EOSState<real_t>::p_e, "p⁻");
    test_T_deriv<1>(&EOSState<real_t>::dpp_dT, &EOSState<real_t>::p_pos, "p⁺");
    test_T_deriv<1>(&EOSState<real_t>::dp_dT, &EOSState<real_t>::p, "p");

    // test ∂²p/∂ρ²
    test_rho2_deriv<2>(&EOSState<real_t>::d2pe_drho2, &EOSState<real_t>::p_e, "p⁻");
    test_rho2_deriv<2>(&EOSState<real_t>::d2pp_drho2, &EOSState<real_t>::p_pos, "p⁺");
    test_rho2_deriv<2>(&EOSState<real_t>::d2p_drho2, &EOSState<real_t>::p, "p");

    // test ∂²p/∂ρ∂T
    test_rhoT_deriv<2>(&EOSState<real_t>::d2pe_drhodT, &EOSState<real_t>::dpe_dT, "p⁻");
    test_rhoT_deriv<2>(&EOSState<real_t>::d2pp_drhodT, &EOSState<real_t>::dpp_dT, "p⁺");
    test_rhoT_deriv<2>(&EOSState<real_t>::d2p_drhodT, &EOSState<real_t>::dp_dT, "p");

    // test ∂²p/∂T²
    test_T2_deriv<2>(&EOSState<real_t>::d2pe_dT2, &EOSState<real_t>::p_e, "p⁻");
    test_T2_deriv<2>(&EOSState<real_t>::d2pp_dT2, &EOSState<real_t>::p_pos, "p⁺");
    test_T2_deriv<2>(&EOSState<real_t>::d2p_dT2, &EOSState<real_t>::p, "p");

    // test ∂³p/∂ρ³
    test_rho3_deriv(&EOSState<real_t>::d3pe_drho3, &EOSState<real_t>::dpe_drho, "p⁻");
    test_rho3_deriv(&EOSState<real_t>::d3pp_drho3, &EOSState<real_t>::dpp_drho, "p⁺");
    test_rho3_deriv(&EOSState<real_t>::d3p_drho3, &EOSState<real_t>::dp_drho, "p");

    // test ∂³p/∂ρ²∂T
    test_rho2T_deriv(&EOSState<real_t>::d3pe_drho2dT, &EOSState<real_t>::dpe_dT, "p⁻");
    test_rho2T_deriv(&EOSState<real_t>::d3pp_drho2dT, &EOSState<real_t>::dpp_dT, "p⁺");
    test_rho2T_deriv(&EOSState<real_t>::d3p_drho2dT, &EOSState<real_t>::dp_dT, "p");

    // test ∂³p/∂ρ∂T²
    test_rhoT2_deriv(&EOSState<real_t>::d3pe_drhodT2, &EOSState<real_t>::d2pe_dT2, "p⁻");
    test_rhoT2_deriv(&EOSState<real_t>::d3pp_drhodT2, &EOSState<real_t>::d2pp_dT2, "p⁺");
    test_rhoT2_deriv(&EOSState<real_t>::d3p_drhodT2, &EOSState<real_t>::d2p_dT2, "p");

    // test ∂³p/∂T³
    test_T3_deriv(&EOSState<real_t>::d3pe_dT3, &EOSState<real_t>::dpe_dT, "p⁻");
    test_T3_deriv(&EOSState<real_t>::d3pp_dT3, &EOSState<real_t>::dpp_dT, "p⁺");
    test_T3_deriv(&EOSState<real_t>::d3p_dT3, &EOSState<real_t>::dp_dT, "p");

}

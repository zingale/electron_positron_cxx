#include <print>

#include "real_type.H"
#include "electron_positron.H"
#include "eos_types.H"
#include "difference_utils.H"
#include "util.H"
#include "mp_math.H"
#include "derivative_comparison.H"


auto main() -> int
{

    // test ∂e/∂ρ
    test_rho_deriv(&EOSState<real_t>::dee_drho, &EOSState<real_t>::e_e, "e⁻");
    test_rho_deriv(&EOSState<real_t>::dep_drho, &EOSState<real_t>::e_pos, "e⁺");
    test_rho_deriv(&EOSState<real_t>::de_drho, &EOSState<real_t>::e, "e");

    // test ∂e/∂T
    test_T_deriv(&EOSState<real_t>::dee_dT, &EOSState<real_t>::e_e, "e⁻");
    test_T_deriv(&EOSState<real_t>::dep_dT, &EOSState<real_t>::e_pos, "e⁺");
    test_T_deriv(&EOSState<real_t>::de_dT, &EOSState<real_t>::e, "e");

    // test ∂²e/∂ρ²
    test_rho2_deriv(&EOSState<real_t>::d2ee_drho2, &EOSState<real_t>::e_e, "e⁻");
    test_rho2_deriv(&EOSState<real_t>::d2ep_drho2, &EOSState<real_t>::e_pos, "e⁺");
    test_rho2_deriv(&EOSState<real_t>::d2e_drho2, &EOSState<real_t>::e, "e");

    // test ∂²e/∂ρ∂T
    test_rhoT_deriv(&EOSState<real_t>::d2ee_drhodT, &EOSState<real_t>::dee_dT, "e⁻");
    test_rhoT_deriv(&EOSState<real_t>::d2ep_drhodT, &EOSState<real_t>::dep_dT, "e⁺");
    test_rhoT_deriv(&EOSState<real_t>::d2e_drhodT, &EOSState<real_t>::de_dT, "e");

    // test ∂²e/∂T²
    test_T2_deriv(&EOSState<real_t>::d2ee_dT2, &EOSState<real_t>::e_e, "e⁻");
    test_T2_deriv(&EOSState<real_t>::d2ep_dT2, &EOSState<real_t>::e_pos, "e⁺");
    test_T2_deriv(&EOSState<real_t>::d2e_dT2, &EOSState<real_t>::e, "e");

    // test ∂³e/∂ρ³
    test_rho3_deriv(&EOSState<real_t>::d3ee_drho3, &EOSState<real_t>::dee_drho, "e⁻");
    test_rho3_deriv(&EOSState<real_t>::d3ep_drho3, &EOSState<real_t>::dep_drho, "e⁺");
    test_rho3_deriv(&EOSState<real_t>::d3e_drho3, &EOSState<real_t>::de_drho, "e");

    // test ∂³e/∂ρ²∂T
    test_rho2T_deriv(&EOSState<real_t>::d3ee_drho2dT, &EOSState<real_t>::dee_dT, "e⁻");
    test_rho2T_deriv(&EOSState<real_t>::d3ep_drho2dT, &EOSState<real_t>::dep_dT, "e⁺");
    test_rho2T_deriv(&EOSState<real_t>::d3e_drho2dT, &EOSState<real_t>::de_dT, "e");

    // test ∂³e/∂ρ∂T²
    test_rhoT2_deriv(&EOSState<real_t>::d3ee_drhodT2, &EOSState<real_t>::d2ee_dT2, "e⁻");
    test_rhoT2_deriv(&EOSState<real_t>::d3ep_drhodT2, &EOSState<real_t>::d2ep_dT2, "e⁺");
    test_rhoT2_deriv(&EOSState<real_t>::d3e_drhodT2, &EOSState<real_t>::d2e_dT2, "e");

    // test ∂³e/∂T³
    test_T3_deriv(&EOSState<real_t>::d3ee_dT3, &EOSState<real_t>::dee_dT, "e⁻");
    test_T3_deriv(&EOSState<real_t>::d3ep_dT3, &EOSState<real_t>::dep_dT, "e⁺");
    test_T3_deriv(&EOSState<real_t>::d3e_dT3, &EOSState<real_t>::de_dT, "e");

}

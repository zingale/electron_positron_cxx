#include "real_type.H"
#include "eos_types.H"
#include "derivative_comparison.H"


auto main() -> int
{

    // test ∂s/∂ρ
    test_rho_deriv(&EOSState<real_t>::dse_drho, &EOSState<real_t>::s_e, "s⁻");
    test_rho_deriv(&EOSState<real_t>::dsp_drho, &EOSState<real_t>::s_pos, "s⁺");
    test_rho_deriv(&EOSState<real_t>::ds_drho, &EOSState<real_t>::s, "s");

    // test ∂s/∂T
    test_T_deriv(&EOSState<real_t>::dse_dT, &EOSState<real_t>::s_e, "s⁻");
    test_T_deriv(&EOSState<real_t>::dsp_dT, &EOSState<real_t>::s_pos, "s⁺");
    test_T_deriv(&EOSState<real_t>::ds_dT, &EOSState<real_t>::s, "s");

    // test ∂²s/∂ρ²
    test_rho2_deriv(&EOSState<real_t>::d2se_drho2, &EOSState<real_t>::s_e, "s⁻");
    test_rho2_deriv(&EOSState<real_t>::d2sp_drho2, &EOSState<real_t>::s_pos, "s⁺");
    test_rho2_deriv(&EOSState<real_t>::d2s_drho2, &EOSState<real_t>::s, "s");

    // test ∂²s/∂ρ∂T
    test_rhoT_deriv(&EOSState<real_t>::d2se_drhodT, &EOSState<real_t>::dse_dT, "s⁻");
    test_rhoT_deriv(&EOSState<real_t>::d2sp_drhodT, &EOSState<real_t>::dsp_dT, "s⁺");
    test_rhoT_deriv(&EOSState<real_t>::d2s_drhodT, &EOSState<real_t>::ds_dT, "s");

    // test ∂²s/∂T²
    test_T2_deriv(&EOSState<real_t>::d2se_dT2, &EOSState<real_t>::s_e, "s⁻");
    test_T2_deriv(&EOSState<real_t>::d2sp_dT2, &EOSState<real_t>::s_pos, "s⁺");
    test_T2_deriv(&EOSState<real_t>::d2s_dT2, &EOSState<real_t>::s, "s");

    // test ∂³s/∂ρ³
    test_rho3_deriv(&EOSState<real_t>::d3se_drho3, &EOSState<real_t>::dse_drho, "s⁻");
    test_rho3_deriv(&EOSState<real_t>::d3sp_drho3, &EOSState<real_t>::dsp_drho, "s⁺");
    test_rho3_deriv(&EOSState<real_t>::d3s_drho3, &EOSState<real_t>::ds_drho, "s");

    // test ∂³s/∂ρ²∂T
    test_rho2T_deriv(&EOSState<real_t>::d3se_drho2dT, &EOSState<real_t>::dse_dT, "s⁻");
    test_rho2T_deriv(&EOSState<real_t>::d3sp_drho2dT, &EOSState<real_t>::dsp_dT, "s⁺");
    test_rho2T_deriv(&EOSState<real_t>::d3s_drho2dT, &EOSState<real_t>::ds_dT, "s");

    // test ∂³s/∂ρ∂T²
    test_rhoT2_deriv(&EOSState<real_t>::d3se_drhodT2, &EOSState<real_t>::d2se_dT2, "s⁻");
    test_rhoT2_deriv(&EOSState<real_t>::d3sp_drhodT2, &EOSState<real_t>::d2sp_dT2, "s⁺");
    test_rhoT2_deriv(&EOSState<real_t>::d3s_drhodT2, &EOSState<real_t>::d2s_dT2, "s");

    // test ∂³s/∂T³
    test_T3_deriv(&EOSState<real_t>::d3se_dT3, &EOSState<real_t>::dse_dT, "s⁻");
    test_T3_deriv(&EOSState<real_t>::d3sp_dT3, &EOSState<real_t>::dsp_dT, "s⁺");
    test_T3_deriv(&EOSState<real_t>::d3s_dT3, &EOSState<real_t>::ds_dT, "s");

}

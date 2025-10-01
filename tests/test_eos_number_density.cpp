#include "real_type.H"
#include "eos_types.H"
#include "derivative_comparison.H"


auto main() -> int
{

    // test ∂n/∂ρ
    test_rho_deriv(&EOSState<real_t>::dne_drho, &EOSState<real_t>::n_e, "n⁻");
    test_rho_deriv(&EOSState<real_t>::dnp_drho, &EOSState<real_t>::n_pos, "n⁺");
    test_rho_deriv(&EOSState<real_t>::dn_drho, &EOSState<real_t>::n, "n");

    // test ∂n/∂T
    test_T_deriv(&EOSState<real_t>::dne_dT, &EOSState<real_t>::n_e, "n⁻");
    test_T_deriv(&EOSState<real_t>::dnp_dT, &EOSState<real_t>::n_pos, "n⁺");
    test_T_deriv(&EOSState<real_t>::dn_dT, &EOSState<real_t>::n, "n");

    // test ∂²n/∂ρ²
    test_rho2_deriv(&EOSState<real_t>::d2ne_drho2, &EOSState<real_t>::n_e, "n⁻");
    test_rho2_deriv(&EOSState<real_t>::d2np_drho2, &EOSState<real_t>::n_pos, "n⁺");
    test_rho2_deriv(&EOSState<real_t>::d2n_drho2, &EOSState<real_t>::n, "n");

    // test ∂²n/∂ρ∂T
    test_rhoT_deriv(&EOSState<real_t>::d2ne_drhodT, &EOSState<real_t>::dne_dT, "n⁻");
    test_rhoT_deriv(&EOSState<real_t>::d2np_drhodT, &EOSState<real_t>::dnp_dT, "n⁺");
    test_rhoT_deriv(&EOSState<real_t>::d2n_drhodT, &EOSState<real_t>::dn_dT, "n");

    // test ∂²n/∂T²
    test_T2_deriv(&EOSState<real_t>::d2ne_dT2, &EOSState<real_t>::n_e, "n⁻");
    test_T2_deriv(&EOSState<real_t>::d2np_dT2, &EOSState<real_t>::n_pos, "n⁺");
    test_T2_deriv(&EOSState<real_t>::d2n_dT2, &EOSState<real_t>::n, "n");

    // test ∂³n/∂ρ³
    test_rho3_deriv(&EOSState<real_t>::d3ne_drho3, &EOSState<real_t>::dne_drho, "n⁻");
    test_rho3_deriv(&EOSState<real_t>::d3np_drho3, &EOSState<real_t>::dnp_drho, "n⁺");
    test_rho3_deriv(&EOSState<real_t>::d3n_drho3, &EOSState<real_t>::dn_drho, "n");

    // test ∂³n/∂ρ²∂T
    test_rho2T_deriv(&EOSState<real_t>::d3ne_drho2dT, &EOSState<real_t>::dne_dT, "n⁻");
    test_rho2T_deriv(&EOSState<real_t>::d3np_drho2dT, &EOSState<real_t>::dnp_dT, "n⁺");
    test_rho2T_deriv(&EOSState<real_t>::d3n_drho2dT, &EOSState<real_t>::dn_dT, "n");

    // test ∂³n/∂ρ∂T²
    test_rhoT2_deriv(&EOSState<real_t>::d3ne_drhodT2, &EOSState<real_t>::d2ne_dT2, "n⁻");
    test_rhoT2_deriv(&EOSState<real_t>::d3np_drhodT2, &EOSState<real_t>::d2np_dT2, "n⁺");
    test_rhoT2_deriv(&EOSState<real_t>::d3n_drhodT2, &EOSState<real_t>::d2n_dT2, "n");

    // test ∂³n/∂T³
    test_T3_deriv(&EOSState<real_t>::d3ne_dT3, &EOSState<real_t>::dne_dT, "n⁻");
    test_T3_deriv(&EOSState<real_t>::d3np_dT3, &EOSState<real_t>::dnp_dT, "n⁺");
    test_T3_deriv(&EOSState<real_t>::d3n_dT3, &EOSState<real_t>::dn_dT, "n");

}

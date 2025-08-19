#ifndef EOS_TYPES_H
#define EOS_TYPES_H

#include "real_type.H"

template <typename T>
class EOSState {

public:
    T rho{};  // density (g/cm^3)
    T temp{};  // temperature (Kelvin)
    T Y_e{};  // electron fraction

    T n_e{};  // number density of electrons (1/cm^3)
    T n_pos{};  // number density of positrons (1/cm^3)

    T dne_drho{};
    T dne_dT{};
    T d2ne_drho2{};
    T d2ne_dT2{};
    T d2ne_drhodT{};

    T dnp_drho{};
    T dnp_dT{};
    T d2np_drho2{};
    T d2np_dT2{};
    T d2np_drhodT{};

    T p_e{};  // electron pressure (erg/cm^3)
    T p_pos{};  // positron pressure (erg/cm?^3)

    T dpe_drho{};
    T dpe_dT{};
    T d2pe_drho2{};
    T d2pe_dT2{};
    T d2pe_drhodT{};

    T dpp_drho{};
    T dpp_dT{};
    T d2pp_drho2{};
    T d2pp_dT2{};
    T d2pp_drhodT{};

    T e_e{};  // electron specific energy (erg/g)
    T e_pos{};  // positron specific energy (erg/g)

    T dee_drho{};
    T dee_dT{};
    T d2ee_drho2{};
    T d2ee_dT2{};
    T d2ee_drhodT{};

    T dep_drho{};
    T dep_dT{};
    T d2ep_drho2{};
    T d2ep_dT2{};
    T d2ep_drhodT{};

    T s_e{};  // electron specific entropy(erg/g/K)
    T s_pos{};  // positron specific entropy (erg/g/K)

    T dse_drho{};
    T dse_dT{};
    T d2se_drho2{};
    T d2se_dT2{};
    T d2se_drhodT{};


    T dsp_drho{};
    T dsp_dT{};
    T d2sp_drho2{};
    T d2sp_dT2{};
    T d2sp_drhodT{};

    T eta{};  // degeneracy parameter
    T beta{};  // k T / rest mass

    T deta_drho{};
    T deta_dT{};
    T d2eta_drho2{};
    T d2eta_drhodT{};
    T d2eta_dT2{};

    EOSState() = default;

};


template <typename T>
class EtaBetaDerivs {

public:
    T deta{};
    T dbeta{};

    T deta2{};
    T detadbeta{};
    T dbeta2{};

    EtaBetaDerivs() = default;

};


template <typename T>
class RhoTDerivs {

public:
    T drho{};
    T dT{};

    T drho2{};
    T drhodT{};
    T dT2{};

    RhoTDerivs() = default;

};


template<typename T>
inline
std::ostream& operator<< (std::ostream& os, EOSState<T>& es) {

    std::string out{"\n"};

    out += std::format(" ρ = {:15.8g}       T  = {:15.8g}      Yₑ  = {:15.8g}\n\n", es.rho, es.temp, es.Y_e);
    out += std::format(" β = {:15.8g}\n", es.beta);
    out += std::format(" η = {:15.8g}    ∂η/∂ρ  = {:15.8g}    ∂η/∂T  = {:15.8g}\n", es.eta, es.deta_drho, es.deta_dT);
    out += std::format("     {:15}   ∂²η/∂ρ² = {:15.8g}   ∂²η/∂T² = {:15.8g}   ∂²η/∂ρ∂T = {:15.8g}\n\n", " ", es.d2eta_drho2, es.d2eta_dT2, es.d2eta_drhodT);
    out += std::format("n⁻ = {:15.8g}   ∂n⁻/∂ρ  = {:15.8g}   ∂n⁻/∂T  = {:15.8g}\n", es.n_e, es.dne_drho, es.dne_dT);
    out += std::format("     {:15}  ∂²n⁻/∂ρ² = {:15.8g}  ∂²n⁻/∂T² = {:15.8g}  ∂²n⁻/∂ρ∂T = {:15.8g}\n", " ", es.d2ne_drho2, es.d2ne_dT2, es.d2ne_drhodT);
    out += std::format("n⁺ = {:15.8g}   ∂n⁺/∂ρ  = {:15.8g}   ∂n⁺/∂T  = {:15.8g}\n", es.n_pos, es.dnp_drho, es.dnp_dT);
    out += std::format("     {:15}  ∂²n⁺/∂ρ² = {:15.8g}  ∂²n⁺/∂T² = {:15.8g}  ∂²n⁺/∂ρ∂T = {:15.8g}\n\n", " ", es.d2np_drho2, es.d2np_dT2, es.d2np_drhodT);
    out += std::format("p⁻ = {:15.8g}   ∂p⁻/∂ρ  = {:15.8g}   ∂p⁻/∂T  = {:15.8g}\n", es.p_e, es.dpe_drho, es.dpe_dT);
    out += std::format("     {:15}  ∂²p⁻/∂ρ² = {:15.8g}  ∂²p⁻/∂T² = {:15.8g}  ∂²p⁻/∂ρ∂T = {:15.8g}\n", " ", es.d2pe_drho2, es.d2pe_dT2, es.d2pe_drhodT);
    out += std::format("p⁺ = {:15.8g}   ∂p⁺/∂ρ  = {:15.8g}   ∂p⁺/∂T  = {:15.8g}\n", es.p_pos, es.dpp_drho, es.dpp_dT);
    out += std::format("     {:15}  ∂²p⁺/∂ρ² = {:15.8g}  ∂²p⁺/∂T² = {:15.8g}  ∂²p⁺/∂ρ∂T = {:15.8g}\n\n", " ", es.d2pp_drho2, es.d2pp_dT2, es.d2pp_drhodT);
    out += std::format("e⁻ = {:15.8g}   ∂e⁻/∂ρ  = {:15.8g}   ∂e⁻/∂T  = {:15.8g}\n", es.e_e, es.dee_drho, es.dee_dT);
    out += std::format("     {:15}  ∂²e⁻/∂ρ² = {:15.8g}  ∂²e⁻/∂T² = {:15.8g}  ∂²e⁻/∂ρ∂T = {:15.8g}\n", " ", es.d2ee_drho2, es.d2ee_dT2, es.d2ee_drhodT);
    out += std::format("e⁺ = {:15.8g}   ∂e⁺/∂ρ  = {:15.8g}   ∂e⁺/∂T  = {:15.8g}\n", es.e_pos, es.dep_drho, es.dep_dT);
    out += std::format("     {:15}  ∂²e⁺/∂ρ² = {:15.8g}  ∂²e⁺/∂T² = {:15.8g}  ∂²e⁺/∂ρ∂T = {:15.8g}\n\n", " ", es.d2ep_drho2, es.d2ep_dT2, es.d2ep_drhodT);
    out += std::format("s⁻ = {:15.8g}   ∂s⁻/∂ρ  = {:15.8g}   ∂s⁻/∂T  = {:15.8g}\n", es.s_e, es.dse_drho, es.dse_dT);
    out += std::format("     {:15}  ∂²s⁻/∂ρ² = {:15.8g}  ∂²s⁻/∂T² = {:15.8g}  ∂²s⁻/∂ρ∂T = {:15.8g}\n", " ", es.d2se_drho2, es.d2se_dT2, es.d2se_drhodT);
    out += std::format("s⁺ = {:15.8g}   ∂s⁺/∂ρ  = {:15.8g}   ∂s⁺/∂T  = {:15.8g}\n", es.s_pos, es.dsp_drho, es.dsp_dT);
    out += std::format("     {:15}  ∂²s⁺/∂ρ² = {:15.8g}  ∂²s⁺/∂T² = {:15.8g}  ∂²s⁺/∂ρ∂T = {:15.8g}\n\n", " ", es.d2sp_drho2, es.d2sp_dT2, es.d2sp_drhodT);


    os << out;
    return os;
}


#endif

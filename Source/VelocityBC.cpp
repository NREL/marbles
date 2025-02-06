#include "VelocityBC.H"

namespace lbm::bc {
NoOp::NoOp() = default;

Constant::Constant()
{
    amrex::ParmParse pp("velocity_bc_constant");
    pp.query("dir", m_op.dir);
    pp.query("u0", m_op.u0);

    std::string m_model_type;
    pp.query("model_type", m_model_type); 

    if (m_model_type == "energyD3Q27") {
        m_op.m_model_type = 1;
        pp.query("Mach_0", m_op.Mach_0); 
        pp.query(
            "initial_temperature",
            m_op.initialTemperature); 
        pp.query(
            "adiabatic_exponent", m_op.adiabaticExponent); 
        pp.query("mean_molecular_mass", m_op.m_bar);        
        m_op.speedOfSound_Ref = std::sqrt(
            m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
            m_op.initialTemperature); 
        m_op.u0 = m_op.Mach_0 * m_op.speedOfSound_Ref;
    } else {
        // The defaults are used. No action required.
    }
}

Channel::Channel()
{
    amrex::ParmParse pp("velocity_bc_channel");
    pp.query("u_ref", m_op.u_ref);
}

Parabolic::Parabolic()
{
    amrex::ParmParse pp("velocity_bc_parabolic");
    pp.query("normal_dir", m_op.normal_dir);
    pp.query("tangential_dir", m_op.tangential_dir);
    pp.query("um", m_op.um);

    std::string m_model_type;
    pp.query("model_type", m_model_type); 
    
    if (m_model_type == "energyD3Q27") {
        m_op.m_model_type = 1;
        pp.query("Mach_m", m_op.Mach_m); 
        pp.query(
            "initial_temperature",
            m_op.initialTemperature); 
        pp.query(
            "adiabatic_exponent",
            m_op.adiabaticExponent); 
        pp.query("mean_molecular_mass", m_op.m_bar); 
        m_op.speedOfSound_Ref = std::sqrt(
            m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
            m_op.initialTemperature); 
        m_op.um = m_op.Mach_m * m_op.speedOfSound_Ref;
    } else {
        // The defaults are used. No action required.
    }
}

} // namespace lbm::bc

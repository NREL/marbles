#include "VelocityBC.H"

namespace lbm::bc {
NoOp::NoOp() = default;

Constant::Constant()
{
    amrex::ParmParse pp("velocity_bc_constant");
    pp.query("dir", m_op.dir);
    pp.query("u0", m_op.u0);

    pp.query("model_type", m_op.m_model_type); //ns: default is "isothermal". "energyD3Q27" activates product equilibrium, energy equation etc.
    
    if (m_op.m_model_type == "energyD3Q27")
    {
    pp.query("Mach_0", m_op.Mach_0); //ns:
    pp.query("initialTemperature", m_op.initialTemperature); //ns:initial condition temperature 
    //pp.query("adiabaticExponent", m_op.adiabaticExponent);   //ns: reference gamma. safety block. do not enable. not implemented.
    pp.query("meanMolecularMass", m_op.m_bar);               //ns: reference m_bar
    m_op.speedOfSound_Ref=std::sqrt(m_op.adiabaticExponent*(m_op.R_u/m_op.m_bar)*m_op.initialTemperature); //set the actual speed of sound
    m_op.u0 = m_op.Mach_0*m_op.speedOfSound_Ref;
    }
    else
    {
    //ns: The defaults are used. No action required.    
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
    
    pp.query("model_type", m_op.m_model_type); //ns: default is "isothermal". "energyD3Q27" activates product equilibrium, energy equation etc.
    
    if (m_op.m_model_type == "energyD3Q27")
    {
    pp.query("Mach_m", m_op.Mach_m); //ns:
    pp.query("initialTemperature", m_op.initialTemperature); //ns:initial condition temperature 
    //pp.query("adiabaticExponent", m_op.adiabaticExponent);   //ns: reference gamma. safety block. do not enable. not implemented.
    pp.query("meanMolecularMass", m_op.m_bar);               //ns: reference m_bar
    m_op.speedOfSound_Ref=std::sqrt(m_op.adiabaticExponent*(m_op.R_u/m_op.m_bar)*m_op.initialTemperature); //set the actual speed of sound
    m_op.um = m_op.Mach_m*m_op.speedOfSound_Ref;
    }
    else
    {
    //ns: The defaults are used. No action required.    
    }

}

} // namespace lbm::bc

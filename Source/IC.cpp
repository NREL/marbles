#include "IC.H"
namespace lbm::ic {

Constant::Constant()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    amrex::Vector<amrex::Real> MachComponents{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("MachComponents", MachComponents, 0, AMREX_SPACEDIM);
    for (int n = 0; n < MachComponents.size(); n++) {
        m_op.MachComponents[n] = MachComponents[n];
    }

    pp.query("initialTemperature", m_op.initialTemperature); //ns:initial condition temperature 
    //pp.query("adiabaticExponent", m_op.adiabaticExponent);   //ns: reference gamma. safety block. do not enable. not implemented.
    pp.query("meanMolecularMass", m_op.m_bar);               //ns: reference m_bar
    pp.query("model_type", m_op.m_model_type); //ns: default is "isothermal". "energyD3Q27" activates product equilibrium, energy equation etc.

    m_op.speedOfSound_Ref=std::sqrt(m_op.adiabaticExponent*(m_op.R_u/m_op.m_bar)*m_op.initialTemperature); //set the actual speed of sound

    if (m_op.m_model_type == "energyD3Q27")
    {
        for (int n = 0; n < MachComponents.size(); n++) m_op.velocity[n]=m_op.MachComponents[n]*m_op.speedOfSound_Ref;
    }

}

TaylorGreen::TaylorGreen()
{
    amrex::ParmParse pp(identifier());
    pp.query("rho0", m_op.rho0);
    pp.query("v0", m_op.v0);

    amrex::Vector<amrex::Real> omega{AMREX_D_DECL(1.0, 1.0, 1.0)};
    pp.queryarr("omega", omega, 0, AMREX_SPACEDIM);
    for (int n = 0; n < omega.size(); n++) {
        m_op.omega[n] = omega[n];
    }
}

viscosityTest::viscosityTest()
{
    // amrex::ParmParse pp(identifier());
    // pp.query("rho0", m_op.rho0);
    // pp.query("v0", m_op.v0);

    // amrex::Vector<amrex::Real> omega{AMREX_D_DECL(1.0, 1.0, 1.0)};
    // pp.queryarr("omega", omega, 0, AMREX_SPACEDIM);
    // for (int n = 0; n < omega.size(); n++) {
    //     m_op.omega[n] = omega[n];
    // }

    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    amrex::Vector<amrex::Real> MachComponents{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("MachComponents", MachComponents, 0, AMREX_SPACEDIM);
    for (int n = 0; n < MachComponents.size(); n++) {
        m_op.MachComponents[n] = MachComponents[n];
    }

    pp.query("initialTemperature", m_op.initialTemperature); //ns:initial condition temperature 
    //pp.query("adiabaticExponent", m_op.adiabaticExponent);   //ns: reference gamma. safety block. do not enable. not implemented.
    pp.query("meanMolecularMass", m_op.m_bar);               //ns: reference m_bar
    pp.query("model_type", m_op.m_model_type); //ns: default is "isothermal". "energyD3Q27" activates product equilibrium, energy equation etc.

    m_op.speedOfSound_Ref=std::sqrt(m_op.adiabaticExponent*(m_op.R_u/m_op.m_bar)*m_op.initialTemperature); //set the actual speed of sound

    if (m_op.m_model_type == "energyD3Q27")
    {
        for (int n = 0; n < MachComponents.size(); n++) m_op.velocity[n]=m_op.MachComponents[n]*m_op.speedOfSound_Ref;
    }

    pp.query("waveLength", m_op.waveLength); //ns:wavelength of the perturbation
}

} // namespace lbm::ic

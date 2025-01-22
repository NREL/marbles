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

    std::string m_model_type;
    pp.query("model_type", m_model_type); // ns: default is "isothermal".
                                          // "energyD3Q27" activates product
                                          // equilibrium, energy equation etc.

    if (m_model_type == "energyD3Q27") {
        m_op.m_model_type = 1;
        amrex::Vector<amrex::Real> MachComponents{AMREX_D_DECL(0, 0, 0)};
        pp.queryarr("MachComponents", MachComponents, 0, AMREX_SPACEDIM);
        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.MachComponents[n] = MachComponents[n];
        }

        pp.query(
            "initialTemperature",
            m_op.initialTemperature); // ns:initial condition temperature
        pp.query(
            "adiabaticExponent",
            m_op.adiabaticExponent);               // ns: reference gamma.
        pp.query("meanMolecularMass", m_op.m_bar); // ns: reference m_bar
        m_op.speedOfSound_Ref = std::sqrt(
            m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
            m_op.initialTemperature); // set the actual speed of sound

        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
        }
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
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("waveLength", m_op.waveLength); // ns:wavelength of the
                                             // perturbation
    std::string m_model_type;
    pp.query("model_type", m_model_type); // ns: default is "isothermal".
                                          // "energyD3Q27" activates product
                                          // equilibrium, energy equation etc.

    if (m_model_type == "energyD3Q27") {
        m_op.m_model_type = 1;
        amrex::Vector<amrex::Real> MachComponents{AMREX_D_DECL(0, 0, 0)};
        pp.queryarr("MachComponents", MachComponents, 0, AMREX_SPACEDIM);
        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.MachComponents[n] = MachComponents[n];
        }

        pp.query(
            "initialTemperature",
            m_op.initialTemperature); // ns:initial condition temperature
        pp.query(
            "adiabaticExponent",
            m_op.adiabaticExponent); // ns: reference gamma. safety block. do
                                     // not enable. not implemented.
        pp.query("meanMolecularMass", m_op.m_bar); // ns: reference m_bar
        m_op.speedOfSound_Ref = std::sqrt(
            m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
            m_op.initialTemperature); // set the actual speed of sound

        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
        }
    }
}

thermalDiffusivityTest::thermalDiffusivityTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("waveLength", m_op.waveLength); // ns:wavelength of the
                                             // perturbation
    std::string m_model_type;
    pp.query("model_type", m_model_type); // ns: default is "isothermal".
                                          // "energyD3Q27" activates product
                                          // equilibrium, energy equation etc.

    if (m_model_type == "energyD3Q27") {
        m_op.m_model_type = 1;
        amrex::Vector<amrex::Real> MachComponents{AMREX_D_DECL(0, 0, 0)};
        pp.queryarr("MachComponents", MachComponents, 0, AMREX_SPACEDIM);
        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.MachComponents[n] = MachComponents[n];
        }

        pp.query(
            "initialTemperature",
            m_op.initialTemperature); // ns:initial condition temperature
        pp.query(
            "adiabaticExponent",
            m_op.adiabaticExponent);               // ns: reference gamma.
        pp.query("meanMolecularMass", m_op.m_bar); // ns: reference m_bar
        m_op.speedOfSound_Ref = std::sqrt(
            m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
            m_op.initialTemperature); // set the actual speed of sound

        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
        }
    }
}

sodTest::sodTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("densityRatio", m_op.densityRatio);         // ns: R to L
    pp.query("temperatureRatio", m_op.temperatureRatio); // ns: R to L
    pp.query("xDiscontinuity", m_op.xDiscontinuity);     // ns: R to L

    std::string m_model_type;
    pp.query("model_type", m_model_type); // ns: default is "isothermal".
                                          // "energyD3Q27" activates product
                                          // equilibrium, energy equation etc.

    if (m_model_type == "energyD3Q27") {
        m_op.m_model_type = 1;
        amrex::Vector<amrex::Real> MachComponents{AMREX_D_DECL(0, 0, 0)};
        pp.queryarr("MachComponents", MachComponents, 0, AMREX_SPACEDIM);
        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.MachComponents[n] = MachComponents[n];
        }

        pp.query(
            "initialTemperature",
            m_op.initialTemperature); // ns:initial condition temperature
        pp.query(
            "adiabaticExponent",
            m_op.adiabaticExponent);               // ns: reference gamma.
        pp.query("meanMolecularMass", m_op.m_bar); // ns: reference m_bar
        m_op.speedOfSound_Ref = std::sqrt(
            m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
            m_op.initialTemperature); // set the actual speed of sound

        for (int n = 0; n < MachComponents.size(); n++) {
            m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
        }
    }
}

} // namespace lbm::ic

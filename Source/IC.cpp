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
    pp.query("model_type", m_model_type);

    m_op.m_model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.MachComponents[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initialTemperature);
    pp.query("adiabatic_exponent", m_op.adiabaticExponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speedOfSound_Ref = std::sqrt(
        m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
        m_op.initialTemperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
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

ViscosityTest::ViscosityTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("waveLength", m_op.waveLength);
    std::string m_model_type;
    pp.query("model_type", m_model_type);

    m_op.m_model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.MachComponents[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initialTemperature);
    pp.query("adiabatic_exponent", m_op.adiabaticExponent);

    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speedOfSound_Ref = std::sqrt(
        m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
        m_op.initialTemperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
    }
}

ThermalDiffusivityTest::ThermalDiffusivityTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("waveLength", m_op.waveLength);
    std::string m_model_type;
    pp.query("model_type", m_model_type);

    m_op.m_model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.MachComponents[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initialTemperature);
    pp.query("adiabatic_exponent", m_op.adiabaticExponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speedOfSound_Ref = std::sqrt(
        m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
        m_op.initialTemperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
    }
}

SodTest::SodTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("density_ratio", m_op.densityRatio);
    pp.query("temperature_ratio", m_op.temperatureRatio);
    pp.query("x_discontinuity", m_op.xDiscontinuity);

    std::string m_model_type;
    pp.query("model_type", m_model_type);

    m_op.m_model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.MachComponents[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initialTemperature);
    pp.query("adiabatic_exponent", m_op.adiabaticExponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speedOfSound_Ref = std::sqrt(
        m_op.adiabaticExponent * (m_op.R_u / m_op.m_bar) *
        m_op.initialTemperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.MachComponents[n] * m_op.speedOfSound_Ref;
    }
}

} // namespace lbm::ic

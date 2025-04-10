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

    std::string model_type;
    pp.query("model_type", model_type);

    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
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

    pp.query("wave_length", m_op.wave_length);
    std::string model_type;
    pp.query("model_type", model_type);

    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);

    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
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

    pp.query("wave_length", m_op.wave_length);
    std::string model_type;
    pp.query("model_type", model_type);

    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
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

    pp.query("density_ratio", m_op.density_ratio);
    pp.query("temperature_ratio", m_op.temperature_ratio);
    pp.query("x_discontinuity", m_op.x_discontinuity);

    std::string model_type;
    pp.query("model_type", model_type);

    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
    }
}

} // namespace lbm::ic

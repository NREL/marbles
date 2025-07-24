#include "VelocityBC.H"

namespace lbm::bc {
NoOp::NoOp() = default;

Constant::Constant()
{
    amrex::ParmParse pp("velocity_bc_constant");
    pp.query("dir", m_op.dir);
    pp.query("u0", m_op.u0);

    m_op.model_type = 1;
    pp.query("Mach_ref", m_op.Mach_ref);
    pp.query("initial_density", m_op.initial_density);
    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);
    m_op.u0 = m_op.Mach_ref * m_op.speed_of_sound_ref;
}

Channel::Channel()
{
    amrex::ParmParse pp("velocity_bc_channel");
    pp.query("u_ref", m_op.u_ref);

    m_op.model_type = 1;
    pp.query("Mach_ref", m_op.Mach_ref);
    pp.query("initial_density", m_op.initial_density);
    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);
    m_op.u_ref = m_op.Mach_ref * m_op.speed_of_sound_ref;
}

Parabolic::Parabolic()
{
    amrex::ParmParse pp("velocity_bc_parabolic");
    pp.query("normal_dir", m_op.normal_dir);
    pp.query("tangential_dir", m_op.tangential_dir);
    pp.query("um", m_op.um);

    m_op.model_type = 1;
    pp.query("Mach_ref", m_op.Mach_ref);
    pp.query("initial_density", m_op.initial_density);
    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);
    m_op.um = m_op.Mach_ref * m_op.speed_of_sound_ref;
}

} // namespace lbm::bc

#include "VelocityBC.H"

namespace lbm::bc {
NoOp::NoOp() {}

Constant::Constant()
{
    amrex::ParmParse pp("constant");
    pp.query("dir", m_op.dir);
    pp.query("u0", m_op.u0);
}

Channel::Channel()
{
    amrex::ParmParse pp("channel");
    pp.query("u_ref", m_op.u_ref);
}

} // namespace lbm::bc

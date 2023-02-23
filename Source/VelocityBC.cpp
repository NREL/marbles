#include "VelocityBC.H"

namespace lbm::bc {
Channel::Channel()
{
    amrex::ParmParse pp("channel");
    pp.query("u_ref", m_op.u_ref);
}

} // namespace lbm::bc

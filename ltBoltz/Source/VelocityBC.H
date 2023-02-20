#ifndef VELOCITYBC_H
#define VELOCITYBC_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>

namespace lbm::bc {

struct Channel
{
    struct DeviceOp
    {
        amrex::Real u_ref{1.0};

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::RealVect& vel_bc) const
        {
            const int* domhi = geom.Domain().hiVect();
            const amrex::Real c1 = iv[1] * (domhi[1] - iv[1]);
            const amrex::Real c2 = iv[2] * (domhi[2] - iv[2]);
            vel_bc[0] = 16.0 * u_ref * c1 * c2 / std::pow((domhi[2] + 1), 4);
        }
    };

    using DeviceType = DeviceOp;

    explicit Channel();

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace lbm::bc
#endif

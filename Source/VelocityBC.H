#ifndef VELOCITYBC_H
#define VELOCITYBC_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>

namespace lbm::bc {

struct NoOp
{
    struct DeviceOp
    {
        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& /*iv*/,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geom*/,
            amrex::RealVect& /*vel*/) const
        {}
    };

    using DeviceType = DeviceOp;

    explicit NoOp();

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct Constant
{
    struct DeviceOp
    {
        int dir{1};
        amrex::Real u0{0.0};

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& /*iv*/,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geom*/,
            amrex::RealVect& vel) const
        {
            vel[dir] = u0;
        }
    };

    using DeviceType = DeviceOp;

    explicit Constant();

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct Channel
{
    struct DeviceOp
    {
        amrex::Real u_ref{1.0};

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& iv,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& geom,
            amrex::RealVect& vel) const
        {
            const int* domhi = geom.Domain().hiVect();
            const amrex::Real c1 = iv[1] * (domhi[1] - iv[1]);
            const amrex::Real c2 = iv[2] * (domhi[2] - iv[2]);
            vel[0] = 16.0 * u_ref * c1 * c2 / std::pow((domhi[2] + 1), 4);
        }
    };

    using DeviceType = DeviceOp;

    explicit Channel();

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct Parabolic
{
    struct DeviceOp
    {
        int normal_dir{0};
        int tangential_dir{1};
        amrex::Real um{1.0};

        AMREX_GPU_DEVICE
        inline void operator()(
            const amrex::IntVect& iv,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& geom,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* prob_hi = geom.ProbHi();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real height =
                prob_hi[normal_dir] - prob_lo[normal_dir];
            const amrex::Real x =
                prob_lo[normal_dir] + (iv[normal_dir] + 0.5) * dx[normal_dir];
            vel[tangential_dir] =
                4.0 * um * x * (height - x) / (height * height);
        }
    };

    using DeviceType = DeviceOp;

    explicit Parabolic();

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace lbm::bc
#endif

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

        // ns: The initialisation below is no longer necessary. The product form
        // of initialization is always used. For isotermanl, RT = 1/3 as
        // initialised

        // AMREX_GPU_DEVICE
        // void operator()(
        //     const amrex::IntVect& /*iv*/,
        //     const int /*idir*/,
        //     const int /*ndir*/,
        //     const amrex::Real /*time*/,
        //     amrex::GeometryData const& /*geom*/,
        //     amrex::RealVect& /*vel*/) const
        // {}

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& /*iv*/,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geom*/,
            amrex::RealVect& /*vel*/,
            amrex::Real& /*R*/,
            amrex::Real& /*T*/,
            amrex::Real& /*gamma*/) const
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

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit
        amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, always true for monoatomic gas, i.e.
                       // single lattice. Non-dimensional. Safety constant
                       // declaration.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        std::string m_model_type = "isothermal";
        amrex::Real Mach_0{0.05}; // ns: mean Mach number

        // ns: The initialisation below is no longer necessary. The product form
        // of initialization is always used. For isotermanl, RT = 1/3 as
        // initialised

        // AMREX_GPU_DEVICE
        // void operator()(
        //     const amrex::IntVect& /*iv*/,
        //     const int /*idir*/,
        //     const int /*ndir*/,
        //     const amrex::Real /*time*/,
        //     amrex::GeometryData const& /*geom*/,
        //     amrex::RealVect& vel) const
        // {
        //     vel[dir] = u0;
        // }

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& /*iv*/,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geom*/,
            amrex::RealVect& vel,
            amrex::Real& R,
            amrex::Real& T,
            amrex::Real& gamma) const
        {
            vel[dir] = u0;
            // RT=1.0/3.0;
            R = (R_u / m_bar);
            T = initialTemperature;
            gamma = adiabaticExponent;
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

        // ns: The initialisation below is no longer necessary. The product form
        // of initialization is always used. For isotermanl, RT = 1/3 as
        // initialised

        // AMREX_GPU_DEVICE
        // void operator()(
        //     const amrex::IntVect& iv,
        //     const int /*idir*/,
        //     const int /*ndir*/,
        //     const amrex::Real /*time*/,
        //     amrex::GeometryData const& geom,
        //     amrex::RealVect& vel) const
        // {
        //     const int* domhi = geom.Domain().hiVect();
        //     const amrex::Real c1 = iv[1] * (domhi[1] - iv[1]);
        //     const amrex::Real c2 =
        //         AMREX_D_PICK(, c1, iv[2] * (domhi[2] - iv[2]));
        //     vel[0] = 16.0 * u_ref * c1 * c2 / std::pow((domhi[1] + 1), 4);
        // }

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& geom,
            amrex::RealVect& vel,
            amrex::Real& R,
            amrex::Real& T,
            amrex::Real& gamma) const
        {
            const int* domhi = geom.Domain().hiVect();
            const amrex::Real c1 = iv[1] * (domhi[1] - iv[1]);
            const amrex::Real c2 =
                AMREX_D_PICK(, c1, iv[2] * (domhi[2] - iv[2]));
            vel[0] = 16.0 * u_ref * c1 * c2 / std::pow((domhi[1] + 1), 4);
            R = 1.0;           //(R_u/m_bar)*initialTemperature;
            T = 1.0 / 3.0;     //(R_u/m_bar)*initialTemperature;
            gamma = 5.0 / 3.0; // ns: monoatomic
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

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit
        amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, always true for monoatomic gas, i.e.
                       // single lattice. Non-dimensional. Safety constant
                       // declaration.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        std::string m_model_type = "isothermal";
        amrex::Real Mach_m{0.05}; // ns: mean Mach number

        // ns: The initialisation below is no longer necessary. The product form
        // of initialization is always used. For isotermanl, RT = 1/3 as
        // initialised

        // AMREX_GPU_DEVICE
        // void operator()(
        //     const amrex::IntVect& iv,
        //     const int /*idir*/,
        //     const int /*ndir*/,
        //     const amrex::Real /*time*/,
        //     amrex::GeometryData const& geom,
        //     amrex::RealVect& vel) const
        // {
        //     const amrex::Real* prob_lo = geom.ProbLo();
        //     const amrex::Real* prob_hi = geom.ProbHi();
        //     const amrex::Real* dx = geom.CellSize();
        //     const amrex::Real height =
        //         prob_hi[normal_dir] - prob_lo[normal_dir];
        //     const amrex::Real x =
        //         prob_lo[normal_dir] + (iv[normal_dir] + 0.5) *
        //         dx[normal_dir];

        //         vel[tangential_dir] =
        //         4.0 * um * x * (height - x) / (height * height);

        // }

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& geom,
            amrex::RealVect& vel,
            amrex::Real& R,
            amrex::Real& T,
            amrex::Real& gamma) const
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

            R = (R_u / m_bar);
            T = initialTemperature;
            gamma = adiabaticExponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit Parabolic();

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace lbm::bc
#endif

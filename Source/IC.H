#ifndef IC_H
#define IC_H
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include "Stencil_D3Q27.H"
#include "Utilities.H"

namespace lbm::ic {

struct Constant
{
    struct DeviceOp
    {
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect MachComponents = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit
        amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, monoatomic gas. 1.4 for air,
                       // non-dimensional.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        // ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium,
        // energy equation etc.
        std::string m_model_type = "isothermal";

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& /*iv*/,
            amrex::GeometryData const& /*geom*/,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            rho = density;
            vel = velocity;
        }

        AMREX_GPU_DEVICE
        void thermalInitialise(
            const amrex::IntVect& /*iv*/,
            amrex::GeometryData const& /*geom*/,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {
            rho = density;
            vel = velocity;

            temperature = initialTemperature;
            R = R_u / m_bar; // ns: Specific gas constant. Implemented.
            gamma = adiabaticExponent; // ns: not implemented yet
        }
    };

    using DeviceType = DeviceOp;

    explicit Constant();

    static std::string identifier() { return "ic_constant"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct TaylorGreen
{
    struct DeviceOp
    {
        amrex::Real rho0{1.0};
        amrex::Real v0{1.0};
        amrex::RealVect omega = {AMREX_D_DECL(1.0, 1.0, 1.0)};
        amrex::Real L{1.0 / M_PI};

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit
        const amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, monoatomic gas. 1.4 for air,
                       // non-dimensional. Safety constant declaration.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        // ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium,
        // energy equation etc.
        std::string m_model_type = "isothermal";

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5) * dx[1];
            const amrex::Real z =
                AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5) * dx[2]);
            const amrex::Real omega0 = omega[0];
            const amrex::Real omega1 = omega[1];
            const amrex::Real omega2 = AMREX_D_PICK(, 0.0, omega[2]);

            rho = rho0 + rho0 * v0 * v0 / 16.0 *
                             (std::cos(2.0 * omega0 * x / L) +
                              std::cos(2.0 * omega1 * y / L)) *
                             (std::cos(2.0 * omega2 * z / L) + 2.0);
            AMREX_D_DECL(
                vel[0] = v0 * std::sin(omega0 * x / L) *
                         std::cos(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[1] = -v0 * std::cos(omega0 * x / L) *
                         std::sin(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[2] = 0.0);
        }

        AMREX_GPU_DEVICE
        void thermalInitialise(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5) * dx[1];
            const amrex::Real z =
                AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5) * dx[2]);
            const amrex::Real omega0 = omega[0];
            const amrex::Real omega1 = omega[1];
            const amrex::Real omega2 = AMREX_D_PICK(, 0.0, omega[2]);

            rho = rho0 + rho0 * v0 * v0 / 16.0 *
                             (std::cos(2.0 * omega0 * x / L) +
                              std::cos(2.0 * omega1 * y / L)) *
                             (std::cos(2.0 * omega2 * z / L) + 2.0);
            AMREX_D_DECL(
                vel[0] = v0 * std::sin(omega0 * x / L) *
                         std::cos(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[1] = -v0 * std::cos(omega0 * x / L) *
                         std::sin(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[2] = 0.0);

            temperature = initialTemperature;
            R = R_u / m_bar; // ns: Specific gas constant. Implemented.
            gamma = adiabaticExponent; // ns: not implemented yet
        }
    };

    using DeviceType = DeviceOp;

    explicit TaylorGreen();

    static std::string identifier() { return "ic_taylorgreen"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct viscosityTest
{
    struct DeviceOp
    {
        amrex::Real waveLength{1.0}; // number of points in Y
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect MachComponents = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit
        amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, monoatomic gas. 1.4 for air,
                       // non-dimensional. Safety constant declaration.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        // ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium,
        // energy equation etc.
        std::string m_model_type = "isothermal";

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            // const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];
            // const amrex::Real z =
            //     AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5 * 0.0) *
            //     dx[2]);

            rho = density;

            const amrex::Real a0 =
                0.010 * speedOfSound_Ref; // ns:the perturbation amplitude

            vel[0] = velocity[0] + a0 * std::sin(2.0 * M_PI * y / waveLength);
            vel[1] = velocity[1];
            if (AMREX_SPACEDIM == 3) vel[2] = velocity[2];
        }

        AMREX_GPU_DEVICE
        void thermalInitialise(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {

            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            // const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];
            // const amrex::Real z =
            //     AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5 * 0.0) *
            //     dx[2]);

            rho = density;

            const amrex::Real a0 =
                0.010 * speedOfSound_Ref; // ns:the perturbation amplitude

            vel[0] = velocity[0] + a0 * std::sin(2.0 * M_PI * y / waveLength);
            vel[1] = velocity[1];
            if (AMREX_SPACEDIM == 3) vel[2] = velocity[2];

            temperature = initialTemperature;
            R = R_u / m_bar; // ns: Specific gas constant. Implemented.
            gamma = adiabaticExponent; // ns: not implemented yet
        }
    };

    using DeviceType = DeviceOp;

    explicit viscosityTest();

    static std::string identifier() { return "ic_viscosityTest"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct thermalDiffusivityTest
{
    struct DeviceOp
    {
        amrex::Real waveLength{1.0}; // number of points in Y
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect MachComponents = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit
        amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, monoatomic gas. 1.4 for air,
                       // non-dimensional. Safety constant declaration.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        // ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium,
        // energy equation etc.
        std::string m_model_type = "isothermal";

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            // const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];
            // const amrex::Real z =
            //     AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5 * 0.0) *
            //     dx[2]);

            const amrex::Real a0 =
                0.0010 * initialTemperature; // ns:the perturbation amplitude

            rho = density + a0 * std::sin(2.0 * M_PI * y / waveLength);

            vel[0] = velocity[0]; //+a0*std::sin(2.0*M_PI*y/waveLength);
            vel[1] = velocity[1];
            if (AMREX_SPACEDIM == 3) vel[2] = velocity[2];
        }

        AMREX_GPU_DEVICE
        void thermalInitialise(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {

            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            // const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];
            // const amrex::Real z =
            //     AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5 * 0.0) *
            //     dx[2]);

            const amrex::Real a0 =
                0.0010 * initialTemperature; // ns:the perturbation amplitude

            R = R_u / m_bar; // ns: Specific gas constant. Implemented.
            const amrex::Real pressure = density * R * initialTemperature;
            rho = density + a0 * std::sin(2.0 * M_PI * y / waveLength);

            vel[0] = velocity[0];
            vel[1] = velocity[1];
            if (AMREX_SPACEDIM == 3) vel[2] = velocity[2];

            gamma = adiabaticExponent; // ns: not implemented yet

            temperature = pressure / (rho * R);
        }
    };

    using DeviceType = DeviceOp;

    explicit thermalDiffusivityTest();

    static std::string identifier() { return "ic_thermalDiffusivityTest"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct sodTest
{
    struct DeviceOp
    {
        amrex::Real densityRatio{1.0}; // ns:
        amrex::Real temperatureRatio{1.0};
        amrex::Real xDiscontinuity{10.0}; // ns: discontinuity location
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect MachComponents = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initialTemperature =
            1.0 / 3.0; // ns: Safe default, lattice temperature, LB unit

        amrex::Real adiabaticExponent =
            5.0 / 3.0; // ns: Safe default, monoatomic gas. 1.4 for air,
                       // non-dimensional. Safety constant declaration.
        const amrex::Real R_u =
            28.96; // ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar =
            28.96; // ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref = std::sqrt(
            adiabaticExponent * (R_u / m_bar) *
            initialTemperature); // the actual speed of sound

        // ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium,
        // energy equation etc.
        std::string m_model_type = "isothermal";

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];
            // const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];
            // const amrex::Real z =
            //     AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5 * 0.0) *
            //     dx[2]);

            vel[0] = velocity[0];
            vel[1] = velocity[1];
            if (AMREX_SPACEDIM == 3) vel[2] = velocity[2];

            rho = density + 0.5 *
                                (1.0 + std::tanh((x - xDiscontinuity) * 3.0)) *
                                (densityRatio * density - density);
        }

        AMREX_GPU_DEVICE
        void thermalInitialise(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {

            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];
            // const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];
            // const amrex::Real z =
            //     AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5 * 0.0) *
            //     dx[2]);

            R = R_u / m_bar; // ns: Specific gas constant. Implemented.

            vel[0] = velocity[0];
            vel[1] = velocity[1];
            if (AMREX_SPACEDIM == 3) vel[2] = velocity[2];

            gamma = adiabaticExponent; // ns: not implemented yet

            rho = density + 0.5 *
                                (1.0 + std::tanh((x - xDiscontinuity) * 3.0)) *
                                (densityRatio * density - density);
            temperature = initialTemperature +
                          0.5 * (1.0 + std::tanh((x - xDiscontinuity) * 3.0)) *
                              (temperatureRatio * initialTemperature -
                               initialTemperature);
        }
    };

    using DeviceType = DeviceOp;

    explicit sodTest();

    static std::string identifier() { return "ic_sodTest"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

class InitializerBase
{
public:
    InitializerBase() = default;

    virtual ~InitializerBase() = default;

    virtual void initialize(const int lev, amrex::GeometryData const& geom) = 0;
};

template <typename ICOp>
class Initializer : public InitializerBase
{
public:
    using ICDeviceOp = typename ICOp::DeviceType;

    Initializer(
        const amrex::Real& mesh_speed,
        const ICOp ic_op,
        amrex::Vector<amrex::MultiFab>& f)
        : m_mesh_speed(mesh_speed), m_op(ic_op), m_f(f)
    {}

    Initializer(
        const amrex::Real& mesh_speed,
        const ICOp ic_op,
        amrex::Vector<amrex::MultiFab>& f,
        amrex::Vector<amrex::MultiFab>& g)
        : m_mesh_speed(mesh_speed), m_op(ic_op), m_f(f), m_g(g)
    {}

    ~Initializer() override = default;

    void initialize(const int lev, amrex::GeometryData const& geom) override
    {
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        const auto& weight = stencil.weights;
        auto const& f_arrs = m_f[lev].arrays();
        auto const& g_arrs = m_g[lev].arrays();
        const amrex::Real l_mesh_speed = m_mesh_speed;
        ICDeviceOp ic(m_op.device_instance());

        amrex::RealVect zeroVec = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::ParallelFor(
            m_f[lev], m_f[lev].nGrowVect(), constants::N_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                const amrex::Real wt = weight[q];
                const auto& ev = evs[q];

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real rho = 1.0;
                amrex::RealVect vel(AMREX_D_DECL(0.0, 0.0, 0.0));
                amrex::Real temperature = 1.0 / 3.0;
                amrex::Real R = ic.R_u / 28.96;
                amrex::Real gamma = 1.667;

                if (ic.m_model_type == "energyD3Q27")
                    ic.thermalInitialise(
                        iv, geom, rho, vel, temperature, R, gamma);
                else
                    ic(iv, geom, rho, vel);

                if (ic.m_model_type == "energyD3Q27") {
                    set_equilibrium_value_D3Q27(
                        rho, vel, R * temperature, l_mesh_speed, wt, ev,
                        f_arrs[nbx](i, j, k, q));

                    amrex::Real Cv = R / (gamma - 1.0);

                    amrex::Real twoRhoE = getEnergy(temperature, rho, vel, Cv);

                    amrex::Real qxEq, qyEq, qzEq, RxxEq, RyyEq, RzzEq, RxyEq,
                        RxzEq, RyzEq;

                    getEquilibriumMoments(
                        rho, vel[0], vel[1], vel[2], twoRhoE, Cv, R, qxEq, qyEq,
                        qzEq, RxxEq, RyyEq, RzzEq, RxyEq, RxzEq, RyzEq);

                    amrex::RealVect heatFlux = {AMREX_D_DECL(qxEq, qyEq, qzEq)};

                    amrex::Vector<amrex::Real> fluxOfHeatFlux = {
                        RxxEq, RyyEq, RzzEq, RxyEq, RxzEq, RyzEq};

                    set_extended_gradExpansion_generic(
                        twoRhoE, heatFlux, fluxOfHeatFlux, l_mesh_speed, wt, ev,
                        stencil.theta0, zeroVec, 1.0, g_arrs[nbx](i, j, k, q));

                } else
                    set_equilibrium_value(
                        rho, vel, l_mesh_speed, wt, ev,
                        f_arrs[nbx](i, j, k, q));
            });
        amrex::Gpu::synchronize();
    }

    const amrex::Real m_mesh_speed;
    const ICOp m_op;
    amrex::Vector<amrex::MultiFab>& m_f;
    amrex::Vector<amrex::MultiFab>& m_g;
};

} // namespace lbm::ic
#endif

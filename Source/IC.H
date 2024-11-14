#ifndef IC_H
#define IC_H
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include "Stencil.H"
#include "Utilities.H"

namespace lbm::ic {

struct Constant
{
    struct DeviceOp
    {
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect MachComponents = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        //amrex::Vector<amrex::Real> MachComponents = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initialTemperature=1.0/3.0; //ns: Safe default, lattice temperature, LB unit
        const amrex::Real adiabaticExponent=5.0/3.0;  //ns: Safe default, monoatomic gas. 1.4 for air, non-dimensional. Safety constant decalaration.
        const amrex::Real R_u=28.96;              //ns: Universal Gas constant, LB unit. Never change.
        amrex::Real m_bar=28.96;                //ns: mean molecular mass, Safe default air. Unit kg/kmol.
        amrex::Real speedOfSound_Ref=std::sqrt(adiabaticExponent*(R_u/m_bar)*initialTemperature); //the actual speed of sound

        //ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium, energy equation etc.
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

        //ns: overloading did not work
        AMREX_GPU_DEVICE
        void thermalInitialise
        (
            const amrex::IntVect& /*iv*/,
            amrex::GeometryData const& /*geom*/,
            amrex::Real& rho,
            amrex::RealVect& vel, amrex::Real& temperature, amrex::Real& R, amrex::Real& gamma) const
        {
            rho = density;
            vel = velocity;
            
            //for (int n = 0; n < MachComponents.size(); n++) velocity[n] = MachComponents[n]*speedOfSound_Ref;
            
            //for (int n = 0; n < MachComponents.size(); n++) vel[n] = MachComponents[n]*speedOfSound_Ref;
            
            temperature = initialTemperature;
            R = R_u/m_bar;             //ns: Specific gas constant. Implemented.
            gamma = adiabaticExponent; //ns: not implemented yet            
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
        //ns: LBM Scheme selector. "energyD3Q27" activates product equilibrium, energy equation etc.
        const amrex::Real R_u=28.96;              //ns: Universal Gas constant, LB unit. Never change.
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
            amrex::RealVect& vel, amrex::Real& temperature, amrex::Real& R, amrex::Real& gamma) const
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

    };

    using DeviceType = DeviceOp;

    explicit TaylorGreen();

    static std::string identifier() { return "ic_taylorgreen"; }

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

    ~Initializer() override = default;

    void initialize(const int lev, amrex::GeometryData const& geom) override
    {
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        const auto& weight = stencil.weights;
        auto const& f_arrs = m_f[lev].arrays();
        const amrex::Real l_mesh_speed = m_mesh_speed;
        ICDeviceOp ic(m_op.device_instance());
        amrex::ParallelFor(
            m_f[lev], m_f[lev].nGrowVect(), constants::N_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                const amrex::Real wt = weight[q];
                const auto& ev = evs[q];

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real rho = 1.0;
                amrex::RealVect vel(AMREX_D_DECL(0.0, 0.0, 0.0));
                amrex::Real temperature = 1.0/3.0;
                amrex::Real R = ic.R_u/28.96;
                amrex::Real gamma = 1.667;
                
                if (ic.m_model_type=="energyD3Q27") ic.thermalInitialise(iv, geom, rho, vel, temperature, R, gamma);
                else ic(iv, geom, rho, vel);

                if (ic.m_model_type=="energyD3Q27") set_equilibrium_value_D3Q27(rho, vel, R*temperature, l_mesh_speed, wt, ev, f_arrs[nbx](i, j, k, q));
                else set_equilibrium_value(rho, vel, l_mesh_speed, wt, ev, f_arrs[nbx](i, j, k, q));

            });
        amrex::Gpu::synchronize();
    }

    const amrex::Real m_mesh_speed;
    const ICOp m_op;
    amrex::Vector<amrex::MultiFab>& m_f;
};

} // namespace lbm::ic
#endif

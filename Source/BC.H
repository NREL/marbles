#ifndef BC_H
#define BC_H
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include "Constants.H"
#include "Utilities.H"
#include "Stencil.H"
#include "VelocityBC.H"

namespace lbm::bc {

static constexpr int periodic = 0;
static constexpr int noslipwall = 1;
static constexpr int velocity = 2;
static constexpr int pressure = 3;
static constexpr int outflow = 4;

template <typename VelBCOp>
struct BCFill
{

    using VelBCOpType = VelBCOp;

    const amrex::GpuArray<amrex::IntVect, constants::n_micro_states> m_evs =
        stencil::Stencil().evs;
    const amrex::GpuArray<int, constants::n_micro_states> m_bounce_dirs =
        stencil::Stencil().bounce_dirs;
    const amrex::GpuArray<amrex::Real, constants::n_micro_states> m_weights =
        stencil::Stencil().weights;

    const amrex::Real m_mesh_speed;
    const amrex::GpuArray<int, AMREX_SPACEDIM* 2> m_bc_type = {0};
    const amrex::IntVect m_ng = {0, 0, 0};
    const VelBCOpType m_vel_bc_op;

    AMREX_GPU_HOST
    constexpr explicit BCFill(
        const amrex::Real& mesh_speed,

        const amrex::GpuArray<int, AMREX_SPACEDIM * 2>& bc_type,
        const amrex::IntVect& ng,
        const VelBCOpType& vel_bc_op)
        : m_mesh_speed(mesh_speed)
        , m_bc_type(bc_type)
        , m_ng(ng)
        , m_vel_bc_op(vel_bc_op)
    {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bounce_back_bc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to bounce back. 2D e.g., for
        // the ghost cell at at (-1, -1), store the
        // value f_(0, 0) associated with e=[-1, -1] and
        // orient it such that it points e=[1, 1]. This
        // reproduces the bounce back effect directly
        // when the streaming step is performed.

        data(iv, q) = data(iv + ev, m_bounce_dirs[q]);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void velocity_bc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int idim,
        const int lohi,
        const amrex::GpuArray<int, 2>& wall_normal_unit,
        const amrex::RealVect& vel_bc,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to a prescribed velocity. This
        // is done by assigning all the values of f_ on
        // this ghost layer to be the equilibrium
        // distributions associated with the desired
        // velocity. They will naturally exert their
        // effect on the interior fluid points during
        // the streaming step.

        // In this case, a value for u, v, and w must be
        // specified, either hard coded or set by a
        // user-specified function
        amrex::Real rho_out = 0.0;
        amrex::Real rho_tan = 0.0;

        const amrex::IntVect ivn(iv + ev);

        // Build the 3-part rho calculation at the
        // (iv+ev) fluid cell
        for (int qq = 0; qq < constants::n_micro_states; qq++) {

            const auto& evn_out = m_evs[qq];
            const int bounce_qq = m_bounce_dirs[qq];
            const auto& evn_in = m_evs[bounce_qq];

            if (evn_in[idim] == -wall_normal_unit[lohi]) {
                rho_out += 2.0 * data(ivn + evn_out, bounce_qq);
            } else if (evn_in[idim] == 0) {
                rho_tan += 1.0 * data(ivn + evn_out, bounce_qq);
            }
        }

        const amrex::Real rho_bc =
            (rho_out + rho_tan) / (1.0 - wall_normal_unit[lohi] * vel_bc[idim]);

        SetEquilibriumValue(
            rho_bc, vel_bc, m_mesh_speed, m_weights[q], ev, data(iv, q));
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void pressure_bc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int idim,
        const int lohi,
        const amrex::GpuArray<int, 2>& wall_normal_unit,
        const amrex::Real rho_bc,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to a prescribed
        // pressure/density. This is done by assigning
        // all the values of f_ on this ghost layer to
        // be the equilibrium distributions associated
        // with the desired density. They will naturally
        // exert their effect on the interior fluid
        // points during the streaming step.

        amrex::Real rho_out = 0.0;
        amrex::Real rho_tan = 0.0;

        const amrex::IntVect ivn(iv + ev);

        // Build the 3-part rho calculation at the
        // (iv+ev) fluid cell
        for (int qq = 0; qq < constants::n_micro_states; qq++) {

            const auto& evn_out = m_evs[qq];
            const int bounce_qq = m_bounce_dirs[qq];
            const auto& evn_in = m_evs[bounce_qq];

            if (evn_in[idim] == -wall_normal_unit[lohi]) {
                rho_out += 2.0 * data(ivn + evn_out, bounce_qq);
            } else if (evn_in[idim] == 0) {
                rho_tan += 1.0 * data(ivn + evn_out, bounce_qq);
            }
        }

        amrex::RealVect vel_bc(0.0, 0.0, 0.0);
        vel_bc[idim] =
            wall_normal_unit[lohi] * (1.0 - (rho_out + rho_tan) / rho_bc);

        SetEquilibriumValue(
            rho_bc, vel_bc, m_mesh_speed, m_weights[q], ev, data(iv, q));
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void outflow_bc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int idim,
        const int lohi,
        const amrex::GpuArray<int, 2>& wall_normal_unit,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to an outflow condition This
        // copies the rate of change from the last two
        // interior cell onto the ghost layer such that
        // when streaming you have the effect of "no
        // change" in the direction across the boundary.
        amrex::IntVect iv_shift(amrex::IntVect::TheDimensionVector(idim));

        if (wall_normal_unit[lohi] * iv_shift[0] == ev[0] &&
            wall_normal_unit[lohi] * iv_shift[1] == ev[1] &&
            wall_normal_unit[lohi] * iv_shift[2] == ev[2]) {

            amrex::Real rho = 0.0;
            amrex::RealVect vel(0.0, 0.0, 0.0);
            for (int qq = 0; qq < constants::n_micro_states; qq++) {
                rho += data(iv + wall_normal_unit[lohi] * iv_shift, qq);
                const auto& ev_qq = m_evs[qq];
                vel[0] +=
                    ev_qq[0] * data(iv + wall_normal_unit[lohi] * iv_shift, qq);
                vel[1] +=
                    ev_qq[1] * data(iv + wall_normal_unit[lohi] * iv_shift, qq);
                vel[2] +=
                    ev_qq[2] * data(iv + wall_normal_unit[lohi] * iv_shift, qq);
            }

            vel *= m_mesh_speed / rho;

            amrex::Real f_eq = 0.0;
            amrex::Real f_star_eq = 0.0;

            const int bounce_q = m_bounce_dirs[q];
            const auto& ev_star = m_evs[bounce_q];

            SetEquilibriumValue(rho, vel, m_mesh_speed, m_weights[q], ev, f_eq);

            SetEquilibriumValue(
                rho, vel, m_mesh_speed, m_weights[bounce_q], ev_star,
                f_star_eq);

            data(iv, q) =
                f_eq + data(iv + wall_normal_unit[lohi] * iv_shift, bounce_q) -
                f_star_eq;

        } else {
            data(iv, q) =
                2.0 * data(iv + 1 * wall_normal_unit[lohi] * iv_shift, q) -
                data(iv + 2 * wall_normal_unit[lohi] * iv_shift, q);
        }
    }

    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& iv,
        amrex::Array4<amrex::Real> const& data,
        const int /*dcomp*/,
        const int /*numcomp*/,
        amrex::GeometryData const& geom,
        const amrex::Real /*time*/,
        const amrex::BCRec* /*bcr*/,
        const int /*bcomp*/,
        const int /*orig_comp*/) const
    {
        // do something for external Dirichlet (BCType::ext_dir)
        const int* domlo = geom.Domain().loVect();
        const int* domhi = geom.Domain().hiVect();

        // Define an inward-pointing normal vector
        // that is, positive on the x-, y-, z- wall and
        // negative on the x+, y+, z+ wall
        const amrex::GpuArray<int, 2> wall_normal_unit = {1, -1};

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            for (int lohi = 0; lohi < 2; lohi++) {

                const int bc_idx = idim + lohi * AMREX_SPACEDIM;

                if ((((lohi == 0) && (iv[idim] < domlo[idim])) ||
                     ((lohi == 1) && (iv[idim] > domhi[idim])))) {

                    const auto& lb = amrex::lbound(data);
                    const auto& ub = amrex::ubound(data);
                    const amrex::IntVect liv(
                        amrex::IntVect(lb.x, lb.y, lb.z) + m_ng);
                    const amrex::IntVect uiv(
                        amrex::IntVect(ub.x, ub.y, ub.z) - m_ng);
                    const amrex::Box inside(liv, uiv);

                    for (int q = 0; q < constants::n_micro_states; q++) {

                        const auto& ev = m_evs[q];
                        const amrex::IntVect ivn(iv + ev);

                        // If this distribution points into (enters) the
                        // domain, and we are accessing a cell within the data
                        // bounds
                        if ((ev[idim] == wall_normal_unit[lohi]) &&
                            (inside.contains(ivn))) {

                            if (m_bc_type[bc_idx] == 1) {
                                bounce_back_bc(iv, q, ev, data);
                            } else if (m_bc_type[bc_idx] == 2) {
                                amrex::RealVect vel_bc(0.0, 0.0, 0.0);
                                m_vel_bc_op(iv, geom, vel_bc);
                                velocity_bc(
                                    iv, q, ev, idim, lohi, wall_normal_unit,
                                    vel_bc, data);
                            } else if (m_bc_type[bc_idx] == 3) {
                                const amrex::Real rho_bc = 1.0;
                                pressure_bc(
                                    iv, q, ev, idim, lohi, wall_normal_unit,
                                    rho_bc, data);
                            } else if (m_bc_type[bc_idx] == 4) {
                                outflow_bc(
                                    iv, q, ev, idim, lohi, wall_normal_unit,
                                    data);
                            }
                        }
                    }
                }
            }
        }
    }
};

template <typename VelBCOp>
struct BCOpCreator
{
    using VelBCOpType = typename VelBCOp::DeviceType;
    using FunctorType = BCFill<VelBCOpType>;

    explicit BCOpCreator(
        const amrex::Real& mesh_speed,
        const amrex::GpuArray<int, AMREX_SPACEDIM * 2>& bc_type,
        const amrex::IntVect& ng)
        : m_mesh_speed(mesh_speed), m_bc_type(bc_type), m_ng(ng), m_vel_bc_op()
    {}

    inline FunctorType operator()() const
    {
        return FunctorType{
            m_mesh_speed, m_bc_type, m_ng, m_vel_bc_op.device_instance()};
    }

    const amrex::Real m_mesh_speed;
    const amrex::GpuArray<int, AMREX_SPACEDIM* 2> m_bc_type = {0};
    const amrex::IntVect m_ng = {0, 0, 0};
    VelBCOp m_vel_bc_op;
};

} // namespace lbm::bc
#endif

#ifndef UTILITIES_H
#define UTILITIES_H

#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "Constants.H"

namespace lbm {

AMREX_GPU_DEVICE AMREX_INLINE void set_equilibrium_value(
    const amrex::Real& rho,
    const amrex::RealVect& vel,
    const amrex::Real& mesh_speed,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    const amrex::Real e_dot_u =
        mesh_speed *
        (AMREX_D_TERM(ev[0] * vel[0], +ev[1] * vel[1], +ev[2] * vel[2]));
    const amrex::Real umag_2 =
        AMREX_D_TERM(vel[0] * vel[0], +vel[1] * vel[1], +vel[2] * vel[2]);
    const amrex::Real cs = mesh_speed / constants::ROOT3;
    const amrex::Real cs_2 = cs * cs;

    // If you make the assumption dx = dt = 1, then
    // c = dx/dt = 1/1 = 1 and
    // mesh_speed_of_sound = c/sqrt(3) = 1/sqrt(3)
    // and these simplifications fall out
    // const amrex::Real c1 = 3.0 * e_dot_u;
    // const amrex::Real c2 = 4.5 * e_dot_u * e_dot_u;
    // const amrex::Real c3 = -1.5 * umag_2;

    const amrex::Real c1 = e_dot_u / cs_2;
    const amrex::Real c2 = e_dot_u * e_dot_u / (2.0 * cs_2 * cs_2);
    const amrex::Real c3 = -umag_2 / (2.0 * cs_2);

    eq = rho * wt * (1.0 + c1 + c2 + c3);
}

AMREX_GPU_DEVICE AMREX_INLINE amrex::Real gradient(
    const int dir,
    const int comp,
    const amrex::IntVect& iv,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& idx,
    const amrex::Box& bx,
    const amrex::Array4<const int>& is_fluid,
    const amrex::Array4<const amrex::Real>& data)
{
    const amrex::IntVect ivp(iv + amrex::IntVect::TheDimensionVector(dir));
    const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
    amrex::Real vp = 0.0, vc = 0.0, vm = 0.0;

    if (((!bx.contains(ivp)) || (is_fluid(ivp, 0) != 1)) &&
        ((!bx.contains(ivm)) || (is_fluid(ivm, 0) != 1))) {
        vp = 0.0;
        vc = 0.0;
        vm = 0.0;
    } else if ((!bx.contains(ivp)) || (is_fluid(ivp, 0) != 1)) {
        vp = 0.0;
        vc = data(iv, comp);
        vm = -data(ivm, comp);
    } else if ((!bx.contains(ivm)) || (is_fluid(ivm, 0) != 1)) {
        vp = data(ivp, comp);
        vc = -data(iv, comp);
        vm = 0.0;
    } else {
        vp = 0.5 * data(ivp, comp);
        vc = 0.0;
        vm = -0.5 * data(ivm, comp);
    }

    return (vp + vc + vm) * idx[dir];
}

// version of amrex_avgdown that masks out certain values
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void masked_avgdown(
    int i,
    int j,
    int k,
    int n,
    amrex::Array4<T> const& crse,
    amrex::Array4<T const> const& fine,
    int ccomp,
    int fcomp,
    amrex::IntVect const& ratio,
    const T masked_val = -1.0) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const int ii = i * facx;
    const int jj = j * facy;
    const int kk = k * facz;
    T c = 0;
    T vol = 0;
    for (int kref = 0; kref < facz; ++kref) {
        for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                const T fval = fine(ii + iref, jj + jref, kk + kref, n + fcomp);
                if (std::abs(fval - masked_val) > constants::SMALL_NUM) {
                    c += fval;
                    vol += T(1.0);
                }
            }
        }
    }
    if (vol > 0.0) {
        crse(i, j, k, n + ccomp) = c / vol;
    }
}

// Version of average_down that also averages down fine ghost cells
void average_down_with_ghosts(
    const amrex::MultiFab& fine,
    amrex::MultiFab& crse,
    const amrex::Geometry& geom,
    const amrex::IntVect ng,
    const amrex::IntVect ref_ratio);

bool file_exists(const std::string& fname);
} // namespace lbm
#endif

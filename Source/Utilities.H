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
        mesh_speed * (ev[0] * vel[0] + ev[1] * vel[1] + ev[2] * vel[2]);
    const amrex::Real umag_2 =
        vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
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

    if (((!bx.contains(ivp)) || (is_fluid(ivp) != 1)) &&
        ((!bx.contains(ivm)) || (is_fluid(ivm) != 1))) {
        vp = 0.0;
        vc = 0.0;
        vm = 0.0;
    } else if ((!bx.contains(ivp)) || (is_fluid(ivp) != 1)) {
        vp = 0.0;
        vc = data(iv, comp);
        vm = -data(ivm, comp);
    } else if ((!bx.contains(ivm)) || (is_fluid(ivm) != 1)) {
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

// Version of average_down that also averages down fine ghost cells
void average_down_with_ghosts(
    const amrex::MultiFab& fine,
    amrex::MultiFab& crse,
    const amrex::IntVect ng,
    const amrex::IntVect ref_ratio);

bool file_exists(const std::string& fname);
} // namespace lbm
#endif

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

AMREX_GPU_DEVICE AMREX_INLINE void set_equilibrium_value_D3Q27(
    const amrex::Real& rho,
    const amrex::RealVect& vel,
    const amrex::Real& mesh_speed,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    //Warning! Achtung! Strictly for D3Q27 only
    //This entropic equilibrium is very stable, mostly positive and relaxes all 27 moments

    // const amrex::Real e_dot_u =
    //     mesh_speed *
    //     (AMREX_D_TERM(ev[0] * vel[0], +ev[1] * vel[1], +ev[2] * vel[2]));
    // const amrex::Real umag_2 =
    //     AMREX_D_TERM(vel[0] * vel[0], +vel[1] * vel[1], +vel[2] * vel[2]);
    // const amrex::Real cs = mesh_speed / constants::ROOT3;
    // const amrex::Real cs_2 = cs * cs;

    // If you make the assumption dx = dt = 1, then
    // c = dx/dt = 1/1 = 1 and
    // mesh_speed_of_sound = c/sqrt(3) = 1/sqrt(3)
    // and these simplifications fall out
    // const amrex::Real c1 = 3.0 * e_dot_u;
    // const amrex::Real c2 = 4.5 * e_dot_u * e_dot_u;
    // const amrex::Real c3 = -1.5 * umag_2;

    // const amrex::Real c1 = e_dot_u / cs_2;
    // const amrex::Real c2 = e_dot_u * e_dot_u / (2.0 * cs_2 * cs_2);
    // const amrex::Real c3 = -umag_2 / (2.0 * cs_2);

    // eq = rho * wt * (1.0 + c1 + c2 + c3);

    //write here
    const amrex::Real temperature = 1.0/3.0; //Warning! Achtung! Placeholder for isothermal and w/o correction only
    const amrex::Real Pxx = vel[0]*vel[0] + temperature;
    const amrex::Real Pyy = vel[1]*vel[1] + temperature;
    const amrex::Real Pzz = vel[2]*vel[2] + temperature;

    const amrex::Real Phix = ev[0]*0.5*vel[0] + std::abs(ev[0])*(1.50*Pxx-1.0)-Pxx+1.0;
    const amrex::Real Phiy = ev[1]*0.5*vel[1] + std::abs(ev[1])*(1.50*Pyy-1.0)-Pyy+1.0;
    const amrex::Real Phiz = ev[2]*0.5*vel[2] + std::abs(ev[2])*(1.50*Pzz-1.0)-Pzz+1.0;

    eq = rho * Phix * Phiy * Phiz;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_extended_equilibrium_value_D3Q27(
    const amrex::Real& rho,
    const amrex::RealVect& vel, const amrex::Real& Pxx, const amrex::Real& Pyy, const amrex::Real& Pzz,
    const amrex::Real& mesh_speed,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    //Warning! Achtung! Strictly for D3Q27 only
    //This entropic equilibrium is very stable, mostly positive and relaxes all 27 moments

    // const amrex::Real e_dot_u =
    //     mesh_speed *
    //     (AMREX_D_TERM(ev[0] * vel[0], +ev[1] * vel[1], +ev[2] * vel[2]));
    // const amrex::Real umag_2 =
    //     AMREX_D_TERM(vel[0] * vel[0], +vel[1] * vel[1], +vel[2] * vel[2]);
    // const amrex::Real cs = mesh_speed / constants::ROOT3;
    // const amrex::Real cs_2 = cs * cs;

    // If you make the assumption dx = dt = 1, then
    // c = dx/dt = 1/1 = 1 and
    // mesh_speed_of_sound = c/sqrt(3) = 1/sqrt(3)
    // and these simplifications fall out
    // const amrex::Real c1 = 3.0 * e_dot_u;
    // const amrex::Real c2 = 4.5 * e_dot_u * e_dot_u;
    // const amrex::Real c3 = -1.5 * umag_2;

    // const amrex::Real c1 = e_dot_u / cs_2;
    // const amrex::Real c2 = e_dot_u * e_dot_u / (2.0 * cs_2 * cs_2);
    // const amrex::Real c3 = -umag_2 / (2.0 * cs_2);

    // eq = rho * wt * (1.0 + c1 + c2 + c3);

    //write here
    //const amrex::Real temperature = 1.0/3.0; //Warning! Achtung! Placeholder for isothermal and w/o correction only
    //const amrex::Real Pxx = vel[0]*vel[0] + temperature;
    //const amrex::Real Pyy = vel[1]*vel[1] + temperature;
    //const amrex::Real Pzz = vel[2]*vel[2] + temperature;

    const amrex::Real Phix = ev[0]*0.5*vel[0] + std::abs(ev[0])*(1.50*Pxx-1.0)-Pxx+1.0;
    const amrex::Real Phiy = ev[1]*0.5*vel[1] + std::abs(ev[1])*(1.50*Pyy-1.0)-Pyy+1.0;
    const amrex::Real Phiz = ev[2]*0.5*vel[2] + std::abs(ev[2])*(1.50*Pzz-1.0)-Pzz+1.0;

    eq = rho * Phix * Phiy * Phiz;
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
    const int facz = AMREX_D_PICK(1, 1, ratio[2]);
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

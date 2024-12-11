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
    const amrex::Real& Rtemperature,
    const amrex::Real& mesh_speed,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    // Warning! Achtung! Strictly for D3Q27 only
    // This entropic equilibrium is very stable. Relaxes all 27 moments.
    // Unfortunately D3Q27 lattice has only 24 independent moments. For the
    // standard lattice, e.g. D3Q27, c_x^3 = c_x , think (-1)^3 = -1. So Qxxx,
    // Qyyy and Qzzz are always wrong. We correct them in Pxx, Pyy and Pzz
    // during collision.

    const amrex::Real Pxx = vel[0] * vel[0] + Rtemperature;
    const amrex::Real Pyy = vel[1] * vel[1] + Rtemperature;
    const amrex::Real Pzz = vel[2] * vel[2] + Rtemperature;

    const amrex::Real Phix =
        ev[0] * 0.5 * vel[0] + std::abs(ev[0]) * (1.50 * Pxx - 1.0) - Pxx + 1.0;
    const amrex::Real Phiy =
        ev[1] * 0.5 * vel[1] + std::abs(ev[1]) * (1.50 * Pyy - 1.0) - Pyy + 1.0;
    const amrex::Real Phiz =
        ev[2] * 0.5 * vel[2] + std::abs(ev[2]) * (1.50 * Pzz - 1.0) - Pzz + 1.0;

    eq = rho * Phix * Phiy * Phiz;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_extended_equilibrium_value_D3Q27(
    const amrex::Real& rho,
    const amrex::RealVect& vel,
    const amrex::Real& Pxx,
    const amrex::Real& Pyy,
    const amrex::Real& Pzz,
    const amrex::Real& mesh_speed,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    // Warning! Achtung! Strictly for D3Q27 only
    const amrex::Real Phix =
        ev[0] * 0.5 * vel[0] + std::abs(ev[0]) * (1.50 * Pxx - 1.0) - Pxx + 1.0;
    const amrex::Real Phiy =
        ev[1] * 0.5 * vel[1] + std::abs(ev[1]) * (1.50 * Pyy - 1.0) - Pyy + 1.0;
    const amrex::Real Phiz =
        ev[2] * 0.5 * vel[2] + std::abs(ev[2]) * (1.50 * Pzz - 1.0) - Pzz + 1.0;

    eq = rho * Phix * Phiy * Phiz;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_extended_gradExpansion_generic(
    const amrex::Real& rho,
    const amrex::RealVect& momentum,
    const amrex::Real& Pxx,
    const amrex::Real& Pyy,
    const amrex::Real& Pzz,
    const amrex::Real& Pxy,
    const amrex::Real& Pxz,
    const amrex::Real& Pyz,
    const amrex::Real& mesh_speed,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    const amrex::Real& theta0,
    const amrex::RealVect& frameVelocity,
    const amrex::Real& s,
    amrex::Real& f)
{
    // Note: This is lattice independent, as long as the correct lattice
    // temeprature theta0 is passed. theta0 = \sum_i w_i cx_i * cx_i

    const amrex::Real jx = momentum[0];
    const amrex::Real jy = momentum[1];
    const amrex::Real jz = momentum[2];

    // Note: ux, uy, uz are not flow velocities.
    // s is a special scaling factor not related to flow temperature
    // Unless absolutely confident about what they do, please use
    // ux=uy=uz=0 and s=1, always. Thank you for your understanding.
    // This is a special grad exapnsion with some future proofing.
    // Reference (A.93) of https://doi.org/10.3929/ethz-b-000607045
    const amrex::Real ux = frameVelocity[0];
    const amrex::Real uy = frameVelocity[1];
    const amrex::Real uz = frameVelocity[2];

    const amrex::Real stheta0 = s * theta0;
    const amrex::Real oneBystheta0 = 1.0 / stheta0;
    const amrex::Real a1x = ((jx)-rho * ux) * oneBystheta0;
    const amrex::Real a1y = ((jy)-rho * uy) * oneBystheta0;
    const amrex::Real a1z = ((jz)-rho * uz) * oneBystheta0;

    const amrex::Real a2xx = ((Pxx)-rho * s * stheta0 - rho * ux * ux -
                              ux * ((jx)-rho * ux) - ux * ((jx)-rho * ux)) *
                             oneBystheta0 * oneBystheta0;
    const amrex::Real a2yy = ((Pyy)-rho * s * stheta0 - rho * uy * uy -
                              uy * ((jy)-rho * uy) - uy * ((jy)-rho * uy)) *
                             oneBystheta0 * oneBystheta0;
    const amrex::Real a2zz = ((Pzz)-rho * s * stheta0 - rho * uz * uz -
                              uz * ((jz)-rho * uz) - uz * ((jz)-rho * uz)) *
                             oneBystheta0 * oneBystheta0;
    const amrex::Real a2xy = ((Pxy)-0 - rho * ux * uy - ux * ((jy)-rho * uy) -
                              uy * ((jx)-rho * ux)) *
                             oneBystheta0 * oneBystheta0;
    const amrex::Real a2xz = ((Pxz)-0 - rho * ux * uz - ux * ((jz)-rho * uz) -
                              uz * ((jx)-rho * ux)) *
                             oneBystheta0 * oneBystheta0;
    const amrex::Real a2yz = ((Pyz)-0 - rho * uy * uz - uy * ((jz)-rho * uz) -
                              uz * ((jy)-rho * uy)) *
                             oneBystheta0 * oneBystheta0;

    f = rho + (a1x * ev[0] + a1y * ev[1] + a1z * ev[2]);
    f += 0.5 *
         ((ev[0] * ev[0] - theta0) * a2xx + (ev[1] * ev[1] - theta0) * a2yy +
          (ev[2] * ev[2] - theta0) * a2zz + 2.0 * (ev[0] * ev[1] - 0) * a2xy +
          2.0 * (ev[0] * ev[2] - 0) * a2xz + 2.0 * (ev[1] * ev[2] - 0) * a2yz);
    f *= wt;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_population_zero(
        amrex::Real& eq)
{
    eq = 0.0;
}

template <class dType>
dType getEnergy(
    const dType T_,
    const dType rho_,
    const dType ux_,
    const dType uy_,
    const dType uz_,
    const dType Cv_)
{
    return rho_ * (2.0 * Cv_ * T_ + (ux_ * ux_ + uy_ * uy_ + uz_ * uz_));
}

template <class dType>
dType getTemperature(
    const dType twoRhoE_,
    const dType rho_,
    const dType ux_,
    const dType uy_,
    const dType uz_,
    const dType Cv_)
{
    return (0.50 / Cv_) *
           ((twoRhoE_ / rho_) - (ux_ * ux_ + uy_ * uy_ + uz_ * uz_));
}

template <class dType>
void getEquilibriumMoments(
    const dType rho,
    const dType ux,
    const dType uy,
    const dType uz,
    const dType totalEnergy,
    const dType Cv,
    const dType R,
    dType& qxEq,
    dType& qyEq,
    dType& qzEq,
    dType& RxxEq,
    dType& RyyEq,
    dType& RzzEq,
    dType& RxyEq,
    dType& RxzEq,
    dType& RyzEq)
{

    dType energy = totalEnergy / (2.0 * rho);
    dType temperature = getTemperature(totalEnergy, rho, ux, uy, uz, Cv);
    dType P = rho * R * temperature;
    dType H = energy + (P / rho);

    qxEq = 2.0 * rho * ux * H;
    qyEq = 2.0 * rho * uy * H;
    qzEq = 2.0 * rho * uz * H;
    RxxEq = 2.0 * rho * ux * ux * (H + (P / rho)) + 2.0 * P * H;
    RyyEq = 2.0 * rho * uy * uy * (H + (P / rho)) + 2.0 * P * H;
    RzzEq = 2.0 * rho * uz * uz * (H + (P / rho)) + 2.0 * P * H;
    RxyEq = 2.0 * rho * ux * uy * (H + (P / rho)) + 0;
    RxzEq = 2.0 * rho * ux * uz * (H + (P / rho)) + 0;
    RyzEq = 2.0 * rho * uy * uz * (H + (P / rho)) + 0;
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

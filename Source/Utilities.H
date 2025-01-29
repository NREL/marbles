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

AMREX_GPU_DEVICE AMREX_INLINE void set_equilibrium_value_d3_q27(
    const amrex::Real& rho,
    const amrex::RealVect& vel,
    const amrex::Real& Rtemperature,
    const amrex::Real& /*mesh_speed*/,
    const amrex::Real& /*wt*/,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    // Warning! Achtung! Strictly for D3Q27 only
    // This entropic equilibrium is very stable. Relaxes all 27 moments.
    // Unfortunately D3Q27 lattice has only 24 independent moments. For the
    // standard lattice, e.g. D3Q27, c_x^3 = c_x , think (-1)^3 = -1. So Qxxx,
    // Qyyy and Qzzz are always wrong. We correct them in pxx, pyy and pzz
    // during collision.

    const amrex::Real pxx = vel[0] * vel[0] + Rtemperature;
    const amrex::Real pyy = vel[1] * vel[1] + Rtemperature;
    amrex::Real pzz = 0.0;
    if (AMREX_SPACEDIM == 3) {
        pzz = vel[2] * vel[2] + Rtemperature;
    }

    const amrex::Real phix =
        ev[0] * 0.5 * vel[0] + std::abs(ev[0]) * (1.50 * pxx - 1.0) - pxx + 1.0;
    const amrex::Real phiy =
        ev[1] * 0.5 * vel[1] + std::abs(ev[1]) * (1.50 * pyy - 1.0) - pyy + 1.0;
    amrex::Real phiz = 1.0;
    if (AMREX_SPACEDIM == 3) {
        phiz = ev[2] * 0.5 * vel[2] + std::abs(ev[2]) * (1.50 * pzz - 1.0) -
               pzz + 1.0;
    }

    eq = rho * phix * phiy * phiz;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_extended_equilibrium_value_d3_q27(
    const amrex::Real& rho,
    const amrex::RealVect& vel,
    const amrex::Real& pxx,
    const amrex::Real& pyy,
    const amrex::Real& pzz,
    const amrex::Real& /*mesh_speed*/,
    const amrex::Real& /*wt*/,
    const amrex::IntVect& ev,
    amrex::Real& eq)
{
    // Warning! Achtung! Strictly for D3Q27 only
    const amrex::Real phix =
        ev[0] * 0.5 * vel[0] + std::abs(ev[0]) * (1.50 * pxx - 1.0) - pxx + 1.0;
    const amrex::Real phiy =
        ev[1] * 0.5 * vel[1] + std::abs(ev[1]) * (1.50 * pyy - 1.0) - pyy + 1.0;
    amrex::Real phiz = 1.0;
    if (AMREX_SPACEDIM == 3) {
        phiz = ev[2] * 0.5 * vel[2] + std::abs(ev[2]) * (1.50 * pzz - 1.0) -
               pzz + 1.0;
    }

    eq = rho * phix * phiy * phiz;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_extended_grad_expansion_generic(
    const amrex::Real& rho,
    const amrex::RealVect& momentum,
    const amrex::GpuArray<amrex::Real, 6>& fluxOfMomentum,
    /*const amrex::Real& pxx,
    const amrex::Real& pyy,
    const amrex::Real& pzz,
    const amrex::Real& pxy,
    const amrex::Real& pxz,
    const amrex::Real& pyz,*/
    const amrex::Real& /*mesh_speed*/,
    const amrex::Real& wt,
    const amrex::IntVect& ev,
    const amrex::Real& theta0,
    const amrex::RealVect& frameVelocity,
    const amrex::Real& s,
    amrex::Real& f)
{
    // Note: This is lattice independent, as long as the correct lattice
    // temperature theta0 is passed. theta0 = \sum_i w_i cx_i * cx_i

    const amrex::Real jx = momentum[0];
    const amrex::Real jy = momentum[1];
    // amrex::Real jz = 0.0;
    AMREX_3D_ONLY(const amrex::Real jz = momentum[2]);

    // Note: ux, uy, uz are not flow velocities.
    // s is a special scaling factor not related to flow temperature
    // Unless absolutely confident about what they do, please use
    // ux=uy=uz=0 and s=1, always. Thank you for your understanding.
    // This is a special grad expansion with some future proofing.
    // Reference (A.93) of https://doi.org/10.3929/ethz-b-000607045
    const amrex::Real ux = frameVelocity[0];
    const amrex::Real uy = frameVelocity[1];
    AMREX_3D_ONLY(const amrex::Real uz = frameVelocity[2]);

    const amrex::Real AMREX_D_DECL(
        pxx = fluxOfMomentum[0], pyy = fluxOfMomentum[1],
        pzz = fluxOfMomentum[2]);

    AMREX_2D_ONLY(const amrex::Real pxy = fluxOfMomentum[3]);

    AMREX_3D_ONLY(const amrex::Real pxy = fluxOfMomentum[3]);
    AMREX_3D_ONLY(const amrex::Real pxz = fluxOfMomentum[4]);
    AMREX_3D_ONLY(const amrex::Real pyz = fluxOfMomentum[5]);

    const amrex::Real stheta0 = s * theta0;
    const amrex::Real one_bystheta0 = 1.0 / stheta0;
    const amrex::Real AMREX_D_DECL(
        a1x = ((jx)-rho * ux) * one_bystheta0,
        a1y = ((jy)-rho * uy) * one_bystheta0,
        a1z = ((jz)-rho * uz) * one_bystheta0);

    const amrex::Real AMREX_D_DECL(
        a2xx = ((pxx)-rho * s * stheta0 - rho * ux * ux - ux * ((jx)-rho * ux) -
                ux * ((jx)-rho * ux)) *
               one_bystheta0 * one_bystheta0,
        a2yy = ((pyy)-rho * s * stheta0 - rho * uy * uy - uy * ((jy)-rho * uy) -
                uy * ((jy)-rho * uy)) *
               one_bystheta0 * one_bystheta0,
        a2zz = ((pzz)-rho * s * stheta0 - rho * uz * uz - uz * ((jz)-rho * uz) -
                uz * ((jz)-rho * uz)) *
               one_bystheta0 * one_bystheta0);

    AMREX_2D_ONLY(
        const amrex::Real a2xy = ((pxy)-0 - rho * ux * uy -
                                  ux * ((jy)-rho * uy) - uy * ((jx)-rho * ux)) *
                                 one_bystheta0 * one_bystheta0);

    AMREX_3D_ONLY(
        const amrex::Real a2xy = ((pxy)-0 - rho * ux * uy -
                                  ux * ((jy)-rho * uy) - uy * ((jx)-rho * ux)) *
                                 one_bystheta0 * one_bystheta0);

    AMREX_3D_ONLY(const amrex::Real a2xz =
                      ((pxz)-0 - rho * ux * uz - ux * ((jz)-rho * uz) -
                       uz * ((jx)-rho * ux)) *
                      one_bystheta0 * one_bystheta0;);

    AMREX_3D_ONLY(const amrex::Real a2yz =
                      ((pyz)-0 - rho * uy * uz - uy * ((jz)-rho * uz) -
                       uz * ((jy)-rho * uy)) *
                      one_bystheta0 * one_bystheta0;);

    f = rho + AMREX_D_TERM(a1x * ev[0], +a1y * ev[1], +a1z * ev[2]);

    f +=
        0.5 *
        (AMREX_D_TERM(
            (ev[0] * ev[0] - theta0) * a2xx,

            +(ev[1] * ev[1] - theta0) * a2yy + 2.0 * (ev[0] * ev[1] - 0) * a2xy,

            +(ev[2] * ev[2] - theta0) * a2zz +
                2.0 * (ev[0] * ev[2] - 0) * a2xz +
                2.0 * (ev[1] * ev[2] - 0) * a2yz));

    f *= wt;
}

AMREX_GPU_DEVICE AMREX_INLINE void set_population_zero(amrex::Real& eq)
{
    eq = 0.0;
}

template <class dType>
AMREX_GPU_DEVICE AMREX_INLINE dType get_energy(
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
AMREX_GPU_DEVICE AMREX_INLINE dType get_energy(
    const dType T_,
    const dType rho_,
    const amrex::RealVect& vel_,
    const dType Cv_)
{
    return rho_ * (2.0 * Cv_ * T_ + AMREX_D_TERM(
                                        vel_[0] * vel_[0], +vel_[1] * vel_[1],
                                        +vel_[2] * vel_[2]));
}

template <class dType>
AMREX_GPU_DEVICE AMREX_INLINE dType get_temperature(
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
AMREX_GPU_DEVICE AMREX_INLINE dType get_temperature(
    const dType twoRhoE_,
    const dType rho_,
    const amrex::RealVect& vel_,
    const dType Cv_)
{
    return (0.50 / Cv_) *
           ((twoRhoE_ / rho_) -
            AMREX_D_TERM(
                vel_[0] * vel_[0], +vel_[1] * vel_[1], +vel_[2] * vel_[2]));
}

template <class dType>
AMREX_GPU_DEVICE AMREX_INLINE void get_equilibrium_moments(
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
    dType p = rho * R * temperature;
    dType h = energy + (p / rho);

    qxEq = 2.0 * rho * ux * h;
    qyEq = 2.0 * rho * uy * h;
    qzEq = 2.0 * rho * uz * h;
    RxxEq = 2.0 * rho * ux * ux * (h + (p / rho)) + 2.0 * p * h;
    RyyEq = 2.0 * rho * uy * uy * (h + (p / rho)) + 2.0 * p * h;
    RzzEq = 2.0 * rho * uz * uz * (h + (p / rho)) + 2.0 * p * h;
    RxyEq = 2.0 * rho * ux * uy * (h + (p / rho)) + 0;
    RxzEq = 2.0 * rho * ux * uz * (h + (p / rho)) + 0;
    RyzEq = 2.0 * rho * uy * uz * (h + (p / rho)) + 0;
}

template <class dType>
AMREX_GPU_DEVICE AMREX_INLINE void get_equilibrium_moments(
    const dType rho,
    const amrex::RealVect& vel_,
    const dType totalEnergy,
    const dType Cv,
    const dType R,
    amrex::RealVect& qEq_,
    dType& RxxEq,
    dType& RyyEq,
    dType& RzzEq,
    dType& RxyEq,
    dType& RxzEq,
    dType& RyzEq)
{

    dType energy = totalEnergy / (2.0 * rho);
    dType temperature = get_temperature(totalEnergy, rho, vel_, Cv);
    dType p = rho * R * temperature;
    dType h = energy + (p / rho);

    qEq_[0] = 2.0 * rho * vel_[0] * h;
    qEq_[1] = 2.0 * rho * vel_[1] * h;
    if (AMREX_SPACEDIM == 3) {
        qEq_[2] = 2.0 * rho * vel_[2] * h;
    }

    RxxEq = 2.0 * rho * vel_[0] * vel_[0] * (h + (p / rho)) + 2.0 * p * h;
    RyyEq = 2.0 * rho * vel_[1] * vel_[1] * (h + (p / rho)) + 2.0 * p * h;
    if (AMREX_SPACEDIM == 3) {
        RzzEq = 2.0 * rho * vel_[2] * vel_[2] * (h + (p / rho)) + 2.0 * p * h;
    }
    RxyEq = 2.0 * rho * vel_[0] * vel_[1] * (h + (p / rho)) + 0;
    if (AMREX_SPACEDIM == 3) {
        RxzEq = 2.0 * rho * vel_[0] * vel_[2] * (h + (p / rho)) + 0;
    }
    if (AMREX_SPACEDIM == 3) {
        RyzEq = 2.0 * rho * vel_[1] * vel_[2] * (h + (p / rho)) + 0;
    }
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

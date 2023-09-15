#ifndef STENCIL_H
#define STENCIL_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include "Constants.H"

namespace lbm::stencil {
struct Stencil
{
    Stencil() = default;

#if AMREX_SPACEDIM == 2
    const amrex::GpuArray<amrex::IntVect, constants::N_MICRO_STATES> evs = {
        amrex::IntVect(0, 0),   amrex::IntVect(1, 0),  amrex::IntVect(-1, 0),
        amrex::IntVect(0, 1),   amrex::IntVect(0, -1), amrex::IntVect(1, 1),
        amrex::IntVect(-1, -1), amrex::IntVect(1, -1), amrex::IntVect(-1, 1)};

    const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_dirs = {
        0, 2, 1, 4, 3, 6, 5, 8, 7};

    const amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> weights = {
        4.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
#else
    const amrex::GpuArray<amrex::IntVect, constants::N_MICRO_STATES> evs = {
        amrex::IntVect(0, 0, 0),    amrex::IntVect(1, 0, 0),
        amrex::IntVect(-1, 0, 0),   amrex::IntVect(0, 1, 0),
        amrex::IntVect(0, -1, 0),   amrex::IntVect(0, 0, 1),
        amrex::IntVect(0, 0, -1),   amrex::IntVect(1, 1, 0),
        amrex::IntVect(-1, -1, 0),  amrex::IntVect(1, -1, 0),
        amrex::IntVect(-1, 1, 0),   amrex::IntVect(1, 0, 1),
        amrex::IntVect(-1, 0, -1),  amrex::IntVect(1, 0, -1),
        amrex::IntVect(-1, 0, 1),   amrex::IntVect(0, 1, 1),
        amrex::IntVect(0, -1, -1),  amrex::IntVect(0, 1, -1),
        amrex::IntVect(0, -1, 1),   amrex::IntVect(1, 1, 1),
        amrex::IntVect(-1, -1, -1), amrex::IntVect(1, 1, -1),
        amrex::IntVect(-1, -1, 1),  amrex::IntVect(1, -1, 1),
        amrex::IntVect(-1, 1, -1),  amrex::IntVect(-1, 1, 1),
        amrex::IntVect(1, -1, -1)};

    const amrex::GpuArray<int, constants::N_MICRO_STATES> bounce_dirs = {
        0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
        13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

    const amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> weights = {
        8.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,
        2.0 / 27.0,  2.0 / 27.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
        1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
        1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 216.0,
        1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
        1.0 / 216.0, 1.0 / 216.0};
#endif
};

void check_stencil();
} // namespace lbm::stencil
#endif

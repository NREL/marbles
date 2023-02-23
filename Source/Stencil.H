#ifndef STENCIL_H
#define STENCIL_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include "Constants.H"

namespace lbm {

namespace stencil {
struct Stencil
{
    Stencil() = default;

    const amrex::GpuArray<amrex::IntVect, constants::n_micro_states> evs = {
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

    const amrex::GpuArray<int, constants::n_micro_states> bounce_dirs = {
        0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
        13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

    const amrex::GpuArray<amrex::Real, constants::n_micro_states> weights = {
        8.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,
        2.0 / 27.0,  2.0 / 27.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
        1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
        1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 216.0,
        1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
        1.0 / 216.0, 1.0 / 216.0};
};

void CheckStencil();
} // namespace stencil

} // namespace lbm
#endif

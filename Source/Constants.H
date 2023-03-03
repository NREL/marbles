#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <AMReX_REAL.H>

namespace lbm::constants {
static constexpr int N_MACRO_STATES = 5;
static constexpr int RHO_IDX = 0;
static constexpr int VELX_IDX = 1;
static constexpr int VELY_IDX = 2;
static constexpr int VELZ_IDX = 3;
static constexpr int VMAG_IDX = 4;

static constexpr int N_MICRO_STATES = 27;

static constexpr amrex::Real ROOT3 =
    static_cast<amrex::Real>(1.732050807568877293527446341);

static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

static constexpr amrex::Real SMALL_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

static constexpr amrex::Real LARGE_NUM =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);
} // namespace lbm::constants

#endif

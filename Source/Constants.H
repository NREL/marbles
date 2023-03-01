#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <AMReX_REAL.H>

namespace lbm::constants {
static constexpr int n_macro_states = 5;
static constexpr int rho_idx = 0;
static constexpr int velx_idx = 1;
static constexpr int vely_idx = 2;
static constexpr int velz_idx = 3;
static constexpr int vmag_idx = 4;

static constexpr int n_micro_states = 27;

static constexpr amrex::Real root3 = 1.732050807568877293527446341;

static constexpr amrex::Real low_num = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

static constexpr amrex::Real small_num = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

static constexpr amrex::Real large_num =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);
} // namespace lbm::constants

#endif

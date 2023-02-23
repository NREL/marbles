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

} // namespace lbm::constants

#endif

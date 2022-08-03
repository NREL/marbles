#ifndef UTILITIES_H
#define UTILITIES_H
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include "constants.H"

namespace lbm {
void set_body_state(
  amrex::MultiFab& mf, const amrex::iMultiFab& is_fluid);
}
#endif

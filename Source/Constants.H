#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <AMReX_REAL.H>

namespace lbm::constants {
// static constexpr int N_MACRO_STATES = 9; //ns: was 5, added 1 for energy and
// 3 for Qcorrections
static constexpr int N_MACRO_STATES =
    19; // ns: was 9, added 6 for pressure tensor and 3 for heat flux
static constexpr int RHO_IDX = 0;
static constexpr int VELX_IDX = 1;
static constexpr int VELY_IDX = 2;
static constexpr int VELZ_IDX = 3;
static constexpr int VMAG_IDX = 4;

static constexpr int TWO_RHO_E_IDX = 5; // ns: This is the total energy
static constexpr int Q_CORR_X_IDX =
    6; // ns: Qxxx deviation from maxwell boltzmann distribution
static constexpr int Q_CORR_Y_IDX =
    7; // ns: Qyyy deviation from maxwell boltzmann distribution
static constexpr int Q_CORR_Z_IDX =
    8; // ns: Qzzz deviation from maxwell boltzmann distribution

static constexpr int PXX_IDX =
    9; // ns: Contains non-equilibrium part of pressure also
static constexpr int PYY_IDX =
    10; // ns: Contains non-equilibrium part of pressure also
static constexpr int PZZ_IDX =
    11; // ns: Contains non-equilibrium part of pressure also
static constexpr int PXY_IDX =
    12; // ns: Contains non-equilibrium part of pressure also
static constexpr int PXZ_IDX =
    13; // ns: Contains non-equilibrium part of pressure also
static constexpr int PYZ_IDX =
    14; // ns: Contains non-equilibrium part of pressure also

static constexpr int QX_IDX =
    15; // ns: Contains non-equilibrium part of heat flux also
static constexpr int QY_IDX =
    16; // ns: Contains non-equilibrium part of heat flux also
static constexpr int QZ_IDX =
    17; // ns: Contains non-equilibrium part of heat flux also

static constexpr int TEMPERATURE_IDX = 18; // ns: Temperature

#if AMREX_SPACEDIM == 2
static constexpr int N_MICRO_STATES = 9;
#else
static constexpr int N_MICRO_STATES = 27;
#endif

static constexpr int N_DERIVED =
    7; // ns: was 4, added 3 more for holding corrections
static constexpr int VORTX_IDX = 0;
static constexpr int VORTY_IDX = 1;
static constexpr int VORTZ_IDX = 2;
static constexpr int VORTM_IDX = 3;

static constexpr int D_Q_CORR_X_IDX =
    4; // ns: Actual correction term for Qxxx after gradient calculation of
       // Q_CORR_X_IDX
static constexpr int D_Q_CORR_Y_IDX =
    5; // ns: Actual correction term for Qyyy after gradient calculation of
       // Q_CORR_Y_IDX
static constexpr int D_Q_CORR_Z_IDX =
    6; // ns: Actual correction term for Qzzz after gradient calculation of
       // Q_CORR_Z_IDX

static constexpr int N_IS_FLUID = 2;

static constexpr amrex::Real ROOT3 =
    static_cast<amrex::Real>(1.732050807568877293527446341);

static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

static constexpr amrex::Real SMALL_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

static constexpr amrex::Real LARGE_NUM =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);

static constexpr int DATWIDTH = 24;
static constexpr int DATPRECISION = 16;

} // namespace lbm::constants

#endif

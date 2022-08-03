#ifndef LBM_PARTICLE_SPECS_H
#define LBM_PARTICLE_SPECS_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <constants.H>

namespace lbm {
class LBMspecs
{
public:
    amrex::Array<amrex::Real, AMREX_SPACEDIM> plo{AMREX_D_DECL(0.0, 0.0, 0.0)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> phi{AMREX_D_DECL(0.0, 0.0, 0.0)};
    amrex::Array<int, AMREX_SPACEDIM> ncells{AMREX_D_DECL(0, 0, 0)};
    amrex::Array<int, AMREX_SPACEDIM> periodic{AMREX_D_DECL(0, 0, 0)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> gravity{
        AMREX_D_DECL(0.0, -9.81, 0.0)};

    amrex::Real timestep = 0.0;
    int num_redist = 10;
    int max_grid_size = 8;
    amrex::Real final_time;
    amrex::Real write_output_time;
    amrex::Real screen_output_time;
    int max_steps = 1;

    amrex::Real CFL = 0.1;
    int print_diagnostics = 0;

    amrex::Vector<int> bclo;
    amrex::Vector<int> bchi;

    amrex::Real mesh_speed = 1.0;
    amrex::Real tau = 0.7;

    void read_lbm_specs()
    {
        bclo.resize(AMREX_SPACEDIM);
        bchi.resize(AMREX_SPACEDIM);

        amrex::ParmParse pp("lbm");
        pp.get("prob_lo", plo);
        pp.get("prob_hi", phi);
        pp.get("ncells", ncells);
        pp.get("is_it_periodic", periodic);
        pp.get("max_grid_size", max_grid_size);

        pp.get("max_steps", max_steps);
        pp.get("final_time", final_time);

        pp.query("gravity", gravity);

        write_output_time = final_time / 10.0;
        pp.query("write_output_time", write_output_time);

        screen_output_time = write_output_time / 10.0;
        pp.query("screen_output_time", screen_output_time);

        pp.query("CFL", CFL);

        // by default it is periodic
        bclo[0] = BC_SLIPWALL;
        bclo[1] = BC_SLIPWALL;
        bclo[2] = BC_SLIPWALL;

        bchi[0] = BC_SLIPWALL;
        bchi[1] = BC_SLIPWALL;
        bchi[2] = BC_SLIPWALL;

        pp.queryarr("bc_lower", bclo);
        pp.queryarr("bc_upper", bchi);

        pp.query("mesh_speed", mesh_speed);
        pp.query("tau", tau);
    }
};
} // namespace lbm
#endif

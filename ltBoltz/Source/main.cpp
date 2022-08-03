#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <typeinfo>
#include "constants.H"
#include "lbm_specs.H"
#include "EB.H"
#include "IO.H"
#include "Utilities.H"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        // Reading input files for the simulation
        lbm::LBMspecs specs;
        specs.read_lbm_specs(); // Read input file

        const int coord = 0; // cartesian
        amrex::RealBox real_box;
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
            real_box.setLo(n, specs.plo[n]);
            real_box.setHi(n, specs.phi[n]);
        }

        const amrex::IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
        const amrex::IntVect domain_hi(AMREX_D_DECL(
            specs.ncells[XDIR] - 1, specs.ncells[YDIR] - 1,
            specs.ncells[ZDIR] - 1));
        const amrex::Box domain(domain_lo, domain_hi);

        amrex::Geometry geom(domain, &real_box, coord, specs.periodic.data());

        amrex::BoxArray ba(domain);
        ba.maxSize(specs.max_grid_size);
        amrex::DistributionMapping dm(ba);

        const int max_level = 1;
        lbm::initialize_eb(geom, max_level);
        const amrex::Vector<int> eb_ng = {5, 5, 5};
        const auto factory = amrex::makeEBFabFactory(
            geom, ba, dm, eb_ng, amrex::EBSupport::basic);
        auto const& flags = factory->getMultiEBCellFlagFab();

        const int ng_cells = 1;
        amrex::MultiFab macrodata(
            ba, dm, NUM_MACRO_STATES, 0, amrex::MFInfo(), *factory);
        amrex::MultiFab f(
            ba, dm, NUM_MICRO_STATES, ng_cells, amrex::MFInfo(), *factory);
        amrex::MultiFab f_star(
            ba, dm, NUM_MICRO_STATES, f.nGrow() + 1, amrex::MFInfo(), *factory);
        amrex::iMultiFab is_fluid(ba, dm, 1, f.nGrow() + 1);
        amrex::MultiFab eq(
            ba, dm, NUM_MICRO_STATES, f.nGrow(), amrex::MFInfo(), *factory);
        macrodata.setVal(0.0);
        lbm::set_body_state(macrodata, is_fluid);

        lbm::IO io(ba, dm, geom, *factory);

        const amrex::GpuArray<int, NUM_MICRO_STATES* 3> e = {
            0, 0, 0,  1,  0,  0, -1, 0,  0,  0,  1,  0,  0, -1, 0,  0,  0,
            1, 0, 0,  -1, 1,  1, 0,  -1, -1, 0,  1,  -1, 0, -1, 1,  0,  1,
            0, 1, -1, 0,  -1, 1, 0,  -1, -1, 0,  1,  0,  1, 1,  0,  -1, -1,
            0, 1, -1, 0,  -1, 1, 1,  1,  1,  -1, -1, -1, 1, 1,  -1, -1, -1,
            1, 1, -1, 1,  -1, 1, -1, -1, 1,  1,  1,  -1, -1};

        const amrex::GpuArray<int, NUM_MICRO_STATES> bounce_dir = {
            0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
            13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

        const amrex::GpuArray<amrex::Real, NUM_MICRO_STATES> weight = {
            8.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,
            2.0 / 27.0,  2.0 / 27.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
            1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,
            1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 216.0,
            1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
            1.0 / 216.0, 1.0 / 216.0};

        // Initialize the fluid point mask and distributions
        // (outside time loop)
        f_star.setVal(0.0, f_star.nGrow());
        auto const& f_arrs = f.arrays();
        auto const& fs_arrs = f_star.arrays();
        auto const& is_fluid_arrs = is_fluid.arrays();
        auto const& flag_arrs = flags.const_arrays();
        amrex::ParallelFor(
          f, f.nGrowVect(), NUM_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                f_arrs[nbx](i, j, k, q) = weight[q];
            });
        amrex::ParallelFor(
          is_fluid, is_fluid.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                is_fluid_arrs[nbx](i, j, k) =
                    ((i < domain_lo[0] || i > domain_hi[0] ||
                      j < domain_lo[1] || j > domain_hi[1] ||
                      k < domain_lo[2] || k > domain_hi[2]) ||
                     (!flag_arrs[nbx](i, j, k).isRegular()))
                        ? 0
                        : 1;
            });
        amrex::Gpu::synchronize();
        f.FillBoundary();
        is_fluid.FillBoundary();

        amrex::Real time = 0.0;
        amrex::Real dt = 1.0;
        const int save_every_n = 10;
        const std::string save_dir = "output";

        int next_save_time = save_every_n;

        io.write(0, save_dir, time, macrodata, f, is_fluid);

        for (int step = 1; step <= specs.max_steps; step++) {

            amrex::Print()
                << "\n==============================================="
                   "==============================="
                << std::endl;
            amrex::Print() << "Step: " << step << " dt : " << dt
                           << " time: " << time << " to " << time + dt
                           << std::endl;

            // Stream the information to the neighbor particles
            amrex::ParallelFor(
              f, f.nGrowVect(), NUM_MICRO_STATES,
                [=] AMREX_GPU_DEVICE(
                    int nbx, int i, int j, int k, int q) noexcept {
                    const amrex::IntVect iv(i, j, k);
                    const amrex::IntVect ivn(
                        i + e[q * 3 + 0], j + e[q * 3 + 1], k + e[q * 3 + 2]);
                    if (is_fluid_arrs[nbx](iv) == 1) {
                        const auto f_arr = f_arrs[nbx];
                        const auto fs_arr = fs_arrs[nbx];
                        if (is_fluid_arrs[nbx](ivn) == 1) {
                            fs_arr(ivn, q) = f_arr(iv, q);
                        } else {
                            fs_arr(iv, bounce_dir[q]) = f_arr(iv, q);
                        }
                    }
                });
            amrex::Gpu::synchronize();

            amrex::MultiFab::Copy(f, f_star, 0, 0, NUM_MICRO_STATES, 0);

            auto const& md_arrs = macrodata.arrays();
            auto const& eq_arrs = eq.arrays();
            const amrex::Real tau = specs.tau;
            const amrex::Real mesh_speed = specs.mesh_speed;
            amrex::ParallelFor(
              f, macrodata.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                    const amrex::IntVect iv(i, j, k);
                    if (is_fluid_arrs[nbx](iv) == 1) {

                      const auto f_arr = f_arrs[nbx];
                      const auto md_arr = md_arrs[nbx];
                    const auto eq_arr = eq_arrs[nbx];

                    // Calculate the macro fluid properties from the
                    // distributions
                    amrex::Real rho_temp = 0.0, u_temp = 0.0, v_temp = 0.0,
                                w_temp = 0.0;
                    for (int q = 0; q < NUM_MICRO_STATES; q++) {
                        rho_temp += f_arr(iv, q);
                        u_temp += mesh_speed * e[q * 3 + 0] * f_arr(iv, q);
                        v_temp += mesh_speed * e[q * 3 + 1] * f_arr(iv, q);
                        w_temp += mesh_speed * e[q * 3 + 2] * f_arr(iv, q);
                    }
                    u_temp /= rho_temp;
                    v_temp /= rho_temp;
                    w_temp /= rho_temp;

                    if (k >= domain_hi[2] - 1) {
                        u_temp = 0.1;
                        v_temp = 0.0;
                        w_temp = 0.0;
                    }

                    md_arr(iv, RHO_INDEX) = rho_temp;
                    md_arr(iv, VELX_INDEX) = u_temp;
                    md_arr(iv, VELY_INDEX) = v_temp;
                    md_arr(iv, VELZ_INDEX) = w_temp;

                    // Calculate the equilibrium state
                    const amrex::Real umag2 =
                        u_temp * u_temp + v_temp * v_temp + w_temp * w_temp;
                    const amrex::Real c3 =
                        -1.5 * umag2 / (mesh_speed * mesh_speed);

                    for (int q = 0; q < NUM_MICRO_STATES; q++) {
                        const amrex::Real e_dot_u = e[q * 3 + 0] * u_temp +
                                                    e[q * 3 + 1] * v_temp +
                                                    e[q * 3 + 2] * w_temp;
                        const amrex::Real e_div_c = e_dot_u / mesh_speed;
                        const amrex::Real c1 = 3.0 * e_div_c;
                        const amrex::Real c2 = 4.5 * (e_div_c * e_div_c);

                        eq_arr(iv, q) =
                            rho_temp * weight[q] * (1.0 + c1 + c2 + c3);
                    }

                    // Relax the particles toward the equilibrium state
                    for (int q = 0; q < NUM_MICRO_STATES; q++) {
                        f_arr(iv, q) -=
                            1.0 / tau * (f_arr(iv, q) - eq_arr(iv, q));
                    }
                    }
                });
            amrex::Gpu::synchronize();

            f.FillBoundary();

            time += dt;
            if (step >= next_save_time) {
                io.write(step, save_dir, time, macrodata, f, is_fluid);

                next_save_time += save_every_n;
            }
        }
    }

    amrex::Finalize();
}

#include "LBM.H"

namespace lbm {
LBM::LBM()
{
    BL_PROFILE("LBM::LBM()");
    ReadParameters();

    AMREX_ASSERT(max_level == 0);

    int nlevs_max = max_level + 1;
    initialize_eb(Geom(maxLevel()), maxLevel());

    const size_t n_zero = 2;
    lbm_varnames.push_back("rho");
    lbm_varnames.push_back("vel_x");
    lbm_varnames.push_back("vel_y");
    lbm_varnames.push_back("vel_z");
    for (int q = 0; q < NUM_MICRO_STATES; q++) {
        const auto num_str = std::to_string(q);
        const auto zero_padded_str =
            std::string(n_zero - std::min(n_zero, num_str.length()), '0') +
            num_str;
        lbm_varnames.push_back("f_" + zero_padded_str);
    }
    lbm_varnames.push_back("is_fluid");

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    macrodata.resize(nlevs_max);
    f_.resize(nlevs_max);
    eq.resize(nlevs_max);
    is_fluid.resize(nlevs_max);
    plt_mf.resize(nlevs_max);

    m_factory.resize(nlevs_max);

    // FIXME all this
    bcs.resize(NUM_MACRO_STATES); // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // lo-side BCs
        if (bc_lo[idim] ==
                amrex::BCType::int_dir || // periodic uses "internal Dirichlet"
            bc_lo[idim] ==
                amrex::BCType::foextrap || // first-order extrapolation
            bc_lo[idim] == amrex::BCType::ext_dir) {
            for (int comp = 0; comp < NUM_MACRO_STATES; comp++) {
                bcs[comp].setLo(idim, bc_lo[idim]);
            }
        } else if (bc_lo[idim] == BC_NOSLIPWALL) {
            for (int comp = 0; comp < NUM_MACRO_STATES; comp++) {
                bcs[comp].setLo(idim, amrex::BCType::reflect_even);
            }
            // bcs[RHOU_INDX].setLo(idim, amrex::BCType::reflect_odd);
            // bcs[RHOV_INDX].setLo(idim, amrex::BCType::reflect_odd);
            // bcs[RHOW_INDX].setLo(idim, amrex::BCType::reflect_odd);
            // bcs[VELX_INDX].setLo(idim, amrex::BCType::reflect_odd);
            // bcs[VELY_INDX].setLo(idim, amrex::BCType::reflect_odd);
            // bcs[VELZ_INDX].setLo(idim, amrex::BCType::reflect_odd);
        } else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] ==
                amrex::BCType::int_dir || // periodic uses "internal Dirichlet"
            bc_hi[idim] ==
                amrex::BCType::foextrap || // first-order extrapolation
            bc_hi[idim] == amrex::BCType::ext_dir) {
            for (int comp = 0; comp < NUM_MACRO_STATES; comp++) {
                bcs[comp].setHi(idim, bc_hi[idim]);
            }
        } else if (bc_hi[idim] == BC_NOSLIPWALL) {
            for (int comp = 0; comp < NUM_MACRO_STATES; comp++) {
                bcs[comp].setHi(idim, amrex::BCType::reflect_even);
            }
            // bcs[RHOU_INDX].setHi(idim, amrex::BCType::reflect_odd);
            // bcs[RHOV_INDX].setHi(idim, amrex::BCType::reflect_odd);
            // bcs[RHOW_INDX].setHi(idim, amrex::BCType::reflect_odd);
            // bcs[VELX_INDX].setHi(idim, amrex::BCType::reflect_odd);
            // bcs[VELY_INDX].setHi(idim, amrex::BCType::reflect_odd);
            // bcs[VELZ_INDX].setHi(idim, amrex::BCType::reflect_odd);
        } else {
            amrex::Abort("Invalid bc_hi");
        }
    }

    flux_reg.resize(nlevs_max + 1);
}

LBM::~LBM() {}

void LBM::InitData()
{
    BL_PROFILE("LBM::InitData()");
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const amrex::Real time = 0.0;
        InitFromScratch(time);
        // AverageDown(); // FIXME

        if (chk_int > 0) {
            amrex::Abort("No support for writing checkpoint yet.");
            // WriteCheckpointFile();
        }

    } else {
        // restart from a checkpoint
        amrex::Abort("No support for restarting from checkpoint yet.");
        // ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
}

void LBM::ReadParameters()
{
    BL_PROFILE("LBM::ReadParameters()");

    {
        amrex::ParmParse pp;
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        amrex::ParmParse pp("amr");

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart", restart_chkfile);
    }

    {
        amrex::ParmParse pp("lbm");
        pp.queryarr("lo_bc", bc_lo, 0, AMREX_SPACEDIM);
        pp.queryarr("hi_bc", bc_hi, 0, AMREX_SPACEDIM);
        pp.query("do_reflux", do_reflux);
        pp.query("tau", m_tau);
        pp.query("mesh_speed", m_mesh_speed);
    }
}

void LBM::Evolve()
{
    BL_PROFILE("LBM::Evolve()");

    amrex::Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step) {
        ComputeDt();

        amrex::Print() << "\n==============================================="
                          "==============================="
                       << std::endl;
        amrex::Print() << "Step: " << step << " dt : " << dt[0]
                       << " time: " << cur_time << " to " << cur_time + dt[0]
                       << std::endl;

        timeStep(0, cur_time, 1);

        cur_time += dt[0];

        // sync up time
        for (int lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step + 1) % plot_int == 0) {
            last_plot_file_step = step + 1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step + 1) % chk_int == 0) {
            amrex::Abort("Chk not supported yet.");
            // FIXME WriteCheckpointFile();
        }

        if (cur_time >= stop_time - 1.e-6 * dt[0]) break;
    }
    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }
}

// advance a level by dt
// includes a recursive call for finer levels
void LBM::timeStep(const int lev, const amrex::Real time, const int iteration)
{
    BL_PROFILE("LBM::timeStep()");
    if (regrid_int > 0) // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static amrex::Vector<int> last_regrid_step(max_level + 1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) {
            if (istep[lev] % regrid_int == 0) {
                // regrid could add newly refine levels (if finest_level <
                // max_level) so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest + 1; k <= finest_level; ++k) {
                    dt[k] = dt[k - 1] / MaxRefRatio(k - 1);
                }
            }
        }
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] + 1
                       << "] ";
        amrex::Print() << "ADVANCE with time = " << t_new[lev]
                       << " dt = " << dt[lev] << std::endl;
    }

    // advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells"
                       << std::endl;
    }

    if (lev < finest_level) {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev + 1]; ++i) {
            timeStep(lev + 1, time + (i - 1) * dt[lev + 1], i);
        }

        if (do_reflux) {
            // update lev based on coarse-fine flux mismatch
            flux_reg[lev + 1]->Reflux(
                macrodata[lev], 1.0, 0, 0, macrodata[lev].nComp(), geom[lev]);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }
}

void LBM::Advance(
    const int lev,
    const amrex::Real time,
    const amrex::Real dt_lev,
    const int iteration,
    const int ncycle)
{
    BL_PROFILE("LBM::Advance()");

    t_old[lev] = t_new[lev]; // old time is now current time (time)
    t_new[lev] += dt_lev;    // new time is ahead

    Stream(lev);

    FToMacrodata(lev);

    MacrodataToEquilibrium(lev);

    RelaxFToEquilibrium(lev);

    f_[lev].FillBoundary(Geom(lev).periodicity()); // FIXME check if needed
}

// Stream the information to the neighbor particles
void LBM::Stream(const int lev)
{
    BL_PROFILE("LBM::Stream()");
    amrex::MultiFab f_star(
        boxArray(lev), DistributionMap(lev), NUM_MICRO_STATES,
        f_[lev].nGrow() + 1, amrex::MFInfo(), *(m_factory[lev]));
    f_star.setVal(0.0);

    auto const& fs_arrs = f_star.arrays();
    auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
    auto const& f_arrs = f_[lev].const_arrays();

    amrex::ParallelFor(
        f_[lev], f_[lev].nGrowVect(), NUM_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
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

    amrex::MultiFab::Copy(f_[lev], f_star, 0, 0, NUM_MICRO_STATES, 0);
}

// convert equilibrium to macrodata
void LBM::MacrodataToEquilibrium(const int lev)
{
    BL_PROFILE("LBM::MacrodataToEquilibrium()");
    auto const& md_arrs = macrodata[lev].const_arrays();
    auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
    auto const& eq_arrs = eq[lev].arrays();
    const amrex::Real mesh_speed = m_mesh_speed;

    amrex::ParallelFor(
        eq[lev], macrodata[lev].nGrowVect(), NUM_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            const amrex::IntVect iv(i, j, k);
            if (is_fluid_arrs[nbx](iv) == 1) {

                const auto md_arr = md_arrs[nbx];
                const auto eq_arr = eq_arrs[nbx];

                const amrex::Real rho = md_arr(iv, RHO_INDEX);
                const amrex::Real u = md_arr(iv, VELX_INDEX);
                const amrex::Real v = md_arr(iv, VELY_INDEX);
                const amrex::Real w = md_arr(iv, VELZ_INDEX);

                const amrex::Real umag2 = u * u + v * v + w * w;
                const amrex::Real c3 = -1.5 * umag2 / (mesh_speed * mesh_speed);

                const amrex::Real e_dot_u =
                    e[q * 3 + 0] * u + e[q * 3 + 1] * v + e[q * 3 + 2] * w;
                const amrex::Real e_div_c = e_dot_u / mesh_speed;
                const amrex::Real c1 = 3.0 * e_div_c;
                const amrex::Real c2 = 4.5 * (e_div_c * e_div_c);

                eq_arr(iv, q) = rho * weight[q] * (1.0 + c1 + c2 + c3);
            }
        });
    amrex::Gpu::synchronize();
}

// Relax the particles toward the equilibrium state
void LBM::RelaxFToEquilibrium(const int lev)
{
    BL_PROFILE("LBM::RelaxFToEquilibrium()");
    auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
    auto const& eq_arrs = eq[lev].const_arrays();
    auto const& f_arrs = f_[lev].arrays();
    const amrex::Real tau = m_tau;
    amrex::ParallelFor(
        f_[lev], macrodata[lev].nGrowVect(), NUM_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            const amrex::IntVect iv(i, j, k);
            if (is_fluid_arrs[nbx](iv) == 1) {
                const auto f_arr = f_arrs[nbx];
                const auto eq_arr = eq_arrs[nbx];
                f_arr(iv, q) -= 1.0 / tau * (f_arr(iv, q) - eq_arr(iv, q));
            }
        });
    amrex::Gpu::synchronize();
}

// calculate the macro fluid properties from the distributions
void LBM::FToMacrodata(const int lev)
{
    BL_PROFILE("LBM::FToMacrodata()");
    auto const& md_arrs = macrodata[lev].arrays();
    auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
    auto const& f_arrs = f_[lev].const_arrays();
    const amrex::Real mesh_speed = m_mesh_speed;
    const int* domhi = Geom(lev).Domain().hiVect();

    amrex::ParallelFor(
        macrodata[lev], macrodata[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::IntVect iv(i, j, k);
            if (is_fluid_arrs[nbx](iv) == 1) {

                const auto f_arr = f_arrs[nbx];
                const auto md_arr = md_arrs[nbx];

                amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0;
                for (int q = 0; q < NUM_MICRO_STATES; q++) {
                    rho += f_arr(iv, q);
                    u += mesh_speed * e[q * 3 + 0] * f_arr(iv, q);
                    v += mesh_speed * e[q * 3 + 1] * f_arr(iv, q);
                    w += mesh_speed * e[q * 3 + 2] * f_arr(iv, q);
                }
                u /= rho;
                v /= rho;
                w /= rho;

                if (k >= domhi[2] - 1) {
                    u = 0.1;
                    v = 0.0;
                    w = 0.0;
                }

                md_arr(iv, RHO_INDEX) = rho;
                md_arr(iv, VELX_INDEX) = u;
                md_arr(iv, VELY_INDEX) = v;
                md_arr(iv, VELZ_INDEX) = w;
            }
        });
    amrex::Gpu::synchronize();
}

// a wrapper for EstTimeStep
void LBM::ComputeDt()
{
    BL_PROFILE("LBM::ComputeDt()");
    amrex::Vector<amrex::Real> dt_tmp(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = EstTimeStep(lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr amrex::Real change_max = 1.1;
    amrex::Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max * dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor * dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const amrex::Real eps = 1.e-3 * dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev - 1] / nsubsteps[lev];
    }
}

// compute dt
amrex::Real LBM::EstTimeStep(const int lev)
{
    BL_PROFILE("LBM::EstTimeStep()");
    return 1.0;
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
void LBM::MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::MakeNewLevelFromCoarse()");
    const int ncomp = macrodata[lev - 1].nComp();
    const int nghost = macrodata[lev - 1].nGrow();

    macrodata[lev].define(ba, dm, ncomp, nghost);
    // FIXME, need the others here?

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(
            new amrex::FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, macrodata[lev], 0, ncomp);
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void LBM::MakeNewLevelFromScratch(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::MakeNewLevelFromScratch()");
    const int nghost = 1;

    m_factory[lev] = amrex::makeEBFabFactory(
        Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);

    macrodata[lev].define(
        ba, dm, NUM_MACRO_STATES, 0, amrex::MFInfo(), *(m_factory[lev]));
    f_[lev].define(
        ba, dm, NUM_MICRO_STATES, nghost, amrex::MFInfo(), *(m_factory[lev]));
    is_fluid[lev].define(ba, dm, 1, f_[lev].nGrow() + 1);
    eq[lev].define(
        ba, dm, NUM_MICRO_STATES, f_[lev].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        amrex::Abort("Not supported yet");
        flux_reg[lev].reset(new amrex::FluxRegister(
            ba, dm, refRatio(lev - 1), lev, macrodata[lev].nComp()));
    }

    amrex::Real cur_time = t_new[lev];

    // Initialize the data
    macrodata[lev].setVal(0.0);

    const int* domlo = Geom(lev).Domain().loVect();
    const int* domhi = Geom(lev).Domain().hiVect();

    auto const& f_arrs = f_[lev].arrays();
    auto const& is_fluid_arrs = is_fluid[lev].arrays();
    auto factory =
        static_cast<amrex::EBFArrayBoxFactory*>(m_factory[lev].get());

    auto const& flags = factory->getMultiEBCellFlagFab();
    auto const& flag_arrs = flags.const_arrays();
    amrex::ParallelFor(
        f_[lev], f_[lev].nGrowVect(), NUM_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            f_arrs[nbx](i, j, k, q) = weight[q];
        });
    amrex::ParallelFor(
        is_fluid[lev], is_fluid[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            is_fluid_arrs[nbx](i, j, k) =
                ((i < domlo[0] || i > domhi[0] || j < domlo[1] ||
                  j > domhi[1] || k < domlo[2] || k > domhi[2]) ||
                 (!flag_arrs[nbx](i, j, k).isRegular()))
                    ? 0
                    : 1;
        });
    amrex::Gpu::synchronize();
    f_[lev].FillBoundary(Geom(lev).periodicity());
    is_fluid[lev].FillBoundary(Geom(lev).periodicity());
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
void LBM::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::RemakeLevel()");
    amrex::MultiFab new_state(
        ba, dm, macrodata[lev].nComp(), macrodata[lev].nGrow());

    FillPatch(lev, time, new_state, 0, new_state.nComp());

    std::swap(new_state, macrodata[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new amrex::FluxRegister(
            ba, dm, refRatio(lev - 1), lev, macrodata[lev].nComp()));
    }
}

// Delete level data
void LBM::ClearLevel(int lev)
{
    BL_PROFILE("LBM::ClearLevel()");
    macrodata[lev].clear();
    f_[lev].clear();
    eq[lev].clear();
    is_fluid[lev].clear();
    flux_reg[lev].reset(nullptr);
}

// compute a new multifab by coping in phi from valid region and filling ghost
// cells works for single level and 2-level cases (fill fine grid ghost by
// interpolating from coarse)
void LBM::FillPatch(
    int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp)
{
    BL_PROFILE("LBM::FillPatch()");
    amrex::Abort("FillPatch not implemented");
    if (lev == 0) {
        amrex::Vector<amrex::MultiFab*> smf;
        amrex::Vector<amrex::Real> stime;
        GetData(0, time, smf, stime);

        amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(amrcore_fill_func);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill>> physbc(
            geom[lev], bcs, gpu_bndry_func);
        amrex::FillPatchSingleLevel(
            mf, time, smf, stime, 0, icomp, ncomp, geom[lev], physbc, 0);
    } else {
        amrex::Vector<amrex::MultiFab*> cmf, fmf;
        amrex::Vector<amrex::Real> ctime, ftime;
        GetData(lev - 1, time, cmf, ctime);
        GetData(lev, time, fmf, ftime);

        amrex::Interpolater* mapper = &amrex::cell_cons_interp;

        amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(amrcore_fill_func);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill>> cphysbc(
            geom[lev - 1], bcs, gpu_bndry_func);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill>> fphysbc(
            geom[lev], bcs, gpu_bndry_func);

        amrex::FillPatchTwoLevels(
            mf, time, cmf, ctime, fmf, ftime, 0, icomp, ncomp, geom[lev - 1],
            geom[lev], cphysbc, 0, fphysbc, 0, refRatio(lev - 1), mapper, bcs,
            0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void LBM::FillCoarsePatch(
    int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp)
{
    BL_PROFILE("LBM::FillCoarsePatch()");
    AMREX_ASSERT(lev > 0);

    amrex::Abort("FillCoarsePatch not implemented");
    amrex::Vector<amrex::MultiFab*> cmf;
    amrex::Vector<amrex::Real> ctime;
    GetData(lev - 1, time, cmf, ctime);
    amrex::Interpolater* mapper = &amrex::cell_cons_interp;

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(amrcore_fill_func);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill>> cphysbc(
        geom[lev - 1], bcs, gpu_bndry_func);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill>> fphysbc(
        geom[lev], bcs, gpu_bndry_func);

    amrex::InterpFromCoarseLevel(
        mf, time, *cmf[0], 0, icomp, ncomp, geom[lev - 1], geom[lev], cphysbc,
        0, fphysbc, 0, refRatio(lev - 1), mapper, bcs, 0);
}

// utility to copy in data from macrodata into another multifab
void LBM::GetData(
    int lev,
    amrex::Real time,
    amrex::Vector<amrex::MultiFab*>& data,
    amrex::Vector<amrex::Real>& datatime)
{
    BL_PROFILE("LBM::GetData()");
    amrex::Abort("GetData not implemented");
    data.clear();
    datatime.clear();

    const amrex::Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps) {
        data.push_back(&macrodata[lev]);
        datatime.push_back(t_new[lev]);
    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps) {
        data.push_back(&macrodata[lev]);
        datatime.push_back(t_old[lev]);
    } else {
        data.push_back(&macrodata[lev]);
        data.push_back(&macrodata[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void LBM::AverageDown()
{
    BL_PROFILE("LBM::AverageDown()");
    amrex::Abort("AverageDown not implemented");
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        amrex::average_down(
            macrodata[lev + 1], macrodata[lev], geom[lev + 1], geom[lev], 0,
            macrodata[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across
// multiple levels
void LBM::AverageDownTo(int crse_lev)
{
    BL_PROFILE("LBM::AverageDownTo()");
    amrex::Abort("AverageDownTo not implemented");
    amrex::average_down(
        macrodata[crse_lev + 1], macrodata[crse_lev], geom[crse_lev + 1],
        geom[crse_lev], 0, macrodata[crse_lev].nComp(), refRatio(crse_lev));
}

// tag cells for refinement
void LBM::ErrorEst(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
    BL_PROFILE("LBM::ErrorEst()");
}

amrex::Vector<std::string> LBM::PlotFileVarNames() const
{
    return lbm_varnames;
}

std::string LBM::PlotFileName(const int step) const
{
    return amrex::Concatenate(plot_file, step, 5);
}

// put together an array of multifabs for writing
amrex::Vector<const amrex::MultiFab*> LBM::PlotFileMF()
{
    amrex::Vector<const amrex::MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {

        plt_mf[i].define(
            boxArray(i), DistributionMap(i), PlotFileVarNames().size(), 0);
        int cnt = 0;
        amrex::MultiFab::Copy(
            plt_mf[i], macrodata[i], 0, cnt, macrodata[i].nComp(), 0);
        cnt += macrodata[i].nComp();
        amrex::MultiFab::Copy(plt_mf[i], f_[i], 0, cnt, f_[i].nComp(), 0);
        cnt += f_[i].nComp();
        auto const& is_fluid_arrs = is_fluid[i].const_arrays();
        auto const& plt_mf_arrs = plt_mf[i].arrays();
        amrex::ParallelFor(
            plt_mf[i], plt_mf[i].nGrowVect(), is_fluid[i].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                plt_mf_arrs[nbx](i, j, k, n + cnt) =
                    is_fluid_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();
        cnt += is_fluid[i].nComp();

        r.push_back(&plt_mf[i]);
    }
    return r;
}

void LBM::WritePlotFile()
{
    BL_PROFILE("LBM::WritePlotFile()");
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plot file " << plotfilename << " at time "
                   << t_new[0] << std::endl;

    amrex::WriteMultiLevelPlotfile(
        plotfilename, finest_level + 1, mf, varnames, Geom(), t_new[0], istep,
        refRatio());
}

} // namespace lbm

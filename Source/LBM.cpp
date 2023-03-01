#include "LBM.H"

namespace lbm {
LBM::LBM()
{
    BL_PROFILE("LBM::LBM()");
    ReadParameters();

    int nlevs_max = max_level + 1;
    initialize_eb(Geom(maxLevel()), maxLevel());

    macrodata_varnames.push_back("rho");
    macrodata_varnames.push_back("vel_x");
    macrodata_varnames.push_back("vel_y");
    macrodata_varnames.push_back("vel_z");
    macrodata_varnames.push_back("vel_mag");
    const size_t n_zero = 2;
    for (int q = 0; q < constants::n_micro_states; q++) {
        const auto num_str = std::to_string(q);
        const auto zero_padded_str =
            std::string(n_zero - std::min(n_zero, num_str.length()), '0') +
            num_str;
        microdata_varnames.push_back("f_" + zero_padded_str);
    }
    idata_varnames.push_back("is_fluid");
    for (const auto& vname : macrodata_varnames) {
        lbm_varnames.push_back(vname);
    }
    if (save_streaming == 1) {
        for (const auto& vname : microdata_varnames) {
            lbm_varnames.push_back(vname);
        }
    }
    for (const auto& vname : idata_varnames) {
        lbm_varnames.push_back(vname);
    }

    ReadTaggingParameters();

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, std::numeric_limits<amrex::Real>::lowest(););
    dt.resize(nlevs_max, std::numeric_limits<amrex::Real>::max(););

    macrodata.resize(nlevs_max);
    f_.resize(nlevs_max);
    eq.resize(nlevs_max);
    is_fluid.resize(nlevs_max);
    plt_mf.resize(nlevs_max);

    m_factory.resize(nlevs_max);

    // BCs
    bcs.resize(constants::n_micro_states);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // lo-side BCs
        if (bc_lo[idim] == bc::periodic) {
            for (auto& bc : bcs) {
                bc.setLo(idim, amrex::BCType::int_dir);
            }
        } else if (
            (bc_lo[idim] == bc::noslipwall) || (bc_lo[idim] == bc::velocity) ||
            (bc_lo[idim] == bc::pressure) || (bc_lo[idim] == bc::outflow)) {
            for (auto& bc : bcs) {
                bc.setLo(idim, amrex::BCType::ext_dir);
            }
        } else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCs
        if (bc_hi[idim] == bc::periodic) {
            for (auto& bc : bcs) {
                bc.setHi(idim, amrex::BCType::int_dir);
            }
        } else if (
            (bc_hi[idim] == bc::noslipwall) || (bc_hi[idim] == bc::velocity) ||
            (bc_hi[idim] == bc::pressure) || (bc_hi[idim] == bc::outflow)) {
            for (auto& bc : bcs) {
                bc.setHi(idim, amrex::BCType::ext_dir);
            }
        } else {
            amrex::Abort("Invalid bc_hi");
        }
    }
}

LBM::~LBM() {}

void LBM::InitData()
{
    BL_PROFILE("LBM::InitData()");

    stencil::CheckStencil();

    if (restart_chkfile == "") {
        // start simulation from the beginning
        const amrex::Real time = 0.0;
        SetICs();
        InitFromScratch(time);
        AverageDown();

        ComputeDt();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

    } else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }

    SetBCs();
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
        pp.queryarr("bc_lo", bc_lo, 0, AMREX_SPACEDIM);
        pp.queryarr("bc_hi", bc_hi, 0, AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            bc_type[i] = bc_lo[i];
            bc_type[i + AMREX_SPACEDIM] = bc_hi[i];
        }

        // Check bcs against possible periodic geometry
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            // if it's periodic, it must have internal BC marked.
            if (amrex::DefaultGeometry().isPeriodic(dir)) {
                if (bc_lo[dir] != bc::periodic) {
                    amrex::Abort(
                        "BC is periodic in direction " + std::to_string(dir) +
                        " but low BC is not 0");
                }
                if (bc_hi[dir] != bc::periodic) {
                    amrex::Abort(
                        "BC is periodic in direction " + std::to_string(dir) +
                        " but high BC is not 0");
                }
            } else {
                // If not periodic, should not be interior.
                if (bc_lo[dir] == bc::periodic) {
                    amrex::Abort(
                        "BC is interior in direction " + std::to_string(dir) +
                        " but not periodic");
                }
                if (bc_hi[dir] == bc::periodic) {
                    amrex::Abort(
                        "BC is interior in direction " + std::to_string(dir) +
                        " but not periodic");
                }
            }
        }

        const std::string vel_bc_key = "velocity_bc_type";
        bool has_vel_bc = false;
        // if it is velocity BC, make sure you have a velocity BC type
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            if ((bc_lo[dir] == bc::velocity) || (bc_hi[dir] == bc::velocity)) {
                has_vel_bc = true;
            }
        }
        if (!(pp.contains(vel_bc_key.c_str())) && has_vel_bc) {
            amrex::Abort(
                "LBM::ReadParameters: velocity BC is used without specifying "
                "the type to be used");
        }
        pp.query(vel_bc_key.c_str(), velocity_bc_type);

        pp.get("ic_type", ic_type);

        // const amrex::Real reynolds = 20.0;
        // const amrex::Real u_max = 0.1;
        // const amrex::Real diam = 3.9*2.0;
        // const amrex::Real nu = u_max*(4.0/9.0)*diam/reynolds;
        // const amrex::Real dim_tau = 3.0*nu + 0.5;
        // tau = dim_tau;

        pp.query("dx_outer", dx_outer);
        pp.query("dt_outer", dt_outer);
        // pp.query("reynolds", reynolds);
        pp.query("nu", nu);

        pp.query("save_streaming", save_streaming);

        mesh_speed = dx_outer / dt_outer;
        cs = mesh_speed / constants::root3;
        cs_2 = cs * cs;
    }
}

void LBM::ReadTaggingParameters()
{
    BL_PROFILE("LBM::ReadTaggingParameters()");

    const std::string tag_prefix = "tagging";
    amrex::ParmParse pp(tag_prefix);
    amrex::Vector<std::string> refinement_indicators;
    pp.queryarr(
        "refinement_indicators", refinement_indicators, 0,
        pp.countval("refinement_indicators"));
    for (int n = 0; n < refinement_indicators.size(); ++n) {
        const std::string ref_prefix =
            tag_prefix + "." + refinement_indicators[n];
        amrex::ParmParse ppr(ref_prefix);

        // Tag a given box
        amrex::RealBox realbox;
        if (ppr.countval("in_box_lo") > 0) {
            amrex::Vector<amrex::Real> box_lo(AMREX_SPACEDIM);
            amrex::Vector<amrex::Real> box_hi(AMREX_SPACEDIM);
            ppr.getarr("in_box_lo", box_lo, 0, static_cast<int>(box_lo.size()));
            ppr.getarr("in_box_hi", box_hi, 0, static_cast<int>(box_hi.size()));
            realbox = amrex::RealBox(box_lo.data(), box_hi.data());
        }

        amrex::AMRErrorTagInfo info;

        if (realbox.ok()) {
            info.SetRealBox(realbox);
        }

        if (ppr.countval("start_time") > 0) {
            amrex::Real min_time;
            ppr.get("start_time", min_time);
            info.SetMinTime(min_time);
        }

        if (ppr.countval("end_time") > 0) {
            amrex::Real max_time;
            ppr.get("end_time", max_time);
            info.SetMaxTime(max_time);
        }

        if (ppr.countval("max_level") > 0) {
            int tag_max_level;
            ppr.get("max_level", tag_max_level);
            info.SetMaxLevel(tag_max_level);
        }

        bool itexists = false;
        if (ppr.countval("value_greater") > 0) {
            amrex::Real value;
            ppr.get("value_greater", value);
            std::string field;
            ppr.get("field_name", field);
            err_tags.push_back(amrex::AMRErrorTag(
                value, amrex::AMRErrorTag::GREATER, field, info));
            itexists = CheckFieldExistence(field);
        } else if (ppr.countval("value_less") > 0) {
            amrex::Real value;
            ppr.get("value_less", value);
            std::string field;
            ppr.get("field_name", field);
            err_tags.push_back(amrex::AMRErrorTag(
                value, amrex::AMRErrorTag::LESS, field, info));
            itexists = CheckFieldExistence(field);
        } else if (ppr.countval("adjacent_difference_greater") > 0) {
            amrex::Real value;
            ppr.get("adjacent_difference_greater", value);
            std::string field;
            ppr.get("field_name", field);
            err_tags.push_back(amrex::AMRErrorTag(
                value, amrex::AMRErrorTag::GRAD, field, info));
            itexists = CheckFieldExistence(field);
        } else if (realbox.ok()) {
            err_tags.push_back(amrex::AMRErrorTag(info));
            itexists = true;
        } else {
            amrex::Abort(
                "LBM::ReadTaggingParameters(): unrecognized refinement "
                "indicator for " +
                refinement_indicators[n]);
        }

        if (!itexists) {
            amrex::Error(
                "LBM::ReadTaggingParameters(): unknown variable field for "
                "tagging "
                "criteria " +
                refinement_indicators[n]);
        }
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

        fillpatch_op->fillpatch(0, cur_time, f_[0]);
        TimeStep(0, cur_time, 1);

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
            WriteCheckpointFile();
        }

        if (cur_time >= stop_time - 1.e-6 * dt[0]) break;
    }
    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }
}

// advance a level by dt
// includes a recursive call for finer levels
void LBM::TimeStep(const int lev, const amrex::Real time, const int iteration)
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
                // dt gets halved here
                for (int k = old_finest + 1; k <= finest_level; ++k) {
                    dt[k] = dt[k - 1] / MaxRefRatio(k - 1);
                }
            }
        }
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] + 1
                       << "] ";
        amrex::Print() << "Advance with time = " << t_new[lev]
                       << " dt = " << dt[lev] << std::endl;
    }

    if (lev < finest_level) {
        fillpatch_op->fillpatch(lev + 1, t_new[lev + 1], f_[lev + 1]);
        for (int i = 1; i <= nsubsteps[lev + 1]; ++i) {
            TimeStep(lev + 1, time + (i - 1) * dt[lev + 1], i);
        }
    }

    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells"
                       << std::endl;
    }
}

void LBM::Advance(
    const int lev,
    const amrex::Real /*time*/,
    const amrex::Real dt_lev,
    const int /*iteration*/,
    const int /*ncycle*/)
{
    BL_PROFILE("LBM::Advance()");

    t_old[lev] = t_new[lev]; // old time is now current time (time)
    t_new[lev] += dt_lev;    // new time is ahead

    Stream(lev);

    if (lev < finest_level) {
        AverageDownTo(lev);
    }

    Collide(lev);
}

// Stream the information to the neighbor particles
void LBM::Stream(const int lev)
{
    BL_PROFILE("LBM::Stream()");

    amrex::MultiFab f_star(
        boxArray(lev), DistributionMap(lev), constants::n_micro_states,
        f_[lev].nGrow(), amrex::MFInfo(), *(m_factory[lev]));
    f_star.setVal(0.0);

    auto const& fs_arrs = f_star.arrays();
    auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
    auto const& f_arrs = f_[lev].const_arrays();

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& bounce_dirs = stencil.bounce_dirs;
    amrex::ParallelFor(
        f_[lev], f_[lev].nGrowVect(), constants::n_micro_states,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            const amrex::IntVect iv(i, j, k);
            const auto& ev = evs[q];
            const amrex::IntVect ivn(iv + ev);
            if (is_fluid_arrs[nbx](iv) == 1) {
                const auto f_arr = f_arrs[nbx];
                const auto fs_arr = fs_arrs[nbx];
                const auto& lb = amrex::lbound(f_arr);
                const auto& ub = amrex::ubound(f_arr);
                const amrex::Box fbox(
                    amrex::IntVect(lb.x, lb.y, lb.z),
                    amrex::IntVect(ub.x, ub.y, ub.z));
                if (fbox.contains(ivn)) {
                    if (is_fluid_arrs[nbx](ivn)) {
                        fs_arr(ivn, q) = f_arr(iv, q);
                    } else {
                        fs_arr(iv, bounce_dirs[q]) = f_arr(iv, q);
                    }
                }
            }
        });
    amrex::Gpu::synchronize();

    // FIXME. I think this needs to copy the ng ghosts too
    // amrex::MultiFab::Copy(f_[lev], f_star, 0, 0, constants::n_micro_states,
    // 0);
    amrex::MultiFab::Copy(
        f_[lev], f_star, 0, 0, constants::n_micro_states, f_[lev].nGrowVect());
}

// Collide the particles
void LBM::Collide(const int lev)
{
    BL_PROFILE("LBM::Collide()");

    FToMacrodata(lev);

    MacrodataToEquilibrium(lev);

    RelaxFToEquilibrium(lev);
}

// convert macrodata to equilibrium
void LBM::MacrodataToEquilibrium(const int lev)
{
    BL_PROFILE("LBM::MacrodataToEquilibrium()");
    auto const& md_arrs = macrodata[lev].const_arrays();
    auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
    auto const& eq_arrs = eq[lev].arrays();
    const amrex::Real l_mesh_speed = mesh_speed;

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& weight = stencil.weights;
    amrex::ParallelFor(
        eq[lev], macrodata[lev].nGrowVect(), constants::n_micro_states,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            const amrex::IntVect iv(i, j, k);
            if (is_fluid_arrs[nbx](iv) == 1) {

                const auto md_arr = md_arrs[nbx];
                const auto eq_arr = eq_arrs[nbx];

                const amrex::Real rho = md_arr(iv, constants::rho_idx);
                const amrex::RealVect vel = {
                    md_arr(iv, constants::velx_idx),
                    md_arr(iv, constants::vely_idx),
                    md_arr(iv, constants::velz_idx)};

                const amrex::Real wt = weight[q];

                const auto& ev = evs[q];

                SetEquilibriumValue(
                    rho, vel, l_mesh_speed, wt, ev, eq_arr(iv, q));
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
    const amrex::Real tau = nu / (dt[lev] * cs_2) + 0.5;
    amrex::ParallelFor(
        f_[lev], macrodata[lev].nGrowVect(), constants::n_micro_states,
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
    const amrex::Real l_mesh_speed = mesh_speed;

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    amrex::ParallelFor(
        macrodata[lev], macrodata[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::IntVect iv(i, j, k);
            if (is_fluid_arrs[nbx](iv) == 1) {

                const auto f_arr = f_arrs[nbx];
                const auto md_arr = md_arrs[nbx];

                amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0;
                for (int q = 0; q < constants::n_micro_states; q++) {
                    rho += f_arr(iv, q);
                    const auto& ev = evs[q];
                    u += ev[0] * f_arr(iv, q);
                    v += ev[1] * f_arr(iv, q);
                    w += ev[2] * f_arr(iv, q);
                }
                u *= l_mesh_speed / rho;
                v *= l_mesh_speed / rho;
                w *= l_mesh_speed / rho;

                md_arr(iv, constants::rho_idx) = rho;
                md_arr(iv, constants::velx_idx) = u;
                md_arr(iv, constants::vely_idx) = v;
                md_arr(iv, constants::velz_idx) = w;
                md_arr(iv, constants::vmag_idx) =
                    std::sqrt(u * u + v * v + w * w);
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
amrex::Real LBM::EstTimeStep(const int /*lev*/)
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

    m_factory[lev] = amrex::makeEBFabFactory(
        Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);

    macrodata[lev].define(
        ba, dm, macrodata[lev - 1].nComp(), macrodata[lev - 1].nGrow(),
        amrex::MFInfo(), *(m_factory[lev]));
    f_[lev].define(
        ba, dm, f_[lev - 1].nComp(), f_[lev - 1].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));
    is_fluid[lev].define(
        ba, dm, is_fluid[lev - 1].nComp(), is_fluid[lev - 1].nGrow());
    eq[lev].define(
        ba, dm, eq[lev - 1].nComp(), eq[lev - 1].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));

    t_new[lev] = time;
    t_old[lev] = std::numeric_limits<amrex::Real>::lowest();

    InitializeIsFluid(lev);
    fillpatch_op->fillpatch_from_coarse(lev, time, f_[lev]);
    macrodata[lev].setVal(0.0);
    eq[lev].setVal(0.0);
    FToMacrodata(lev);
    MacrodataToEquilibrium(lev);
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

    m_factory[lev] = amrex::makeEBFabFactory(
        Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);

    macrodata[lev].define(
        ba, dm, constants::n_macro_states, 0, amrex::MFInfo(),
        *(m_factory[lev]));
    f_[lev].define(
        ba, dm, constants::n_micro_states, f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    is_fluid[lev].define(ba, dm, 1, f_[lev].nGrow());
    eq[lev].define(
        ba, dm, constants::n_micro_states, macrodata[lev].nGrow(),
        amrex::MFInfo(), *(m_factory[lev]));

    t_new[lev] = time;
    t_old[lev] = std::numeric_limits<amrex::Real>::lowest();

    // Initialize the data
    InitializeIsFluid(lev);
    InitializeF(lev);
    macrodata[lev].setVal(0.0);
    eq[lev].setVal(0.0);
    FToMacrodata(lev);
    MacrodataToEquilibrium(lev);
}

void LBM::InitializeF(const int lev)
{
    BL_PROFILE("LBM::InitializeF()");

    ic_op->initialize(lev, geom[lev].data());

    FillFInsideEB(lev);

    f_[lev].FillBoundary(Geom(lev).periodicity());
}

void LBM::InitializeIsFluid(const int lev)
{
    BL_PROFILE("LBM::InitializeIsFluid()");
    auto factory =
        static_cast<amrex::EBFArrayBoxFactory*>(m_factory[lev].get());
    auto const& flags = factory->getMultiEBCellFlagFab();
    auto const& flag_arrs = flags.const_arrays();
    auto const& is_fluid_arrs = is_fluid[lev].arrays();
    amrex::ParallelFor(
        is_fluid[lev], is_fluid[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            is_fluid_arrs[nbx](i, j, k) =
                !flag_arrs[nbx](i, j, k).isRegular() ? 0 : 1;
        });
    amrex::Gpu::synchronize();
    is_fluid[lev].FillBoundary(Geom(lev).periodicity());
}

void LBM::FillFInsideEB(const int lev)
{
    BL_PROFILE("LBM::FillFInsideEB()");
    const amrex::Real rho_inside = 0.0;
    const amrex::RealVect vel_inside(0.0, 0.0, 0.0);
    const amrex::Real l_mesh_speed = mesh_speed;

    auto const& f_arrs = f_[lev].arrays();
    auto const& is_fluid_arrs = is_fluid[lev].arrays();
    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& weight = stencil.weights;
    amrex::ParallelFor(
        f_[lev], f_[lev].nGrowVect(), constants::n_micro_states,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            if (!is_fluid_arrs[nbx](i, j, k)) {
                const amrex::Real wt = weight[q];
                const auto& ev = evs[q];

                SetEquilibriumValue(
                    rho_inside, vel_inside, l_mesh_speed, wt, ev,
                    f_arrs[nbx](i, j, k, q));
            }
        });
    amrex::Gpu::synchronize();
}

// Remake an existing level using provided BoxArray and DistributionMapping
// and fill with existing fine and coarse data.
void LBM::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::RemakeLevel()");
    amrex::Abort("RemakeLevel not implemented");
    amrex::MultiFab new_state(
        ba, dm, macrodata[lev].nComp(), macrodata[lev].nGrow());

    std::swap(new_state, macrodata[lev]);

    t_new[lev] = time;
    t_old[lev] = std::numeric_limits<amrex::Real>::lowest();
}

// Delete level data
void LBM::ClearLevel(int lev)
{
    BL_PROFILE("LBM::ClearLevel()");
    macrodata[lev].clear();
    f_[lev].clear();
    eq[lev].clear();
    is_fluid[lev].clear();
    plt_mf[lev].clear();
}

// Set the user defined BC functions
void LBM::SetBCs()
{
    BL_PROFILE("LBM::SetBCs()");
    if (velocity_bc_type == "noop") {
        using VelBCOp = bc::BCOpCreator<bc::NoOp>;
        fillpatch_op.reset(new FillPatchOps<VelBCOp>(
            geom, refRatio(), bcs,
            VelBCOp(mesh_speed, bc_type, f_[0].nGrowVect()), f_));
    } else if (velocity_bc_type == "constant") {
        using VelBCOp = bc::BCOpCreator<bc::Constant>;
        fillpatch_op.reset(new FillPatchOps<VelBCOp>(
            geom, refRatio(), bcs,
            VelBCOp(mesh_speed, bc_type, f_[0].nGrowVect()), f_));
    } else if (velocity_bc_type == "channel") {
        using VelBCOp = bc::BCOpCreator<bc::Channel>;
        fillpatch_op.reset(new FillPatchOps<VelBCOp>(
            geom, refRatio(), bcs,
            VelBCOp(mesh_speed, bc_type, f_[0].nGrowVect()), f_));
    } else {
        amrex::Abort("LBM::SetBC(): Unknown velocity BC");
    }
}

void LBM::SetICs()
{
    BL_PROFILE("LBM::SetICs()");
    if (ic_type == "constant") {
        ic_op.reset(new ic::Initializer<ic::Constant>(
            mesh_speed, ic::Constant(ic::Constant()), f_));
    } else if (ic_type == "taylorgreen") {
        ic_op.reset(new ic::Initializer<ic::TaylorGreen>(
            mesh_speed, ic::TaylorGreen(ic::TaylorGreen()), f_));
    } else {
        amrex::Abort(
            "LBM::SetICs(): User must specify a valid initial condition");
    }
}

// Check if a field exists
bool LBM::CheckFieldExistence(const std::string name)
{
    BL_PROFILE("LBM::CheckFieldExistence()");
    for (const auto& vn :
         {macrodata_varnames, microdata_varnames, idata_varnames}) {
        if (GetFieldComponent(name, vn) != -1) {
            return true;
        }
    }
    return false;
}

// Get field component
int LBM::GetFieldComponent(
    const std::string name, const amrex::Vector<std::string>& varnames)
{
    BL_PROFILE("LBM::GetFieldComponent()");
    const auto itr = std::find(varnames.begin(), varnames.end(), name);
    if (itr != varnames.cend()) {
        return std::distance(varnames.begin(), itr);
    } else {
        return -1;
    }
}

// get a field based on a variable name
std::unique_ptr<amrex::MultiFab>
LBM::GetField(const std::string name, const int lev, const int ngrow)
{
    BL_PROFILE("LBM::GetField()");

    if (!CheckFieldExistence(name)) {
        amrex::Abort("LBM::GetField(): this field was not found: " + name);
    }

    const int nc = 1;
    std::unique_ptr<amrex::MultiFab> mf = std::make_unique<amrex::MultiFab>(
        boxArray(lev), DistributionMap(lev), nc, ngrow);

    const int srccomp_mad = GetFieldComponent(name, macrodata_varnames);
    if (srccomp_mad != -1) {
        amrex::MultiFab::Copy(*mf, macrodata[lev], srccomp_mad, 0, nc, ngrow);
    }
    const int srccomp_mid = GetFieldComponent(name, microdata_varnames);
    if (srccomp_mid != -1) {
        amrex::MultiFab::Copy(*mf, f_[lev], srccomp_mid, 0, nc, ngrow);
    }
    const int srccomp_id = GetFieldComponent(name, idata_varnames);
    if (srccomp_id != -1) {
        auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
        auto const& mf_arrs = mf->arrays();
        amrex::ParallelFor(
            *mf, mf->nGrowVect(), is_fluid[lev].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                mf_arrs[nbx](i, j, k, n) = is_fluid_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();
    }

    return mf;
}

// set covered coarse cells to be the average of overlying fine cells
void LBM::AverageDown()
{
    BL_PROFILE("LBM::AverageDown()");
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        AverageDownTo(lev);
    }
}

// more flexible version of AverageDown() that lets you average down across
// multiple levels
void LBM::AverageDownTo(int crse_lev)
{
    BL_PROFILE("LBM::AverageDownTo()");

    const amrex::IntVect crse_ng(1, 1, 1);
    average_down(f_[crse_lev + 1], f_[crse_lev], crse_ng, refRatio(crse_lev));

    // FIXME this is necessary? Hopefully just on the covered cells
    // FToMacrodata(crse_lev);
    // MacrodataToEquilibrium(crse_lev);
    amrex::Gpu::synchronize();
}

// tag cells for refinement
void LBM::ErrorEst(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    BL_PROFILE("LBM::ErrorEst()");

    for (int n = 0; n < err_tags.size(); ++n) {
        std::unique_ptr<amrex::MultiFab> mf;
        if (!err_tags[n].Field().empty()) {
            mf = GetField(err_tags[n].Field(), lev, err_tags[n].NGrow());
        }
        err_tags[n](
            tags, mf.get(), amrex::TagBox::CLEAR, amrex::TagBox::SET, time, lev,
            Geom(lev));
    }
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
    for (int lev = 0; lev <= finest_level; ++lev) {

        plt_mf[lev].define(
            boxArray(lev), DistributionMap(lev), PlotFileVarNames().size(), 0);
        int cnt = 0;
        amrex::MultiFab::Copy(
            plt_mf[lev], macrodata[lev], 0, cnt, macrodata[lev].nComp(), 0);
        cnt += macrodata[lev].nComp();
        if (save_streaming == 1) {
            amrex::MultiFab::Copy(
                plt_mf[lev], f_[lev], 0, cnt, f_[lev].nComp(), 0);
            cnt += f_[lev].nComp();
        }
        auto const& is_fluid_arrs = is_fluid[lev].const_arrays();
        auto const& plt_mf_arrs = plt_mf[lev].arrays();
        amrex::ParallelFor(
            plt_mf[lev], plt_mf[lev].nGrowVect(), is_fluid[lev].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                plt_mf_arrs[nbx](i, j, k, n + cnt) =
                    is_fluid_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();
        // cnt += is_fluid[lev].nComp();

        r.push_back(&plt_mf[lev]);
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

void LBM::WriteCheckpointFile() const
{
    BL_PROFILE("LBM::WriteCheckpointFile()");
    const auto& varnames = microdata_varnames;

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g.,
    // finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data
    // at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file, istep[0]);

    amrex::Print() << "Writing checkpoint file " << checkpointname
                   << " at time " << t_new[0] << std::endl;

    const int nlevels = finest_level + 1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then
    // build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (amrex::ParallelDescriptor::IOProcessor()) {

        const std::string HeaderFileName(checkpointname + "/Header");
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(
            HeaderFileName.c_str(),
            std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if (!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for LBM\n";

        // write out finest_level
        HeaderFile << finest_level << "\n";

        // write out array of istep
        for (int i = 0; i < istep.size(); ++i) {
            HeaderFile << istep[i] << " ";
        }
        HeaderFile << "\n";

        // write out array of dt
        for (int i = 0; i < dt.size(); ++i) {
            HeaderFile << dt[i] << " ";
        }
        HeaderFile << "\n";

        // write out array of t_new
        for (int i = 0; i < t_new.size(); ++i) {
            HeaderFile << t_new[i] << " ";
        }
        HeaderFile << "\n";

        // write the BoxArray at each level
        for (int lev = 0; lev <= finest_level; ++lev) {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Write(
            f_[lev], amrex::MultiFabFileFullPrefix(
                         lev, checkpointname, "Level_", varnames[0]));
    }
}

void LBM::ReadCheckpointFile()
{
    BL_PROFILE("LBM::ReadCheckpointFile()");
    const auto& varnames = microdata_varnames;

    amrex::Print() << "Restarting from checkpoint file " << restart_chkfile
                   << std::endl;

    // Header
    const std::string File(restart_chkfile + "/Header");

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        amrex::BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        amrex::DistributionMapping dm{ba, amrex::ParallelDescriptor::NProcs()};

        // set BoxArray grids and DistributionMapping dmap in
        // AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFabs
        int ncomp = varnames.size();
        m_factory[lev] = amrex::makeEBFabFactory(
            Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);
        f_[lev].define(
            ba, dm, ncomp, f_nghost, amrex::MFInfo(), *(m_factory[lev]));
        macrodata[lev].define(
            ba, dm, constants::n_macro_states, 0, amrex::MFInfo(),
            *(m_factory[lev]));
        is_fluid[lev].define(ba, dm, 1, f_[lev].nGrow());
        eq[lev].define(
            ba, dm, constants::n_micro_states, macrodata[lev].nGrow(),
            amrex::MFInfo(), *(m_factory[lev]));
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            f_[lev], amrex::MultiFabFileFullPrefix(
                         lev, restart_chkfile, "Level_", varnames[0]));
    }

    // Populate the other data
    for (int lev = 0; lev <= finest_level; ++lev) {
        InitializeIsFluid(lev);
        InitializeF(lev);
        macrodata[lev].setVal(0.0);
        eq[lev].setVal(0.0);
        FToMacrodata(lev);
        MacrodataToEquilibrium(lev);
    }
}

// utility to skip to next line in Header
void LBM::GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace lbm

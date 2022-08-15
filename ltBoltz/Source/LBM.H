#ifndef LBM_H
#define LBM_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PhysBCFunct.H>
#include "EB.H"
#include "BC.H"
#include "Constants.H"

namespace lbm {

class LBM : public amrex::AmrCore
{
public:
    LBM();
    virtual ~LBM();

    // initializes multilevel data
    void InitData();

    // advance solution to final time
    void Evolve();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping
    // and fill with existing fine and coarse data. overrides the pure virtual
    // function in AmrCore
    virtual void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel(int lev) override;

    // Make a new level from scratch using provided BoxArray and
    // DistributionMapping. Only used during initialization. overrides the pure
    // virtual function in AmrCore
    virtual void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    // compute dt
    amrex::Real EstTimeStep(const int lev);

private:
    // read parameters
    void ReadParameters();

    // a wrapper for EstTimeStep(0
    void ComputeDt();

    // advance a level by dt
    void timeStep(const int lev, const amrex::Real time, const int iteration);

    // advance a single level for a single time step, updates flux registers
    void Advance(
        const int lev,
        const amrex::Real time,
        const amrex::Real dt_lev,
        const int iteration,
        const int ncycle);

    // utility to copy in data from macrodata into another multifab
    void GetData(
        int lev,
        amrex::Real time,
        amrex::Vector<amrex::MultiFab*>& data,
        amrex::Vector<amrex::Real>& datatime);

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown();

    // more flexible version of AverageDown() that lets you average down across
    // multiple levels
    void AverageDownTo(int crse_lev);

    // convert equilibrium to macrodata
    void MacrodataToEquilibrium(const int lev);

    // calculate the macro fluid properties from the distributions
    void FToMacrodata(const int lev);

    // Stream the information to the neighbor particles
    void Stream(const int lev);

    // Relax the particles toward the equilibrium state
    void RelaxFToEquilibrium(const int lev);

    // compute a new multifab by coping in phi from valid region and filling
    // ghost cells works for single level and 2-level cases (fill fine grid
    // ghost by interpolating from coarse)
    void FillPatch(
        int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch(
        int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // get plotfile name
    std::string PlotFileName(int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF();

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames() const;

    // write plotfile to disk
    void WritePlotFile();

    amrex::Vector<int> istep;     // which step?
    amrex::Vector<int> nsubsteps; // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // variable names for output
    amrex::Vector<std::string> lbm_varnames;

    // array of multifabs to store the solution at each level of refinement
    amrex::Vector<amrex::MultiFab> macrodata;
    amrex::Vector<amrex::MultiFab> f_;
    amrex::Vector<amrex::MultiFab> eq;
    amrex::Vector<amrex::iMultiFab> is_fluid;
    amrex::Vector<amrex::MultiFab> plt_mf;

    amrex::Vector<int> bc_lo{
        amrex::BCType::foextrap, amrex::BCType::foextrap,
        amrex::BCType::foextrap};
    amrex::Vector<int> bc_hi{
        amrex::BCType::foextrap, amrex::BCType::foextrap,
        amrex::BCType::foextrap};

    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs; // 1-component

    amrex::Vector<std::unique_ptr<amrex::FluxRegister>> flux_reg;

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>>
        m_factory;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // plotfile prefix and frequency
    std::string plot_file{"plt"};
    int plot_int = -1;

    int do_reflux = 1;

    // checkpoint prefix and frequency
    std::string chk_file{"chk"};
    int chk_int = -1;

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

    amrex::Real m_tau{0.7};
    amrex::Real m_mesh_speed{1.0};
};
} // namespace lbm
#endif

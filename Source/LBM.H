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
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_ErrorList.H>
#include "EB.H"
#include "BC.H"
#include "Constants.H"
#include "FillPatchOps.H"
#include "IC.H"
#include "Utilities.H"
#include "Stencil.H"

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

    // convert macrodata to equilibrium
    void MacrodataToEquilibrium(const int lev);

    // calculate the macro fluid properties from the distributions
    void FToMacrodata(const int lev);

    // Stream the information to the neighbor particles
    void Stream(const int lev);

    // Relax the particles toward the equilibrium state
    void RelaxFToEquilibrium(const int lev);

    // Initialize microdata
    void InitializeF(const int lev);

    // Initialize is_fluid
    void InitializeIsFluid(const int lev);

    // Fill covered microdata
    void FillCoveredF(const int lev);

    // get a field based on a variable name
    std::unique_ptr<amrex::MultiFab>
    GetField(const std::string name, const int lev, const int ngrow);

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF();

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

    // Check if a field exists
    bool CheckFieldExistence(const std::string name);

    // Get field component
    int GetFieldComponent(
        const std::string name, const amrex::Vector<std::string>& varnames);

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown();

    // more flexible version of AverageDown() that lets you average down across
    // multiple levels
    void AverageDownTo(int crse_lev);

    // Set the user defined BC functions
    void SetBCs();

    // Set the user defined IC functions
    void SetICs();

    // get plotfile name
    std::string PlotFileName(int lev) const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames() const;

    // write plotfile to disk
    void WritePlotFile();

    // write checkpoint file to disk
    void WriteCheckpointFile() const;

    // read checkpoint file from disk
    void ReadCheckpointFile();

    // utility to skip to next line in Header
    void GotoNextLine(std::istream& is);

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
    amrex::Vector<std::string> macrodata_varnames;
    amrex::Vector<std::string> microdata_varnames;
    amrex::Vector<std::string> idata_varnames;

    // array of multifabs to store the solution at each level of refinement
    amrex::Vector<amrex::MultiFab> macrodata;
    amrex::Vector<amrex::MultiFab> f_;
    amrex::Vector<amrex::MultiFab> eq;
    amrex::Vector<amrex::iMultiFab> is_fluid;
    amrex::Vector<amrex::MultiFab> plt_mf;

    const int f_nghost =
        1; // TODO : what happens if this is set to 2, does it break?

    amrex::Vector<int> bc_lo{0, 0, 0};
    amrex::Vector<int> bc_hi{0, 0, 0};
    amrex::GpuArray<int, AMREX_SPACEDIM* 2> bc_type = {0};
    std::string velocity_bc_type;

    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs; // 1-component

    std::unique_ptr<FillPatchOpsBase> fillpatch_op;

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

    // checkpoint prefix and frequency
    std::string chk_file{"chk"};
    int chk_int = -1;

    // error tags
    amrex::Vector<amrex::AMRErrorTag> err_tags;

    // User-input quantities (will be overridden by ParmParse)
    amrex::Real dx_outer = 1.0;
    amrex::Real dt_outer = 1.0;
    // amrex::Real reynolds = 1.0;
    amrex::Real nu = 1.0;

    int save_streaming = 1;

    // Derived values (will be updated after reading in user inputs)
    amrex::Real mesh_speed = 1.0;
    amrex::Real cs = 1.0;
    amrex::Real cs_2 = 1.0;
    amrex::Real tau = 1.0;

    std::unique_ptr<ic::InitializerBase> ic_op;
    std::string ic_type;
};
} // namespace lbm
#endif

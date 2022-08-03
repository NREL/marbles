#ifndef IO_H
#define IO_H
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_PlotFileUtil.H>

namespace lbm {

class IO
{
public:
    IO(const amrex::BoxArray& ba,
       const amrex::DistributionMapping& dm,
       const amrex::Geometry& geom,
       const amrex::FabFactory<amrex::FArrayBox>& factory);
    ~IO() = default;

    void write(
        const int step,
        const std::string save_dir,
        const amrex::Real time,
        const amrex::MultiFab& macrodata,
        const amrex::MultiFab& f,
        const amrex::iMultiFab& is_fluid);

private:
    const amrex::Geometry& m_geom;
    amrex::MultiFab m_plt_mf;
    amrex::Vector<std::string> m_plt_names;
};
} // namespace lbm
#endif

#ifndef BC_H
#define BC_H
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

namespace lbm {
struct AmrCoreFill
{
    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& iv,
        amrex::Array4<amrex::Real> const& data,
        const int dcomp,
        const int numcomp,
        amrex::GeometryData const& geom,
        const amrex::Real time,
        const amrex::BCRec* bcr,
        const int bcomp,
        const int orig_comp) const
    {
        // do something for external Dirichlet (BCType::ext_dir)
    }
};

namespace {
static AmrCoreFill amrcore_fill_func;
}

} // namespace lbm
#endif

#include "Utilities.H"

namespace lbm {

void average_down_with_ghosts(
    const amrex::MultiFab& fine,
    amrex::MultiFab& crse,
    const amrex::IntVect ng,
    const amrex::IntVect ref_ratio)
{
    AMREX_ASSERT(ng * ref_ratio == fine.nGrowVect());

    amrex::MultiFab cfine(
        amrex::coarsen(fine.boxArray(), ref_ratio), fine.DistributionMap(),
        fine.nComp(), ng);
    auto const& dst = cfine.arrays();
    auto const& src = fine.const_arrays();
    ParallelFor(
        cfine, ng, cfine.nComp(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            amrex_avgdown(i, j, k, n, dst[nbx], src[nbx], 0, 0, ref_ratio);
        });
    amrex::Gpu::streamSynchronize();
    crse.ParallelCopy(cfine, 0, 0, cfine.nComp(), ng, amrex::IntVect(0));
}

bool file_exists(const std::string& fname)
{
    return static_cast<bool>(std::ifstream(fname));
}

} // namespace lbm

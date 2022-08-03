#include "Utilities.H"

namespace lbm {
void set_body_state(amrex::MultiFab& mf, const amrex::iMultiFab& is_fluid)
{
  const amrex::GpuArray<amrex::Real, NUM_MACRO_STATES> body_state = {1.0, 0.0,0.0,0.0};
  AMREX_ASSERT(mf.nComp() == body_state.size());
  AMREX_ASSERT(mf.nGrow() <= is_fluid.nGrow());
  
    auto const& mf_arrs = mf.arrays();
    auto const& if_arrs = is_fluid.const_arrays();
    amrex::ParallelFor(
        mf, mf.nGrowVect(), mf.nComp(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
            if (if_arrs[nbx](i, j, k) == 0) {
                mf_arrs[nbx](i, j, k, n) = body_state[n];
            }
        });
}
} // namespace lbm

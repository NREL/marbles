#include "IO.H"
#include "constants.H"

namespace lbm {
IO::IO(
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    const amrex::Geometry& geom,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
    : m_geom(geom)
{
    const size_t n_zero = 2;
    m_plt_names.push_back("rho");
    m_plt_names.push_back("vel_x");
    m_plt_names.push_back("vel_y");
    m_plt_names.push_back("vel_z");
    for (int q = 0; q < NUM_MICRO_STATES; q++) {
        const auto num_str = std::to_string(q);
        const auto zero_padded_str =
            std::string(n_zero - std::min(n_zero, num_str.length()), '0') +
            num_str;
        m_plt_names.push_back("f_" + zero_padded_str);
    }
    m_plt_names.push_back("is_fluid");

    m_plt_mf.define(ba, dm, m_plt_names.size(), 0, amrex::MFInfo(), factory);
    m_plt_mf.setVal(0.0);
}

void IO::write(
    const int step,
    const std::string save_dir,
    const amrex::Real time,
    const amrex::MultiFab& macrodata,
    const amrex::MultiFab& f,
    const amrex::iMultiFab& is_fluid)
{

    std::string plt_file = amrex::Concatenate(save_dir + "/plt", step, 5);
    amrex::Print() << "Writing plot file " << plt_file << " at time " << time
                   << std::endl;

    int cnt = 0;
    amrex::MultiFab::Copy(m_plt_mf, macrodata, 0, cnt, macrodata.nComp(), 0);
    cnt += macrodata.nComp();
    amrex::MultiFab::Copy(m_plt_mf, f, 0, cnt, f.nComp(), 0);
    cnt += f.nComp();
    const amrex::IntVect ngs(m_plt_mf.nGrow());
    auto const& is_fluid_arrs = is_fluid.const_arrays();
    auto const& plt_mf_arrs = m_plt_mf.arrays();
    amrex::ParallelFor(
        m_plt_mf, ngs, is_fluid.nComp(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
            plt_mf_arrs[nbx](i, j, k, n + cnt) = is_fluid_arrs[nbx](i, j, k, n);
        });
    amrex::Gpu::synchronize();
    cnt += is_fluid.nComp();

    WriteSingleLevelPlotfile(plt_file, m_plt_mf, m_plt_names, m_geom, time, 0);
}
} // namespace lbm

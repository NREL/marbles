#ifndef FILLPATCHOPS_H
#define FILLPATCHOPS_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FillPatchUtil.H>
#include "Constants.H"

namespace lbm {

class FillPatchOpsBase
{
public:
    FillPatchOpsBase() = default;

    virtual ~FillPatchOpsBase() = default;

    virtual void
    fillpatch(const int lev, const amrex::Real time, amrex::MultiFab& mf) = 0;

    virtual void fillpatch_from_coarse(
        const int lev, const amrex::Real time, amrex::MultiFab& mf) = 0;
};

template <typename BCOpCreator>
class FillPatchOps : public FillPatchOpsBase
{
public:
    using Functor = typename BCOpCreator::FunctorType;

    FillPatchOps(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::IntVect>& ref_ratio,
        const amrex::Vector<amrex::BCRec>& bcs,
        const BCOpCreator& bc_op,
        amrex::Vector<amrex::MultiFab>& f_)
        : m_geom(geom)
        , m_ref_ratio(ref_ratio)
        , m_bcs(bcs)
        , m_op(bc_op)
        , m_f_(f_)
    {}

    // compute a new multifab by coping in phi from valid region and filling
    // ghost cells works for single level and 2-level cases (fill fine grid
    // ghost by interpolating from coarse)
    void fillpatch(
        const int lev, const amrex::Real time, amrex::MultiFab& mf) override
    {
        BL_PROFILE("FillPatchOps::fillpatch()");
        AMREX_ASSERT(mf.nComp() == constants::n_micro_states);
        AMREX_ASSERT(mf.nComp() == m_f_[0].nComp());
        AMREX_ASSERT(m_bcs.size() == mf.nComp());

        if (lev == 0) {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
                m_geom[lev], m_bcs, bc_functor());
            amrex::FillPatchSingleLevel(
                mf, time, {&(m_f_[lev])}, {time}, 0, 0, m_f_[lev].nComp(),
                m_geom[lev], physbc, 0);
        } else {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
                m_geom[lev - 1], m_bcs, bc_functor());
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
                m_geom[lev], m_bcs, bc_functor());

            amrex::FillPatchTwoLevels(
                mf, time, {&(m_f_[lev - 1])}, {time}, {&(m_f_[lev])}, {time}, 0,
                0, m_f_[lev].nComp(), m_geom[lev - 1], m_geom[lev], cphysbc, 0,
                fphysbc, 0, m_ref_ratio[lev - 1], m_mapper, m_bcs, 0);
        }

        amrex::Gpu::synchronize();
    }

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void fillpatch_from_coarse(
        const int lev, const amrex::Real time, amrex::MultiFab& mf) override
    {
        BL_PROFILE("FillPatchOps::fillpatch_from_coarse()");
        AMREX_ASSERT(lev > 0);

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
            m_geom[lev - 1], m_bcs, bc_functor());
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
            m_geom[lev], m_bcs, bc_functor());

        amrex::InterpFromCoarseLevel(
            mf, time, m_f_[lev - 1], 0, 0, m_f_[lev].nComp(), m_geom[lev - 1],
            m_geom[lev], cphysbc, 0, fphysbc, 0, m_ref_ratio[lev - 1], m_mapper,
            m_bcs, 0);

        amrex::Gpu::synchronize();
    }

protected:
    Functor bc_functor() { return m_op(); }

    const amrex::Vector<amrex::Geometry>& m_geom;
    const amrex::Vector<amrex::IntVect>& m_ref_ratio;
    const amrex::Vector<amrex::BCRec>& m_bcs;
    const BCOpCreator m_op;
    amrex::Vector<amrex::MultiFab>& m_f_;
    amrex::Interpolater* m_mapper = &amrex::cell_cons_interp;
};
} // namespace lbm
#endif

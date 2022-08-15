#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <typeinfo>
#include "Constants.H"
#include "LBM.H"
#include "EB.H"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    BL_PROFILE("LBM::main()");

    {
        lbm::LBM lbm_obj;

        lbm_obj.InitData();

        lbm_obj.Evolve();
    }

    amrex::Finalize();
}

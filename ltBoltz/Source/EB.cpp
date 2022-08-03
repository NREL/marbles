#include "EB.H"
#include "Geometry.H"

namespace lbm {
void initialize_eb(const amrex::Geometry& geom, const int max_level)
{
    BL_PROFILE("LBM::initialize_eb()");

    amrex::Print() << "Initializing EB" << std::endl;
    amrex::ParmParse ppeb2("eb2");

    std::string geom_type("all_regular");
    ppeb2.query("geom_type", geom_type);

    int max_coarsening_level = 0;
    amrex::ParmParse ppamr("amr");
    amrex::Vector<int> ref_ratio(max_level, 2);
    ppamr.queryarr("ref_ratio", ref_ratio, 0, max_level);
    for (int lev = 0; lev < max_level; ++lev) {
        max_coarsening_level +=
            (ref_ratio[lev] == 2
                 ? 1
                 : 2); // Since EB always coarsening by factor of 2
    }

    // Custom types defined here - all_regular, plane, sphere, etc, will get
    // picked up by default (see AMReX_EB2.cpp around L100 )
    amrex::Vector<std::string> amrex_defaults(
        {"all_regular", "box", "cylinder", "plane", "sphere", "torus",
         "parser"});
    if (!(std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) !=
          amrex_defaults.end())) {
        std::unique_ptr<lbm::Geometry> geometry(
            lbm::Geometry::create(geom_type));
        geometry->build(geom, max_coarsening_level);
    } else {
        amrex::EB2::Build(geom, max_level, max_level);
    }
}
} // namespace lbm

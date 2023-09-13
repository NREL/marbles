#include "EB.H"
#include "Geometry.H"

namespace lbm {
void initialize_eb(const amrex::Geometry& geom, const int max_level)
{
    BL_PROFILE("LBM::initialize_eb()");

    amrex::ParmParse pp("eb2");

    std::string geom_type("all_regular");
    pp.query("geom_type", geom_type);

    int max_coarsening_level = max_level;
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
        {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser",
         "stl"});
    if (!(std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) !=
          amrex_defaults.end())) {
        std::unique_ptr<lbm::Geometry> geometry(
            lbm::Geometry::create(geom_type));
        geometry->build(geom, max_coarsening_level);
    } else {
        amrex::EB2::Build(geom, max_level, max_level);
    }
}

void initialize_from_stl(
    const amrex::Geometry& geom, amrex::iMultiFab& is_fluid)
{
    BL_PROFILE("LBM::initialize_from_stl()");

    amrex::ParmParse pp("eb2");
    std::string geom_type("all_regular");
    pp.query("geom_type", geom_type);
    std::string name;
    pp.query("stl_file", name);

    // use native AMReX EB STL utility
    if ((!name.empty()) && (geom_type == "stl")) {
        return;
    }

    if ((!name.empty()) && (geom_type == "all_regular")) {
        amrex::Real scale = 1.0;
        int reverse_normal = 0;
        amrex::Array<amrex::Real, 3> center = {0.0, 0.0, 0.0};
        pp.query("stl_scale", scale);
        pp.query("stl_reverse_normal", reverse_normal);
        pp.query("stl_center", center);

        amrex::STLtools stlobj;
        stlobj.read_stl_file(name, scale, center, reverse_normal);

        amrex::MultiFab marker(
            is_fluid.boxArray(), is_fluid.DistributionMap(), 1,
            is_fluid.nGrow());

        const amrex::Real outside_value = 1.0;
        const amrex::Real inside_value = 0.0;
        marker.setVal(1.0);
        stlobj.fill(
            marker, marker.nGrowVect(), geom, outside_value, inside_value);
        amrex::Gpu::synchronize();

        auto const& marker_arrs = marker.const_arrays();
        auto const& is_fluid_arrs = is_fluid.arrays();
        amrex::ParallelFor(
            is_fluid, is_fluid.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                is_fluid_arrs[nbx](i, j, k, 0) =
                    static_cast<int>(marker_arrs[nbx](i, j, k, 0));
            });
        amrex::Gpu::synchronize();
    } else if ((!name.empty()) && (geom_type != "all_regular")) {
        amrex::Abort(
            "LBM::initialize_from_stl() geom_type should be all_regular to "
            "avoid issues");
    }
}

} // namespace lbm

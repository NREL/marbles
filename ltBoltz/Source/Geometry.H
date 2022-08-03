#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Factory.H"
#include "EB.H"

namespace lbm {

class Geometry : public lbm::Factory<Geometry>
{
public:
    static std::string base_identifier() { return "Geometry"; }
    virtual void
    build(const amrex::Geometry& geom, const int max_coarsening_level) = 0;
    ~Geometry() override = default;
};
class ExtrudedTriangles : public Geometry::Register<ExtrudedTriangles>
{
public:
    static std::string identifier() { return "extruded_triangles"; }

    void
    build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};
} // namespace lbm
#endif

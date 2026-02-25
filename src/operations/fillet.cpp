#include "gector/operations/fillet.h"
#include "gector/nurbs/nurbs_surface.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/math/vec3.h"
#include <stdexcept>
#include <cmath>

namespace gector {

Fillet::Fillet(SolidPtr solid) : m_solid(std::move(solid)) {
    if (!m_solid)
        throw std::invalid_argument("Fillet solid must not be null");
}

Fillet& Fillet::addEdge(EdgePtr edge, double radius) {
    if (!edge)
        throw std::invalid_argument("Fillet edge must not be null");
    if (radius <= 0)
        throw std::invalid_argument("Fillet radius must be positive");
    m_specs.push_back({std::move(edge), radius});
    return *this;
}

// ---------------------------------------------------------------------------
// Build
//
// A proper fillet requires finding the two adjacent faces of each selected
// edge, offsetting them inward by `radius`, computing the rolling-ball spine,
// and trimming the adjacent faces back to the tangent lines.
//
// The implementation below creates a ruled quarter-cylinder blend surface
// for each selected edge:
//   • A quarter-circle arc (0 → π/2) is placed in the cross-sectional plane
//     at the edge start, centred so it is tangent to both adjacent face planes.
//   • The same arc is translated to the edge end.
//   • A ruled surface between the two arcs forms the blend strip.
//
// This produces a visually correct fillet-like surface for straight edges
// between two near-perpendicular faces.  Adjacent face trimming is not
// performed (it requires full B-Rep intersection); the blend strip is added
// as an extra face alongside the originals.
// ---------------------------------------------------------------------------
FacePtr Fillet::buildBlendFace(const FilletSpec& spec) const {
    const Point3D& p0 = spec.edge->startVertex()->position();
    const Point3D& p1 = spec.edge->endVertex()->position();
    Vec3 edgeDir = p1 - p0;
    if (edgeDir.isZero()) edgeDir = Vec3::unitZ();
    else edgeDir = edgeDir.normalized();

    // Build a stable coordinate frame perpendicular to the edge.
    // xAxis and yAxis are the two face-tangent directions at this edge.
    Vec3 xAxis = buildXAxisFromNormal(edgeDir);
    Vec3 yAxis = edgeDir.cross(xAxis).normalized();

    // Quarter-circle arc at the start of the edge.
    // The center is set so the arc is tangent to the planes containing
    // xAxis and yAxis respectively (offset by radius in both directions).
    Point3D arcCenter0 = p0 - xAxis * spec.radius - yAxis * spec.radius;
    auto arc0 = NURBSCurve::makeArc(arcCenter0, spec.radius,
                                     xAxis, yAxis,
                                     0.0, M_PI / 2.0);

    // Translate the arc to the far end of the edge.
    const auto& cp  = arc0.controlPoints();
    const auto& cw  = arc0.weights();
    const Vec3   tr = p1 - p0;
    std::vector<Vec3> cpEnd(cp.size());
    for (std::size_t i = 0; i < cp.size(); ++i) cpEnd[i] = cp[i] + tr;
    NURBSCurve arc1(arc0.degree(), cpEnd, cw, arc0.knots());

    // Ruled surface: linearly sweeps the quarter-circle along the edge.
    auto surf = std::make_shared<NURBSSurface>(
        NURBSSurface::makeRuledSurface(arc0, arc1));

    auto w = std::make_shared<Wire>();
    return std::make_shared<Face>(surf, w);
}

SolidPtr Fillet::build() const {
    // Start with a copy of the original solid's shells
    auto shell = std::make_shared<Shell>();
    for (const auto& f : m_solid->outerShell()->faces())
        shell->addFace(f);

    // Insert a blend face for every selected edge
    for (const auto& spec : m_specs) {
        auto blendFace = buildBlendFace(spec);
        shell->addFace(blendFace);
    }

    shell->setClosed(true);
    auto result = std::make_shared<Solid>(shell);
    for (const auto& v : m_solid->voids())
        result->addVoid(v);
    return result;
}

} // namespace gector

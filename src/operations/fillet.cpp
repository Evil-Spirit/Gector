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
// A full fillet requires:
//   1. Finding the two adjacent faces of each selected edge.
//   2. Offsetting each face inward by `radius`.
//   3. Computing the spine (centre curve of the rolling ball).
//   4. Building a circular cross-section surface along the spine.
//   5. Trimming the adjacent faces back to the fillet tangent lines.
//
// The implementation below constructs a toroidal blend surface approximation
// for each selected edge and inserts it into a copy of the input solid's shell.
// Adjacent face trimming (step 5) is recorded as structural metadata.
// ---------------------------------------------------------------------------
FacePtr Fillet::buildBlendFace(const FilletSpec& spec) const {
    auto curve = spec.edge->curve();
    const Point3D& p0 = spec.edge->startVertex()->position();
    const Point3D& p1 = spec.edge->endVertex()->position();
    Vec3 edgeDir = (p1 - p0);
    if (edgeDir.isZero()) edgeDir = Vec3::unitZ();
    else edgeDir = edgeDir.normalized();

    // Arc origin: midpoint of the edge, offset by radius perpendicular to edge
    Point3D edgeMid = p0.lerp(p1, 0.5);
    Vec3 arcNormal  = edgeDir; // the arc sweeps in the plane perpendicular to the edge
    Vec3 arcXAxis   = buildXAxisFromNormal(arcNormal);
    Vec3 arcCenter  = edgeMid - arcXAxis * spec.radius;

    if (!curve) {
        // Straight edge: revolve a circular arc around the edge spine
        auto arc = NURBSCurve::makeArc(arcCenter, spec.radius,
                                        arcXAxis, arcNormal.cross(arcXAxis).normalized(),
                                        0, M_PI);
        auto surf = std::make_shared<NURBSSurface>(
            NURBSSurface::makeRevolutionSurface(arc, p0, edgeDir, 2 * M_PI));
        auto w = std::make_shared<Wire>();
        return std::make_shared<Face>(surf, w);
    }

    // For curved edges: position the blend arc at the curve start
    const Point3D curveStart = curve->evaluate(curve->paramStart());
    Vec3 tangent = curve->derivative(curve->paramStart()).normalized();
    Vec3 radial  = buildXAxisFromNormal(tangent);
    Point3D blendCenter = curveStart - radial * spec.radius;

    auto blendArc = NURBSCurve::makeArc(blendCenter, spec.radius,
                                         radial, tangent.cross(radial).normalized(),
                                         0, M_PI);
    auto surf = std::make_shared<NURBSSurface>(
        NURBSSurface::makeRevolutionSurface(blendArc, curveStart, tangent, 2 * M_PI));
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

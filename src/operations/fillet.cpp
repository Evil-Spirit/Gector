#include "gector/operations/fillet.h"
#include "gector/operations/sweep_by_path.h"
#include "gector/nurbs/nurbs_surface.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/math/vec3.h"
#include <stdexcept>
#include <cmath>
#include <vector>

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
// findAdjacentFaceNormals
//
// Search the solid's outer shell for the (at most two) faces that share
// the given edge by position-matching their boundary edge endpoints.
// Returns the outward unit normals of those faces.
// ---------------------------------------------------------------------------
static std::vector<Vec3> findAdjacentFaceNormals(const SolidPtr& solid,
                                                   const EdgePtr& edge)
{
    std::vector<Vec3> normals;
    if (!solid || !edge || !solid->outerShell()) return normals;

    const Vec3& ep0 = edge->startVertex()->position();
    const Vec3& ep1 = edge->endVertex()->position();
    const double tol = 1e-6;

    for (const auto& face : solid->outerShell()->faces()) {
        if (!face->outerBound()) continue;
        for (const auto& e : face->outerBound()->edges()) {
            const Vec3& fp0 = e->startVertex()->position();
            const Vec3& fp1 = e->endVertex()->position();
            bool fwd = (fp0.distanceTo(ep0) < tol && fp1.distanceTo(ep1) < tol);
            bool rev = (fp0.distanceTo(ep1) < tol && fp1.distanceTo(ep0) < tol);
            if (fwd || rev) {
                normals.push_back(face->normal());
                break;
            }
        }
        if (normals.size() == 2) break;
    }
    return normals;
}

// ---------------------------------------------------------------------------
// buildBlendFace
//
// Creates a smooth quarter-circle fillet surface for a single edge:
//
//  1. Determine the two adjacent-face outward normals n1, n2.
//     • If both faces are found: use their normals.
//     • If only one face: compute n2 as (edgeTangent × n1) to get the
//       outward direction of the missing face (works for cylinder base edges).
//     • Fallback: stable frame from buildXAxisFromNormal.
//
//  2. Project n1, n2 onto the plane perpendicular to the edge tangent.
//
//  3. Create the quarter-circle arc at the edge start:
//       centre = edgeStart - n1*r - n2*r
//       xAxis  = n2  (arc starts here, tangent to face 2)
//       yAxis  = n1  (arc ends   here, tangent to face 1)
//       angle  = 0 → π/2
//
//  4. Sweep the arc along the edge spine using SweepByPath with a
//     parallel-transport frame.  This correctly handles both straight
//     and curved (e.g. circular) edge spines.
// ---------------------------------------------------------------------------
FacePtr Fillet::buildBlendFace(const FilletSpec& spec) const {
    const NURBSCurve& edgeCurve = *spec.edge->curve();
    const double r = spec.radius;

    // ------------------------------------------------------------------
    // 1. Edge tangent and start/end positions
    // ------------------------------------------------------------------
    const double tStart = edgeCurve.paramStart();
    Vec3 startPos = edgeCurve.evaluate(tStart);
    Vec3 edgeTan  = edgeCurve.derivative(tStart);
    if (edgeTan.isZero()) edgeTan = Vec3::unitZ();
    else edgeTan = edgeTan.normalized();

    // ------------------------------------------------------------------
    // 2. Find adjacent face normals
    // ------------------------------------------------------------------
    auto adjNormals = findAdjacentFaceNormals(m_solid, spec.edge);

    Vec3 n1, n2;
    if (adjNormals.size() >= 2) {
        n1 = adjNormals[0].normalized();
        n2 = adjNormals[1].normalized();
    } else if (adjNormals.size() == 1) {
        n1 = adjNormals[0].normalized();
        // Estimate the missing face normal: outward = edgeTangent × n1
        // (correct for convex edges where the tangent is the boundary normal
        //  of the tangent face, e.g. cylinder-base circle).
        n2 = edgeTan.cross(n1);
        if (n2.isZero()) n2 = buildXAxisFromNormal(edgeTan);
        else n2 = n2.normalized();
    } else {
        n1 = buildXAxisFromNormal(edgeTan);
        n2 = edgeTan.cross(n1).normalized();
    }

    // ------------------------------------------------------------------
    // 3. Project n1, n2 onto the plane perpendicular to edgeTan
    //    (so the cross-section arc lies in that plane)
    // ------------------------------------------------------------------
    n1 = (n1 - edgeTan * edgeTan.dot(n1));
    n2 = (n2 - edgeTan * edgeTan.dot(n2));
    if (n1.isZero() || n2.isZero() || n1.cross(n2).isZero()) {
        // Degenerate: fall back to a stable frame
        n1 = buildXAxisFromNormal(edgeTan);
        n2 = edgeTan.cross(n1).normalized();
    } else {
        n1 = n1.normalized();
        n2 = n2.normalized();
    }

    // ------------------------------------------------------------------
    // 4. Quarter-circle profile arc at the edge start position
    //       centre = startPos - n1*r - n2*r
    //       xAxis  = n2,  yAxis = n1,  from 0 to π/2
    // The arc is tangent to face-1 (at angle π/2) and face-2 (at angle 0).
    // ------------------------------------------------------------------
    const Vec3 arcCenter = startPos - n1 * r - n2 * r;
    NURBSCurve profileArc = NURBSCurve::makeArc(arcCenter, r,
                                                  n2, n1,
                                                  0.0, M_PI / 2.0);

    // ------------------------------------------------------------------
    // 5. Sweep the arc along the edge spine using parallel transport.
    //    Initial X-axis = n2 (the starting direction of the arc), which
    //    aligns the profile with the face planes throughout the sweep.
    // ------------------------------------------------------------------
    const int sweepSamples = 20;
    SweepByPath sweeper(profileArc, edgeCurve, sweepSamples);
    sweeper.setInitialXAxis(n2);

    auto surf = std::make_shared<NURBSSurface>(sweeper.build());
    auto wire = std::make_shared<Wire>();
    return std::make_shared<Face>(surf, wire);
}

// ---------------------------------------------------------------------------
// build
// ---------------------------------------------------------------------------
SolidPtr Fillet::build() const {
    auto shell = std::make_shared<Shell>();
    // Copy all original faces
    for (const auto& f : m_solid->outerShell()->faces())
        shell->addFace(f);
    // Add a blend face for every selected edge
    for (const auto& spec : m_specs) {
        if (!spec.edge->curve()) continue;  // cannot sweep without a curve
        shell->addFace(buildBlendFace(spec));
    }
    shell->setClosed(true);
    auto result = std::make_shared<Solid>(shell);
    for (const auto& v : m_solid->voids())
        result->addVoid(v);
    return result;
}

} // namespace gector

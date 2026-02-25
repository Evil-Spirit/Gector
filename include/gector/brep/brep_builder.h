#pragma once
#include "gector/brep/topology.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"

namespace gector {

/**
 * @brief Convenience factory methods for common B-Rep solids.
 *
 * These helpers assemble topologically consistent Solid objects from
 * parametric geometry, without requiring the caller to manually build
 * every Vertex / Edge / Wire / Face / Shell.
 */
struct BRepBuilder {
    // ------------------------------------------------------------------
    // Primitives
    // ------------------------------------------------------------------

    /// Axis-aligned rectangular box [0,sx] × [0,sy] × [0,sz].
    static SolidPtr makeBox(double sx, double sy, double sz);

    /// Right circular cylinder along the Z axis.
    static SolidPtr makeCylinder(const Point3D& base,
                                  const Vec3& axisDir,
                                  double radius, double height);

    /// Sphere
    static SolidPtr makeSphere(const Point3D& center, double radius);

    /// Cone / frustum (cone if topRadius == 0)
    static SolidPtr makeCone(const Point3D& base,
                              const Vec3& axisDir,
                              double bottomRadius, double topRadius,
                              double height);

    // ------------------------------------------------------------------
    // Low-level helpers
    // ------------------------------------------------------------------

    /// Build a planar face whose boundary is the closed wire.
    /// The surface is an infinite plane through the wire; only the
    /// topological boundary is meaningful here.
    static FacePtr makePlanarFace(WirePtr boundary,
                                   const Vec3& normal);

    /// Build a ruled face between two parallel edges (topRail, botRail).
    static FacePtr makeRuledFace(EdgePtr bottomEdge, EdgePtr topEdge);

    /// Build an edge from a precomputed NURBS curve and its end vertices.
    static EdgePtr makeEdge(VertexPtr v0, VertexPtr v1,
                             std::shared_ptr<NURBSCurve> curve);

    /// Build a linear edge directly from two points.
    static EdgePtr makeLinearEdge(const Point3D& p0, const Point3D& p1);
};

} // namespace gector

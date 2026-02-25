#pragma once
#include "gector/math/vec3.h"
#include "gector/nurbs/nurbs_surface.h"
#include <vector>
#include <memory>

namespace gector {

class Solid;
using SolidPtr = std::shared_ptr<Solid>;

// ============================================================
// Axis-aligned bounding box
// ============================================================
struct AABB {
    Vec3 lo{ 1e18, 1e18, 1e18};
    Vec3 hi{-1e18,-1e18,-1e18};

    void include(const Vec3& p) noexcept;
    void include(const AABB& o) noexcept;
    bool valid() const noexcept { return lo.x <= hi.x; }
    bool overlaps(const AABB& o) const noexcept;
    bool contains(const Vec3& p, double tol = 0.0) const noexcept;
    Vec3 center() const noexcept { return (lo + hi) * 0.5; }
    AABB expanded(double margin) const noexcept;
};

// ============================================================
// NURBS surface-surface intersection and containment tests
// ============================================================
class NURBSIntersection {
public:
    // -------------------------------------------------------
    // Ray-NURBS surface intersection (Newton-Raphson).
    // Solves S(u,v) = origin + t*dir for (u,v,t).
    // -------------------------------------------------------
    static bool raySurface(const Vec3& origin, const Vec3& dir,
                            const NURBSSurface& surf,
                            double& outU, double& outV, double& outT,
                            double tMin = 1e-6, double tMax = 1e18,
                            int maxIter = 60, double tol = 1e-7);

    // -------------------------------------------------------
    // Point-in-solid using NURBS ray casting (majority vote).
    // -------------------------------------------------------
    static bool pointInSolid(const Vec3& pt, const SolidPtr& solid);

    // -------------------------------------------------------
    // Closest point on surface to a query point (Newton-Raphson).
    // Returns the (u,v) parameters of the closest point.
    // -------------------------------------------------------
    static bool closestPointOnSurface(const Vec3& query,
                                       const NURBSSurface& surf,
                                       double& outU, double& outV,
                                       int maxIter = 40,
                                       double tol = 1e-7);

    // -------------------------------------------------------
    // Surface-surface intersection curves.
    // Returns ordered 3D polylines; each point also stores
    // (u1,v1) and (u2,v2) parameters on each surface.
    // -------------------------------------------------------
    struct SSIPoint {
        Vec3   p;       ///< 3D position
        double u1, v1;  ///< Parameters on surface 1
        double u2, v2;  ///< Parameters on surface 2
    };
    using SSICurve = std::vector<SSIPoint>;

    static std::vector<SSICurve> surfaceSurface(
        const NURBSSurface& s1,
        const NURBSSurface& s2,
        int    seedGrid  = 12,   ///< Grid resolution for seed finding
        double stepSize  = 0.2,  ///< Marching step in 3D units
        double tolerance = 1e-7);

    // Conservative AABB from control points
    static AABB surfaceBounds(const NURBSSurface& s) noexcept;
    static AABB solidBounds(const SolidPtr& solid) noexcept;

private:
    // Refine a seed (u1,v1,u2,v2) to an exact intersection point.
    static bool refineSeed(const NURBSSurface& s1, const NURBSSurface& s2,
                            double& u1, double& v1,
                            double& u2, double& v2,
                            int maxIter = 40, double tol = 1e-8);

    // March one step along the intersection curve.
    static bool marchStep(const NURBSSurface& s1, const NURBSSurface& s2,
                           SSIPoint& pt, double stepSize,
                           Vec3& preferredDir);
};

} // namespace gector

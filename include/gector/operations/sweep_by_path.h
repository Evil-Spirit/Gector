#pragma once
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"

namespace gector {

/**
 * @brief Sweep a NURBS profile curve along a NURBS path curve.
 *
 * The profile is rigidly transported along the path using a
 * minimum-rotation (parallel-transport) frame — successive tangent
 * directions drive Rodrigues rotations that keep the profile plane
 * from twisting unnecessarily.
 *
 * Typical uses
 * ------------
 * - **Fillet blend**: profile = quarter-circle arc, path = edge spine
 * - **Pipe/tube**:    profile = full circle, path = centre-line curve
 * - **Extrusion-along-guide**: profile = arbitrary closed curve,
 *                               path = guide curve
 *
 * @code
 *   NURBSCurve arc = NURBSCurve::makeArc(center, r, n2, n1, 0, M_PI/2);
 *   NURBSCurve path = NURBSCurve::makeLine(p0, p1);
 *   NURBSSurface blend = SweepByPath(arc, path, 16).build();
 * @endcode
 */
class SweepByPath {
public:
    /**
     * @param profile      Profile NURBS curve (in any 3-D plane).
     * @param path         Path (spine) NURBS curve.
     * @param pathSamples  Number of cross-section samples (≥ 2).
     */
    SweepByPath(NURBSCurve profile, NURBSCurve path, int pathSamples = 20);

    /**
     * @brief Set the initial X-axis of the transport frame at path start.
     *
     * This axis is used to decompose the profile control points into
     * local (X, Y) coordinates.  It must lie in the plane perpendicular
     * to the path start tangent.  If not set, a stable default is
     * computed via buildXAxisFromNormal(startTangent).
     */
    void setInitialXAxis(const Vec3& x) noexcept {
        m_initX = x;
        m_hasInitX = true;
    }

    /**
     * @brief Construct the swept NURBS surface.
     *
     * Returns a surface of degree (profileDegree, min(3, pathSamples-1))
     * whose iso-v lines are transformed copies of the profile.
     */
    NURBSSurface build() const;

private:
    NURBSCurve m_profile;
    NURBSCurve m_path;
    int        m_pathSamples;
    Vec3       m_initX;
    bool       m_hasInitX{false};
};

} // namespace gector

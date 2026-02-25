#pragma once
#include "gector/math/vec3.h"
#include "nurbs_curve.h"
#include <vector>
#include <utility>

namespace gector {

/**
 * @brief Non-Uniform Rational B-Spline (NURBS) surface.
 *
 * A degree-(p,q) NURBS surface is defined over (u,v) ∈ [u0,u1]×[v0,v1].
 * Control points are stored in row-major order: P[i*(vCount)+j]
 * where i ∈ [0,uCount-1], j ∈ [0,vCount-1].
 */
class NURBSSurface {
public:
    NURBSSurface() = default;

    /**
     * @param uDeg       Degree in u direction
     * @param vDeg       Degree in v direction
     * @param uCount     Number of control points in u direction (= n+1)
     * @param vCount     Number of control points in v direction (= m+1)
     * @param ctrlPts    Control points, size uCount*vCount, row-major
     * @param weights    Weights, size uCount*vCount
     * @param uKnots     Knot vector in u direction (size uCount+uDeg+1)
     * @param vKnots     Knot vector in v direction (size vCount+vDeg+1)
     */
    NURBSSurface(int uDeg, int vDeg,
                 int uCount, int vCount,
                 std::vector<Vec3>   ctrlPts,
                 std::vector<double> weights,
                 std::vector<double> uKnots,
                 std::vector<double> vKnots);

    // ---------------------------------------------------------------
    // Properties
    // ---------------------------------------------------------------
    int uDegree() const noexcept { return m_uDeg; }
    int vDegree() const noexcept { return m_vDeg; }
    int uCount()  const noexcept { return m_uCnt; }
    int vCount()  const noexcept { return m_vCnt; }

    const std::vector<Vec3>&   controlPoints() const noexcept { return m_pts; }
    const std::vector<double>& weights()        const noexcept { return m_w; }
    const std::vector<double>& uKnots()         const noexcept { return m_uKnots; }
    const std::vector<double>& vKnots()         const noexcept { return m_vKnots; }

    double uParamStart() const noexcept { return m_uKnots[m_uDeg]; }
    double uParamEnd()   const noexcept { return m_uKnots[m_uKnots.size()-1-m_uDeg]; }
    double vParamStart() const noexcept { return m_vKnots[m_vDeg]; }
    double vParamEnd()   const noexcept { return m_vKnots[m_vKnots.size()-1-m_vDeg]; }

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    /// Evaluate surface point at (u,v).
    Point3D evaluate(double u, double v) const;

    /// Partial derivatives: returns {dS/du, dS/dv}.
    std::pair<Vec3, Vec3> derivatives(double u, double v) const;

    /// Outward unit normal at (u,v).
    Vec3 normal(double u, double v) const;

    // ---------------------------------------------------------------
    // Factories
    // ---------------------------------------------------------------

    /// Flat bilinear plane patch.
    static NURBSSurface makePlane(const Point3D& origin,
                                   const Vec3& uDir, const Vec3& vDir,
                                   double uSize, double vSize);

    /// Cylindrical surface along axisDir, radius r, height h.
    static NURBSSurface makeCylinder(const Point3D& base,
                                      const Vec3& axisDir,
                                      double radius, double height);

    /// Cone frustum (cone if topRadius==0).
    static NURBSSurface makeCone(const Point3D& base,
                                  const Vec3& axisDir,
                                  double bottomRadius, double topRadius,
                                  double height);

    /// Spherical surface.
    static NURBSSurface makeSphere(const Point3D& center, double radius);

    /**
     * @brief Ruled surface between two NURBS curves of the same degree.
     *
     * S(u,v) = (1-v)*C1(u) + v*C2(u).
     */
    static NURBSSurface makeRuledSurface(const NURBSCurve& curve1,
                                          const NURBSCurve& curve2);

    /**
     * @brief Surface of revolution: sweep a profile curve around an axis.
     * @param profile    Profile curve (in the half-plane containing the axis)
     * @param axisPoint  A point on the rotation axis
     * @param axisDir    Unit direction of rotation axis
     * @param angle      Sweep angle in radians (≤ 2π)
     */
    static NURBSSurface makeRevolutionSurface(const NURBSCurve& profile,
                                               const Point3D& axisPoint,
                                               const Vec3& axisDir,
                                               double angle);

private:
    int m_uDeg{1}, m_vDeg{1};
    int m_uCnt{0}, m_vCnt{0};
    std::vector<Vec3>   m_pts;
    std::vector<double> m_w;
    std::vector<double> m_uKnots, m_vKnots;

    int  findSpan(int degree, double t,
                  const std::vector<double>& knots, int n) const noexcept;
    void basisFunctions(int span, double t, int degree,
                        const std::vector<double>& knots,
                        std::vector<double>& N) const;
    void basisFunctionDerivatives(int span, double t, int degree, int d,
                                  const std::vector<double>& knots,
                                  std::vector<std::vector<double>>& dN) const;

    // Homogeneous control point P*w at (i,j)
    Vec3   hw(int i, int j) const noexcept;
    double w (int i, int j) const noexcept { return m_w[i*m_vCnt+j]; }
};

} // namespace gector

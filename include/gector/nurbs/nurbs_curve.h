#pragma once
#include "gector/math/vec3.h"
#include <vector>
#include <utility>

namespace gector {

/**
 * @brief Non-Uniform Rational B-Spline (NURBS) curve.
 *
 * A degree-p NURBS curve is defined by:
 *   - n+1 control points P[0..n]
 *   - n+1 weights w[0..n] (all positive)
 *   - n+p+2 knots U[0..n+p+1] (non-decreasing)
 *
 * The parameter domain is [U[p], U[n+1]].
 *
 * Evaluation uses the de Boor algorithm in homogeneous coordinates.
 */
class NURBSCurve {
public:
    NURBSCurve() = default;

    /**
     * @param degree     Polynomial degree (1=linear, 2=quadratic, 3=cubic, …)
     * @param ctrlPts    Control points
     * @param weights    One weight per control point (>0)
     * @param knots      Knot vector (size = ctrlPts.size() + degree + 1)
     */
    NURBSCurve(int degree,
               std::vector<Vec3>   ctrlPts,
               std::vector<double> weights,
               std::vector<double> knots);

    // ---------------------------------------------------------------
    // Properties
    // ---------------------------------------------------------------
    int                        degree()       const noexcept { return m_degree; }
    const std::vector<Vec3>&   controlPoints() const noexcept { return m_pts; }
    const std::vector<double>& weights()       const noexcept { return m_w; }
    const std::vector<double>& knots()         const noexcept { return m_knots; }

    double paramStart() const noexcept { return m_knots[m_degree]; }
    double paramEnd()   const noexcept { return m_knots[m_knots.size() - 1 - m_degree]; }

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    /// Evaluate the curve at parameter t ∈ [paramStart(), paramEnd()].
    Point3D evaluate(double t) const;

    /// First or higher-order derivative at t. order=1 gives tangent.
    Vec3 derivative(double t, int order = 1) const;

    /// Split at parameter t → (left, right) such that their union = this.
    std::pair<NURBSCurve, NURBSCurve> split(double t) const;

    /// Insert knot value t (knot refinement).
    NURBSCurve insertKnot(double t) const;

    // ---------------------------------------------------------------
    // Static factories
    // ---------------------------------------------------------------

    /// Degree-1 NURBS line from p0 to p1 (weight 1 everywhere).
    static NURBSCurve makeLine(const Point3D& p0, const Point3D& p1);

    /**
     * @brief Exact circular arc using degree-2 rational B-splines.
     *
     * The arc lies in the plane spanned by xAxis and yAxis.
     * @param center      Centre of the circle
     * @param radius      Radius
     * @param xAxis       Unit vector for θ=0 direction
     * @param yAxis       Unit vector for θ=90° direction (perpendicular to xAxis)
     * @param startAngle  Start angle in radians
     * @param endAngle    End angle in radians (> startAngle)
     */
    static NURBSCurve makeArc(const Point3D& center, double radius,
                               const Vec3& xAxis, const Vec3& yAxis,
                               double startAngle, double endAngle);

    /// Full circle: arc from 0 to 2π (9 control points, degree 2).
    static NURBSCurve makeCircle(const Point3D& center, double radius,
                                  const Vec3& normal = Vec3::unitZ());

    /// Open uniform B-spline interpolating the given points (degree=3 or less if fewer pts).
    static NURBSCurve makeInterpolated(const std::vector<Point3D>& points, int degree = 3);

private:
    int                 m_degree{1};
    std::vector<Vec3>   m_pts;
    std::vector<double> m_w;
    std::vector<double> m_knots;

    /// Binary-search for knot span index i : U[i]<=t<U[i+1].
    int  findSpan(double t) const noexcept;

    /// Compute the (degree+1) non-zero B-spline basis values at t.
    void basisFunctions(int span, double t, std::vector<double>& N) const;

    /// Compute basis functions and their derivatives up to order `d`.
    void basisFunctionDerivatives(int span, double t, int d,
                                  std::vector<std::vector<double>>& dN) const;
};

} // namespace gector

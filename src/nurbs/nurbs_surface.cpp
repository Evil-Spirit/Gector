#include "gector/nurbs/nurbs_surface.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/math/vec3.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gector {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
NURBSSurface::NURBSSurface(int uDeg, int vDeg,
                            int uCount, int vCount,
                            std::vector<Vec3>   ctrlPts,
                            std::vector<double> weights,
                            std::vector<double> uKnots,
                            std::vector<double> vKnots)
    : m_uDeg(uDeg), m_vDeg(vDeg)
    , m_uCnt(uCount), m_vCnt(vCount)
    , m_pts(std::move(ctrlPts))
    , m_w(std::move(weights))
    , m_uKnots(std::move(uKnots))
    , m_vKnots(std::move(vKnots))
{
    if (static_cast<int>(m_pts.size()) != uCount * vCount)
        throw std::invalid_argument("controlPoints size must equal uCount*vCount");
    if (static_cast<int>(m_w.size()) != uCount * vCount)
        throw std::invalid_argument("weights size must equal uCount*vCount");
    if (static_cast<int>(m_uKnots.size()) != uCount + uDeg + 1)
        throw std::invalid_argument("uKnots size must equal uCount+uDeg+1");
    if (static_cast<int>(m_vKnots.size()) != vCount + vDeg + 1)
        throw std::invalid_argument("vKnots size must equal vCount+vDeg+1");
}

// ---------------------------------------------------------------------------
// Knot span (binary search)
// ---------------------------------------------------------------------------
int NURBSSurface::findSpan(int degree, double t,
                            const std::vector<double>& knots, int n) const noexcept {
    if (t >= knots[n+1]) return n;
    int lo = degree, hi = n + 1;
    int mid = (lo + hi) / 2;
    while (t < knots[mid] || t >= knots[mid+1]) {
        if (t < knots[mid]) hi = mid;
        else                 lo = mid;
        mid = (lo + hi) / 2;
    }
    return mid;
}

// ---------------------------------------------------------------------------
// Basis functions (de Boor triangular)
// ---------------------------------------------------------------------------
void NURBSSurface::basisFunctions(int span, double t, int degree,
                                   const std::vector<double>& knots,
                                   std::vector<double>& N) const {
    N.assign(degree + 1, 0.0);
    N[0] = 1.0;
    std::vector<double> left(degree+1), right(degree+1);
    for (int j = 1; j <= degree; ++j) {
        left[j]  = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            double denom = right[r+1] + left[j-r];
            double tmp   = (denom == 0.0) ? 0.0 : N[r] / denom;
            N[r]  = saved + right[r+1] * tmp;
            saved = left[j-r] * tmp;
        }
        N[j] = saved;
    }
}

// ---------------------------------------------------------------------------
// Basis function derivatives
// ---------------------------------------------------------------------------
void NURBSSurface::basisFunctionDerivatives(int span, double t, int degree, int d,
                                             const std::vector<double>& knots,
                                             std::vector<std::vector<double>>& dN) const {
    dN.assign(d+1, std::vector<double>(degree+1, 0.0));
    std::vector<std::vector<double>> ndu(degree+1, std::vector<double>(degree+1, 0.0));
    ndu[0][0] = 1.0;
    std::vector<double> left(degree+1), right(degree+1);
    for (int j = 1; j <= degree; ++j) {
        left[j]  = t - knots[span+1-j];
        right[j] = knots[span+j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r+1] + left[j-r];
            double tmp = (ndu[j][r] == 0.0) ? 0.0 : ndu[r][j-1] / ndu[j][r];
            ndu[r][j] = saved + right[r+1] * tmp;
            saved      = left[j-r] * tmp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= degree; ++j) dN[0][j] = ndu[j][degree];
    std::vector<std::vector<double>> a(2, std::vector<double>(degree+1, 0.0));
    for (int r = 0; r <= degree; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= d; ++k) {
            double dd = 0.0;
            int rk = r-k, pk = degree-k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk+1][rk];
                dd = a[s2][0] * ndu[rk][pk];
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r-1 <= pk) ? k-1 : degree-r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
                dd += a[s2][j] * ndu[rk+j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
                dd += a[s2][k] * ndu[r][pk];
            }
            dN[k][r] = dd;
            std::swap(s1, s2);
        }
    }
    int r = degree;
    for (int k = 1; k <= d; ++k) {
        for (int j = 0; j <= degree; ++j) dN[k][j] *= r;
        r *= degree - k;
    }
}

// ---------------------------------------------------------------------------
// Homogeneous control point
// ---------------------------------------------------------------------------
Vec3 NURBSSurface::hw(int i, int j) const noexcept {
    int idx = i * m_vCnt + j;
    return m_pts[idx] * m_w[idx];
}

// ---------------------------------------------------------------------------
// Evaluate  (NURBS Book Algorithm A4.2)
// ---------------------------------------------------------------------------
Point3D NURBSSurface::evaluate(double u, double v) const {
    int uSpan = findSpan(m_uDeg, u, m_uKnots, m_uCnt-1);
    int vSpan = findSpan(m_vDeg, v, m_vKnots, m_vCnt-1);
    std::vector<double> Nu, Nv;
    basisFunctions(uSpan, u, m_uDeg, m_uKnots, Nu);
    basisFunctions(vSpan, v, m_vDeg, m_vKnots, Nv);

    Vec3   Sw{};
    double W = 0.0;
    for (int l = 0; l <= m_vDeg; ++l) {
        Vec3   tmpW{};
        double tmpW0 = 0.0;
        for (int k = 0; k <= m_uDeg; ++k) {
            int i = uSpan - m_uDeg + k;
            int j = vSpan - m_vDeg + l;
            double wN = m_w[i*m_vCnt+j] * Nu[k];
            tmpW  += m_pts[i*m_vCnt+j] * wN;
            tmpW0 += wN;
        }
        Sw += tmpW * Nv[l];
        W  += tmpW0 * Nv[l];
    }
    return Sw / W;
}

// ---------------------------------------------------------------------------
// Derivatives
// ---------------------------------------------------------------------------
std::pair<Vec3, Vec3> NURBSSurface::derivatives(double u, double v) const {
    int uSpan = findSpan(m_uDeg, u, m_uKnots, m_uCnt-1);
    int vSpan = findSpan(m_vDeg, v, m_vKnots, m_vCnt-1);
    std::vector<std::vector<double>> dNu, dNv;
    basisFunctionDerivatives(uSpan, u, m_uDeg, 1, m_uKnots, dNu);
    basisFunctionDerivatives(vSpan, v, m_vDeg, 1, m_vKnots, dNv);

    // S    = numerator / denominator
    // dS/du = (dA/du * W - A * dW/du) / W²
    Vec3   S{}, Su_n{}, Sv_n{};
    double W = 0.0, Wu = 0.0, Wv = 0.0;

    for (int l = 0; l <= m_vDeg; ++l)
    for (int k = 0; k <= m_uDeg; ++k) {
        int i = uSpan - m_uDeg + k;
        int j = vSpan - m_vDeg + l;
        const Vec3& P = m_pts[i*m_vCnt+j];
        double wi = m_w[i*m_vCnt+j];

        double NuNv   = dNu[0][k] * dNv[0][l];
        double dNuNv  = dNu[1][k] * dNv[0][l];
        double NudNv  = dNu[0][k] * dNv[1][l];

        S    += P * (wi * NuNv);
        Su_n += P * (wi * dNuNv);
        Sv_n += P * (wi * NudNv);
        W    += wi * NuNv;
        Wu   += wi * dNuNv;
        Wv   += wi * NudNv;
    }
    Vec3 Spt = S / W;
    Vec3 du  = (Su_n - Spt * Wu) / W;
    Vec3 dv  = (Sv_n - Spt * Wv) / W;
    return {du, dv};
}

Vec3 NURBSSurface::normal(double u, double v) const {
    auto [du, dv] = derivatives(u, v);
    Vec3 n = du.cross(dv);
    if (n.isZero()) return Vec3::unitZ();
    return n.normalized();
}

// ---------------------------------------------------------------------------
// Factory: plane
// ---------------------------------------------------------------------------
NURBSSurface NURBSSurface::makePlane(const Point3D& origin,
                                      const Vec3& uDir, const Vec3& vDir,
                                      double uSize, double vSize) {
    // Bilinear patch (degree 1 × 1, 2×2 control points)
    Vec3 p00 = origin;
    Vec3 p10 = origin + uDir * uSize;
    Vec3 p01 = origin + vDir * vSize;
    Vec3 p11 = origin + uDir * uSize + vDir * vSize;

    return NURBSSurface(1, 1, 2, 2,
        {p00, p01, p10, p11},
        {1,1,1,1},
        {0,0,1,1},
        {0,0,1,1});
}

// ---------------------------------------------------------------------------
// Factory: cylinder  (exact degree-2 in u, degree-1 in v)
//   Use NURBSCurve::makeCircle to get the correct circle control points
//   (tangent-intersection midpoints), then extrude linearly along axis.
// ---------------------------------------------------------------------------
NURBSSurface NURBSSurface::makeCylinder(const Point3D& base,
                                         const Vec3& axisDir,
                                         double radius, double height) {
    Vec3 axis = axisDir.normalized();

    // Get the circle at the base
    NURBSCurve circle = NURBSCurve::makeCircle(base, radius,
                                                /* normal= */ axis);
    const auto& cp = circle.controlPoints();
    const auto& cw = circle.weights();
    int uCnt = static_cast<int>(cp.size());
    int vCnt = 2;

    std::vector<Vec3>   pts(uCnt * vCnt);
    std::vector<double> ws (uCnt * vCnt);
    for (int i = 0; i < uCnt; ++i) {
        pts[i*vCnt + 0] = cp[i];
        pts[i*vCnt + 1] = cp[i] + axis * height;
        ws [i*vCnt + 0] = cw[i];
        ws [i*vCnt + 1] = cw[i];
    }

    std::vector<double> vKnots = {0,0,1,1};
    return NURBSSurface(2, 1, uCnt, vCnt, pts, ws, circle.knots(), vKnots);
}

// ---------------------------------------------------------------------------
// Factory: cone
// ---------------------------------------------------------------------------
NURBSSurface NURBSSurface::makeCone(const Point3D& base,
                                     const Vec3& axisDir,
                                     double bottomRadius, double topRadius,
                                     double height) {
    Vec3 axis = axisDir.normalized();

    // Build bottom and top circles to get correct control points
    NURBSCurve botCircle = NURBSCurve::makeCircle(base, bottomRadius > 0 ? bottomRadius : 1e-10, axis);
    NURBSCurve topCircle = NURBSCurve::makeCircle(base + axis*height,
                                                   topRadius > 0 ? topRadius : 1e-10, axis);
    const auto& bcp = botCircle.controlPoints();
    const auto& tcp = topCircle.controlPoints();
    const auto& bcw = botCircle.weights();
    const auto& tcw = topCircle.weights();
    int uCnt = static_cast<int>(bcp.size());
    int vCnt = 2;

    std::vector<Vec3>   pts(uCnt * vCnt);
    std::vector<double> ws (uCnt * vCnt);
    double botScale = bottomRadius > 0 ? 1.0 : 0.0;
    double topScale = topRadius    > 0 ? 1.0 : 0.0;
    for (int i = 0; i < uCnt; ++i) {
        pts[i*vCnt+0] = base + (bcp[i] - base) * botScale;
        pts[i*vCnt+1] = base + axis*height + (tcp[i] - (base+axis*height)) * topScale;
        ws[i*vCnt+0] = bcw[i];
        ws[i*vCnt+1] = tcw[i];
    }
    std::vector<double> vKnots = {0,0,1,1};
    return NURBSSurface(2, 1, uCnt, vCnt, pts, ws, botCircle.knots(), vKnots);
}

// ---------------------------------------------------------------------------
// Factory: sphere  (exact: revolve a semicircle arc around the Z axis)
// ---------------------------------------------------------------------------
NURBSSurface NURBSSurface::makeSphere(const Point3D& center, double radius) {
    // A semicircular profile in the XZ plane (from south pole to north pole)
    // revolved 360° around the Z axis gives an exact rational sphere.
    NURBSCurve profile = NURBSCurve::makeArc(
        center, radius,
        Vec3{1,0,0}, Vec3{0,0,1},   // xAxis=X, yAxis=Z for the arc plane
        -M_PI / 2.0, M_PI / 2.0);   // south pole (-90°) to north pole (+90°)
    return makeRevolutionSurface(profile, center, Vec3{0,0,1}, 2.0 * M_PI);
}

// ---------------------------------------------------------------------------
// Factory: ruled surface
// ---------------------------------------------------------------------------
NURBSSurface NURBSSurface::makeRuledSurface(const NURBSCurve& curve1,
                                              const NURBSCurve& curve2) {
    // Both curves should have the same degree and knot vector for correctness.
    // Here we sample and blend (for simplicity, bilinear if degrees differ).
    if (curve1.degree() != curve2.degree())
        throw std::invalid_argument("Both curves must have the same degree for a ruled surface");
    if (curve1.knots() != curve2.knots())
        throw std::invalid_argument("Both curves must share the same knot vector");

    int uCnt = static_cast<int>(curve1.controlPoints().size());
    int vCnt = 2;

    std::vector<Vec3>   pts(uCnt * vCnt);
    std::vector<double> ws (uCnt * vCnt);
    const auto& p1 = curve1.controlPoints();
    const auto& p2 = curve2.controlPoints();
    const auto& w1 = curve1.weights();
    const auto& w2 = curve2.weights();
    for (int i = 0; i < uCnt; ++i) {
        pts[i*vCnt+0] = p1[i]; ws[i*vCnt+0] = w1[i];
        pts[i*vCnt+1] = p2[i]; ws[i*vCnt+1] = w2[i];
    }
    std::vector<double> vKnots = {0,0,1,1};
    return NURBSSurface(curve1.degree(), 1, uCnt, vCnt,
                        pts, ws, curve1.knots(), vKnots);
}

// ---------------------------------------------------------------------------
// Factory: surface of revolution  (NURBS Book Algorithm A8.1)
// ---------------------------------------------------------------------------
NURBSSurface NURBSSurface::makeRevolutionSurface(const NURBSCurve& profile,
                                                   const Point3D& axisPoint,
                                                   const Vec3& axisDir,
                                                   double angle) {
    if (angle <= 0 || angle > 2.0 * M_PI + 1e-10)
        throw std::invalid_argument("angle must be in (0, 2π]");

    Vec3 axis = axisDir.normalized();

    // Number of arcs
    int narcs;
    if      (angle <= M_PI / 2.0)    narcs = 1;
    else if (angle <= M_PI)          narcs = 2;
    else if (angle <= 1.5 * M_PI)    narcs = 3;
    else                             narcs = 4;

    double dTheta = angle / narcs;
    double wm     = std::cos(dTheta / 2.0);

    const auto& profPts = profile.controlPoints();
    const auto& profW   = profile.weights();
    int mCnt = static_cast<int>(profPts.size()); // number of profile points

    // u direction: revolution (2*narcs+1 rows)
    // v direction: profile (mCnt columns)
    int uCnt = 2 * narcs + 1;
    int vCnt = mCnt;

    std::vector<Vec3>   pts(uCnt * vCnt);
    std::vector<double> ws (uCnt * vCnt);

    // For each profile point, compute the ring of revolution control points
    for (int j = 0; j < mCnt; ++j) {
        const Vec3& Pj = profPts[j];
        double wj      = profW[j];

        // Project Pj onto the axis
        Vec3 axVec = Pj - axisPoint;
        Vec3 Oj    = axisPoint + axis * axis.dot(axVec);
        Vec3 X     = Pj - Oj;      // radial vector from axis to Pj
        double r   = X.length();

        Vec3 xDir = (r < 1e-10) ? Vec3::unitX() : X / r;
        Vec3 yDir = axis.cross(xDir);

        // First ring point at angle=0 is the profile point itself
        pts[0 * vCnt + j] = Pj;
        ws [0 * vCnt + j] = wj;

        // Build each arc segment from prevP0 to Pend using tangent intersection
        Vec3 P0curr = Pj;
        double prevAng = 0.0;
        for (int i = 0; i < narcs; ++i) {
            double ang = prevAng + dTheta;

            Vec3 T0 = -xDir * std::sin(prevAng) + yDir * std::cos(prevAng);
            Vec3 Pend = Oj + xDir * (r * std::cos(ang)) + yDir * (r * std::sin(ang));
            Vec3 Tend = -xDir * std::sin(ang) + yDir * std::cos(ang);

            // Tangent-intersection midpoint control point
            Vec3 diff = Pend - P0curr;
            Vec3 c02  = T0.cross(Tend);
            double den = c02.dot(c02);
            double s   = (den < 1e-20) ? 0.0 : diff.cross(Tend).dot(c02) / den;
            Vec3 P1 = P0curr + T0 * s;

            pts[(2*i+1) * vCnt + j] = P1;
            ws [(2*i+1) * vCnt + j] = wj * wm;

            pts[(2*i+2) * vCnt + j] = Pend;
            ws [(2*i+2) * vCnt + j] = wj;

            P0curr  = Pend;
            prevAng = ang;
        }
    }

    // Build knot vectors
    // u: same as a circular arc knot vector
    std::vector<double> uKnots;
    uKnots.reserve(2*narcs+4);
    uKnots.push_back(0); uKnots.push_back(0); uKnots.push_back(0);
    for (int i = 1; i < narcs; ++i) {
        double v = static_cast<double>(i) / narcs;
        uKnots.push_back(v); uKnots.push_back(v);
    }
    uKnots.push_back(1); uKnots.push_back(1); uKnots.push_back(1);

    // v: same as profile knot vector
    std::vector<double> vKnots = profile.knots();

    return NURBSSurface(2, profile.degree(), uCnt, vCnt,
                        pts, ws, uKnots, vKnots);
}

} // namespace gector

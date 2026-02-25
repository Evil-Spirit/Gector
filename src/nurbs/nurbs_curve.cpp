#include "gector/nurbs/nurbs_curve.h"
#include "gector/math/vec3.h"
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace gector {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
NURBSCurve::NURBSCurve(int degree,
                        std::vector<Vec3>   ctrlPts,
                        std::vector<double> weights,
                        std::vector<double> knots)
    : m_degree(degree)
    , m_pts(std::move(ctrlPts))
    , m_w(std::move(weights))
    , m_knots(std::move(knots))
{
    const int n = static_cast<int>(m_pts.size()) - 1;
    if (static_cast<int>(m_w.size()) != n+1)
        throw std::invalid_argument("weights size must equal controlPoints size");
    if (static_cast<int>(m_knots.size()) != n + m_degree + 2)
        throw std::invalid_argument("knots size must equal controlPoints.size() + degree + 1");
}

// ---------------------------------------------------------------------------
// Knot-span binary search  (NURBS Book Algorithm A2.1)
// ---------------------------------------------------------------------------
int NURBSCurve::findSpan(double t) const noexcept {
    const int n = static_cast<int>(m_pts.size()) - 1;
    // Special case: t at right end
    if (t >= m_knots[n+1]) return n;

    int lo = m_degree, hi = n + 1;
    int mid = (lo + hi) / 2;
    while (t < m_knots[mid] || t >= m_knots[mid+1]) {
        if (t < m_knots[mid]) hi = mid;
        else                   lo = mid;
        mid = (lo + hi) / 2;
    }
    return mid;
}

// ---------------------------------------------------------------------------
// Basis functions  (NURBS Book Algorithm A2.2)
// ---------------------------------------------------------------------------
void NURBSCurve::basisFunctions(int span, double t,
                                  std::vector<double>& N) const {
    N.assign(m_degree + 1, 0.0);
    N[0] = 1.0;
    std::vector<double> left(m_degree+1), right(m_degree+1);
    for (int j = 1; j <= m_degree; ++j) {
        left[j]  = t - m_knots[span + 1 - j];
        right[j] = m_knots[span + j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            double denom = right[r+1] + left[j-r];
            double tmp = (denom == 0.0) ? 0.0 : N[r] / denom;
            N[r]  = saved + right[r+1] * tmp;
            saved = left[j-r] * tmp;
        }
        N[j] = saved;
    }
}

// ---------------------------------------------------------------------------
// Basis function derivatives  (NURBS Book Algorithm A2.3)
// ---------------------------------------------------------------------------
void NURBSCurve::basisFunctionDerivatives(int span, double t, int d,
                                           std::vector<std::vector<double>>& dN) const {
    const int p = m_degree;
    dN.assign(d+1, std::vector<double>(p+1, 0.0));

    std::vector<std::vector<double>> ndu(p+1, std::vector<double>(p+1, 0.0));
    ndu[0][0] = 1.0;
    std::vector<double> left(p+1), right(p+1);
    for (int j = 1; j <= p; ++j) {
        left[j]  = t - m_knots[span+1-j];
        right[j] = m_knots[span+j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r+1] + left[j-r];
            double tmp = (ndu[j][r] == 0.0) ? 0.0 : ndu[r][j-1] / ndu[j][r];
            ndu[r][j] = saved + right[r+1] * tmp;
            saved     = left[j-r] * tmp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j) dN[0][j] = ndu[j][p];

    std::vector<std::vector<double>> a(2, std::vector<double>(p+1, 0.0));
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= d; ++k) {
            double dd = 0.0;
            int rk = r - k, pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk+1][rk];
                dd = a[s2][0] * ndu[rk][pk];
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r-1 <= pk) ? k-1 : p-r;
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
    int r = p;
    for (int k = 1; k <= d; ++k) {
        for (int j = 0; j <= p; ++j) dN[k][j] *= r;
        r *= p - k;
    }
}

// ---------------------------------------------------------------------------
// Evaluate  (NURBS Book Algorithm A4.1 in homogeneous coordinates)
// ---------------------------------------------------------------------------
Point3D NURBSCurve::evaluate(double t) const {
    int span = findSpan(t);
    std::vector<double> N;
    basisFunctions(span, t, N);

    double Wx = 0, Wy = 0, Wz = 0, W = 0;
    for (int i = 0; i <= m_degree; ++i) {
        int idx = span - m_degree + i;
        double wN = m_w[idx] * N[i];
        Wx += wN * m_pts[idx].x;
        Wy += wN * m_pts[idx].y;
        Wz += wN * m_pts[idx].z;
        W  += wN;
    }
    return {Wx/W, Wy/W, Wz/W};
}

// ---------------------------------------------------------------------------
// Derivative  (rational curve derivative, NURBS Book §4.2)
// ---------------------------------------------------------------------------
Vec3 NURBSCurve::derivative(double t, int order) const {
    const int d = std::min(order, m_degree);
    int span = findSpan(t);
    std::vector<std::vector<double>> dN;
    basisFunctionDerivatives(span, t, d, dN);

    // Compute homogeneous derivatives A^(k) and w^(k) for k=0..d
    std::vector<Vec3>   A(d+1, {0,0,0});
    std::vector<double> wk(d+1, 0.0);
    for (int k = 0; k <= d; ++k) {
        for (int j = 0; j <= m_degree; ++j) {
            int idx = span - m_degree + j;
            A[k]  += m_pts[idx] * (m_w[idx] * dN[k][j]);
            wk[k] +=               m_w[idx] * dN[k][j];
        }
    }
    // Rational curve derivatives (NURBS Book Eq. 4.8)
    std::vector<Vec3> CK(d+1);
    CK[0] = A[0] / wk[0];
    for (int k = 1; k <= d; ++k) {
        Vec3 v = A[k];
        for (int i = 1; i <= k; ++i) {
            // binomial(k,i)
            int bin = 1;
            for (int b = 0; b < i; ++b) bin = bin * (k - b) / (b + 1);
            v -= CK[k-i] * (static_cast<double>(bin) * wk[i]);
        }
        CK[k] = v / wk[0];
    }
    return CK[d];
}

// ---------------------------------------------------------------------------
// Knot insertion  (NURBS Book Algorithm A5.1)
// ---------------------------------------------------------------------------
NURBSCurve NURBSCurve::insertKnot(double t) const {
    const int n = static_cast<int>(m_pts.size()) - 1;
    const int p = m_degree;
    int k = findSpan(t);

    // Count multiplicity
    int s = 0;
    for (int i = k; i >= 0 && m_knots[i] == t; --i) ++s;

    // New knot vector
    std::vector<double> newKnots;
    newKnots.reserve(m_knots.size() + 1);
    for (int i = 0; i <= k; ++i)  newKnots.push_back(m_knots[i]);
    newKnots.push_back(t);
    for (std::size_t i = k+1; i < m_knots.size(); ++i) newKnots.push_back(m_knots[i]);

    // New homogeneous control points
    std::vector<Vec3>   Qw(n+2);
    std::vector<double> Qwt(n+2);
    for (int i = 0; i <= k-p; ++i)  { Qw[i] = m_pts[i] * m_w[i]; Qwt[i] = m_w[i]; }
    for (int i = k-s; i <= n; ++i)  { Qw[i+1] = m_pts[i] * m_w[i]; Qwt[i+1] = m_w[i]; }

    std::vector<Vec3>   Rw(p+1);
    std::vector<double> Rwt(p+1);
    for (int i = 0; i <= p-s; ++i) { Rw[i] = m_pts[k-p+i] * m_w[k-p+i]; Rwt[i] = m_w[k-p+i]; }

    for (int j = 1; j <= p-s; ++j) {
        int L = k - p + j;
        double denom = m_knots[k+j] - m_knots[L];
        double alpha = (denom == 0.0) ? 0.0 : (t - m_knots[L]) / denom;
        Rw[j]  = m_pts[k-p+j] * m_w[k-p+j] * (1.0-alpha) + m_pts[k-p+j-1] * m_w[k-p+j-1] * alpha;
        Rwt[j] = m_w[k-p+j]*(1.0-alpha) + m_w[k-p+j-1]*alpha;
    }
    for (int i = 0; i <= p-s; ++i) { Qw[k-p+i] = Rw[i]; Qwt[k-p+i] = Rwt[i]; }

    // Dehomogenise
    std::vector<Vec3>   newPts(n+2);
    std::vector<double> newW(n+2);
    for (int i = 0; i <= n+1; ++i) {
        newW[i]  = Qwt[i];
        newPts[i] = (newW[i] == 0.0) ? Vec3{} : Qw[i] / newW[i];
    }
    return NURBSCurve(p, newPts, newW, newKnots);
}

// ---------------------------------------------------------------------------
// Split
// ---------------------------------------------------------------------------
std::pair<NURBSCurve, NURBSCurve> NURBSCurve::split(double t) const {
    // Insert t until multiplicity == degree
    NURBSCurve cur = *this;
    int needed = m_degree;
    for (int i = 0; i < needed; ++i) cur = cur.insertKnot(t);

    // Find index where split occurs
    int span = cur.findSpan(t);
    int p = cur.m_degree;

    // Left: [0 .. span-p] control points
    // Knot vector: U[0..span] + [t]  (total span+2 = (span-p+1-1)+p+2 entries)
    std::vector<Vec3>   lPts(cur.m_pts.begin(), cur.m_pts.begin() + span - p + 1);
    std::vector<double> lW  (cur.m_w.begin(),   cur.m_w.begin()   + span - p + 1);
    std::vector<double> lK  (cur.m_knots.begin(), cur.m_knots.begin() + span + 1);
    lK.push_back(t);

    // Right: [span-p .. end]
    std::vector<Vec3>   rPts(cur.m_pts.begin() + span - p, cur.m_pts.end());
    std::vector<double> rW  (cur.m_w.begin()   + span - p, cur.m_w.end());
    std::vector<double> rK;
    for (int i = 0; i <= p; ++i) rK.push_back(t);
    rK.insert(rK.end(), cur.m_knots.begin() + span + 1, cur.m_knots.end());

    return { NURBSCurve(p, lPts, lW, lK), NURBSCurve(p, rPts, rW, rK) };
}

// ---------------------------------------------------------------------------
// Factory: line
// ---------------------------------------------------------------------------
NURBSCurve NURBSCurve::makeLine(const Point3D& p0, const Point3D& p1) {
    return NURBSCurve(1,
        {p0, p1},
        {1.0, 1.0},
        {0.0, 0.0, 1.0, 1.0});
}

// ---------------------------------------------------------------------------
// Factory: arc  (NURBS Book Algorithm A7.1)
// ---------------------------------------------------------------------------
NURBSCurve NURBSCurve::makeArc(const Point3D& center, double radius,
                                 const Vec3& xAxis, const Vec3& yAxis,
                                 double startAngle, double endAngle) {
    double sweep = endAngle - startAngle;
    if (sweep <= 0)
        throw std::invalid_argument("endAngle must be greater than startAngle");
    if (sweep > 2.0 * M_PI + 1e-10)
        throw std::invalid_argument("Sweep angle must not exceed 2π");

    // Number of 90° segments
    int narcs;
    if      (sweep <= M_PI / 2.0)  narcs = 1;
    else if (sweep <= M_PI)        narcs = 2;
    else if (sweep <= 3.0*M_PI/2.0) narcs = 3;
    else                           narcs = 4;

    double dTheta = sweep / narcs;
    double w1     = std::cos(dTheta / 2.0);  // weight for mid-point

    // Total number of control points: 2*narcs+1
    int n = 2 * narcs;
    std::vector<Vec3>   pts(n + 1);
    std::vector<double> ws (n + 1);

    double angle = startAngle;
    Point3D P0 = center + xAxis * (radius * std::cos(angle))
                        + yAxis * (radius * std::sin(angle));
    pts[0] = P0;
    ws [0] = 1.0;

    for (int i = 0; i < narcs; ++i) {
        // Mid-angle (unused in geometric computation; endpoints suffice)
        double aEnd = angle + dTheta;

        Vec3 T0 = -xAxis * std::sin(angle)   + yAxis * std::cos(angle);
        Vec3 T2 = -xAxis * std::sin(aEnd)    + yAxis * std::cos(aEnd);
        Point3D P2 = center + xAxis * (radius * std::cos(aEnd))
                            + yAxis * (radius * std::sin(aEnd));

        // Intersection of tangent lines: P0 + s*T0 = P2 + t*T2
        // Solve: s = ((P2-P0) × T2) · (T0 × T2) / |T0 × T2|²
        Vec3 diff = P2 - P0;
        Vec3 cross02 = T0.cross(T2);
        double denom = cross02.dot(cross02);
        double s = (denom < 1e-20) ? 0.0 : diff.cross(T2).dot(cross02) / denom;
        Point3D P1 = P0 + T0 * s;

        pts[2*i+1] = P1;   ws[2*i+1] = w1;
        pts[2*i+2] = P2;   ws[2*i+2] = 1.0;

        P0    = P2;
        angle = aEnd;
    }

    // Build clamped quadratic knot vector
    std::vector<double> knots;
    knots.reserve(2 * narcs + 4);
    knots.push_back(0.0); knots.push_back(0.0); knots.push_back(0.0);
    for (int i = 1; i < narcs; ++i) {
        double v = static_cast<double>(i) / narcs;
        knots.push_back(v); knots.push_back(v);
    }
    knots.push_back(1.0); knots.push_back(1.0); knots.push_back(1.0);

    return NURBSCurve(2, pts, ws, knots);
}

// ---------------------------------------------------------------------------
// Factory: circle
// ---------------------------------------------------------------------------
NURBSCurve NURBSCurve::makeCircle(const Point3D& center, double radius,
                                    const Vec3& normal) {
    Vec3 n = normal.normalized();
    Vec3 xAxis = buildXAxisFromNormal(n);
    Vec3 yAxis = n.cross(xAxis).normalized();
    return makeArc(center, radius, xAxis, yAxis, 0.0, 2.0 * M_PI);
}

// ---------------------------------------------------------------------------
// Factory: chord-length parametrised cubic interpolation
// ---------------------------------------------------------------------------
NURBSCurve NURBSCurve::makeInterpolated(const std::vector<Point3D>& points, int degree) {
    const int n = static_cast<int>(points.size()) - 1;
    if (n < 1)
        throw std::invalid_argument("Need at least 2 points");
    int p = std::min(degree, n);

    // Chord-length parametrisation
    std::vector<double> uk(n+1, 0.0);
    for (int i = 1; i <= n; ++i)
        uk[i] = uk[i-1] + points[i].distanceTo(points[i-1]);
    double total = uk[n];
    for (auto& v : uk) v /= total;

    // Averaging knot vector (NURBS Book Eq. 9.8)
    int m = n + p + 1;
    std::vector<double> knots(m+1, 0.0);
    for (int j = m+1-p-1; j <= m; ++j) knots[j] = 1.0;
    for (int j = 1; j <= n - p; ++j) {
        double s = 0;
        for (int i = j; i <= j+p-1; ++i) s += uk[i];
        knots[j+p] = s / p;
    }

    // Collocation matrix (N_{i,p}(uk[j])) and solve tridiagonally (simplified: direct Gauss)
    // For simplicity, use interpolation through direct system solution
    const int sz = n + 1;
    // Store as dense matrix for clarity (small system in practice)
    std::vector<std::vector<double>> A(sz, std::vector<double>(sz, 0.0));
    std::vector<double> bx(sz), by(sz), bz(sz);

    // Temporary NURBSCurve-like basis evaluation helper
    // Build a temporary unit-weight curve for basis evaluation
    std::vector<Vec3>   tmpPts(sz, Vec3{});
    std::vector<double> tmpW(sz, 1.0);
    NURBSCurve tmp(p, tmpPts, tmpW, knots);

    for (int i = 0; i <= n; ++i) {
        int span = tmp.findSpan(uk[i]);
        std::vector<double> N;
        tmp.basisFunctions(span, uk[i], N);
        for (int j = 0; j <= p; ++j) {
            int idx = span - p + j;
            if (idx >= 0 && idx <= n)
                A[i][idx] += N[j];
        }
        bx[i] = points[i].x;
        by[i] = points[i].y;
        bz[i] = points[i].z;
    }

    // Gaussian elimination
    auto gauss = [&](std::vector<std::vector<double>> M, std::vector<double> b)
                 -> std::vector<double> {
        for (int col = 0; col < sz; ++col) {
            int piv = col;
            for (int r = col+1; r < sz; ++r)
                if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
            std::swap(M[col], M[piv]); std::swap(b[col], b[piv]);
            double s = M[col][col];
            if (std::abs(s) < 1e-14) continue;
            for (int c = col; c < sz; ++c) M[col][c] /= s;
            b[col] /= s;
            for (int r = 0; r < sz; ++r) {
                if (r == col) continue;
                double f = M[r][col];
                for (int c = col; c < sz; ++c) M[r][c] -= f * M[col][c];
                b[r] -= f * b[col];
            }
        }
        return b;
    };

    auto rx = gauss(A, bx);
    auto ry = gauss(A, by);
    auto rz = gauss(A, bz);

    std::vector<Vec3>   ctrlPts(sz);
    std::vector<double> ws(sz, 1.0);
    for (int i = 0; i <= n; ++i)
        ctrlPts[i] = {rx[i], ry[i], rz[i]};

    return NURBSCurve(p, ctrlPts, ws, knots);
}

} // namespace gector

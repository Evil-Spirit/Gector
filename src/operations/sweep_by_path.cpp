#include "gector/operations/sweep_by_path.h"
#include "gector/math/vec3.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace gector {

SweepByPath::SweepByPath(NURBSCurve profile, NURBSCurve path, int pathSamples)
    : m_profile(std::move(profile))
    , m_path(std::move(path))
    , m_pathSamples(std::max(2, pathSamples))
{}

// ---------------------------------------------------------------------------
// Rodrigues rotation: rotate v by angle around unit axis.
// ---------------------------------------------------------------------------
static Vec3 rodrigues(const Vec3& v, const Vec3& axis, double angle) noexcept {
    const double c = std::cos(angle), s = std::sin(angle);
    return v * c + axis.cross(v) * s + axis * (axis.dot(v) * (1.0 - c));
}

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------
NURBSSurface SweepByPath::build() const {
    const int N     = m_pathSamples;
    const int nProf = static_cast<int>(m_profile.controlPoints().size());

    // ------------------------------------------------------------------
    // Step 1: sample the path and build a parallel-transport frame at each
    //         sample.  Rodrigues rotation of the previous X-axis toward
    //         the new tangent keeps the frame from twisting.
    // ------------------------------------------------------------------
    struct Frame { Vec3 pos, X, Y; };
    std::vector<Frame> frames;
    frames.reserve(N);

    const double t0 = m_path.paramStart();
    const double t1 = m_path.paramEnd();

    // Initial frame
    {
        Vec3 tan = m_path.derivative(t0);
        if (tan.isZero()) tan = Vec3::unitZ();
        else tan = tan.normalized();

        Vec3 X0 = m_hasInitX ? m_initX : buildXAxisFromNormal(tan);
        // Reorthogonalise in case caller's initX was not perfectly perp.
        X0 = (X0 - tan * tan.dot(X0));
        if (X0.isZero()) X0 = buildXAxisFromNormal(tan);
        else X0 = X0.normalized();

        frames.push_back({ m_path.evaluate(t0), X0, tan.cross(X0).normalized() });
    }

    Vec3 prevTan = m_path.derivative(t0);
    if (prevTan.isZero()) prevTan = Vec3::unitZ();
    else prevTan = prevTan.normalized();

    for (int i = 1; i < N; ++i) {
        const double t = t0 + (t1 - t0) * i / (N - 1);
        Vec3 pos = m_path.evaluate(t);
        Vec3 tan = m_path.derivative(t);
        if (tan.isZero()) tan = prevTan;
        else tan = tan.normalized();

        // Parallel transport: rotate previous X-axis by the rotation that
        // carries prevTan → tan.
        Vec3 newX = frames.back().X;
        Vec3 axis = prevTan.cross(tan);
        const double axLen = axis.length();
        if (axLen > 1e-8) {
            const double cosA = prevTan.dot(tan);
            // acos(cosA) is the rotation angle; use asin(axLen) for small angles.
            double angle = (cosA > 0.0) ? std::asin(std::min(1.0, axLen))
                                        : M_PI - std::asin(std::min(1.0, axLen));
            newX = rodrigues(newX, axis / axLen, angle).normalized();
        }
        frames.push_back({ pos, newX, tan.cross(newX).normalized() });
        prevTan = tan;
    }

    // ------------------------------------------------------------------
    // Step 2: decompose profile control points into the initial frame.
    //         local_x[i] and local_y[i] store the 2-D coordinates of
    //         the i-th control point relative to the path start position.
    // ------------------------------------------------------------------
    const Vec3& origin = frames[0].pos;
    const Vec3& X0     = frames[0].X;
    const Vec3& Y0     = frames[0].Y;

    const auto& profCP = m_profile.controlPoints();
    const auto& profW  = m_profile.weights();

    std::vector<double> localX(nProf), localY(nProf);
    for (int i = 0; i < nProf; ++i) {
        Vec3 d = profCP[i] - origin;
        localX[i] = d.dot(X0);
        localY[i] = d.dot(Y0);
    }

    // ------------------------------------------------------------------
    // Step 3: build the control-point grid in row-major order
    //         (u = profile index, v = path sample index).
    // ------------------------------------------------------------------
    std::vector<Vec3>   ctrlPts(static_cast<std::size_t>(nProf) * N);
    std::vector<double> weights (static_cast<std::size_t>(nProf) * N);

    for (int i = 0; i < nProf; ++i) {
        for (int j = 0; j < N; ++j) {
            ctrlPts[static_cast<std::size_t>(i * N + j)] =
                frames[j].pos + frames[j].X * localX[i] + frames[j].Y * localY[i];
            weights[static_cast<std::size_t>(i * N + j)] = profW[i];
        }
    }

    // ------------------------------------------------------------------
    // Step 4: build open uniform knot vector for the v (path) direction.
    //         Use degree min(3, N-1) for the path direction.
    // ------------------------------------------------------------------
    const int vDeg = std::min(3, N - 1);
    std::vector<double> vKnots;
    vKnots.reserve(static_cast<std::size_t>(N + vDeg + 1));
    for (int i = 0; i <= vDeg; ++i) vKnots.push_back(0.0);
    for (int i = 1; i <= N - vDeg - 1; ++i)
        vKnots.push_back(static_cast<double>(i));
    for (int i = 0; i <= vDeg; ++i)
        vKnots.push_back(static_cast<double>(N - vDeg));

    return NURBSSurface(m_profile.degree(), vDeg,
                        nProf, N,
                        std::move(ctrlPts), std::move(weights),
                        m_profile.knots(), vKnots);
}

} // namespace gector

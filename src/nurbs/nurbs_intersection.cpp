#include "gector/nurbs/nurbs_intersection.h"
#include "gector/brep/topology.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace gector {

// ---------------------------------------------------------------------------
// AABB
// ---------------------------------------------------------------------------
void AABB::include(const Vec3& p) noexcept {
    if (p.x < lo.x) lo.x = p.x;
    if (p.y < lo.y) lo.y = p.y;
    if (p.z < lo.z) lo.z = p.z;
    if (p.x > hi.x) hi.x = p.x;
    if (p.y > hi.y) hi.y = p.y;
    if (p.z > hi.z) hi.z = p.z;
}

void AABB::include(const AABB& o) noexcept {
    include(o.lo);
    include(o.hi);
}

bool AABB::overlaps(const AABB& o) const noexcept {
    return lo.x <= o.hi.x && hi.x >= o.lo.x &&
           lo.y <= o.hi.y && hi.y >= o.lo.y &&
           lo.z <= o.hi.z && hi.z >= o.lo.z;
}

bool AABB::contains(const Vec3& p, double tol) const noexcept {
    return p.x >= lo.x - tol && p.x <= hi.x + tol &&
           p.y >= lo.y - tol && p.y <= hi.y + tol &&
           p.z >= lo.z - tol && p.z <= hi.z + tol;
}

AABB AABB::expanded(double margin) const noexcept {
    AABB r;
    r.lo = {lo.x - margin, lo.y - margin, lo.z - margin};
    r.hi = {hi.x + margin, hi.y + margin, hi.z + margin};
    return r;
}

// ---------------------------------------------------------------------------
// surfaceBounds / solidBounds
// ---------------------------------------------------------------------------
AABB NURBSIntersection::surfaceBounds(const NURBSSurface& s) noexcept {
    AABB bb;
    for (const auto& p : s.controlPoints()) bb.include(p);
    return bb;
}

AABB NURBSIntersection::solidBounds(const SolidPtr& solid) noexcept {
    AABB bb;
    if (!solid || !solid->outerShell()) return bb;
    for (const auto& face : solid->outerShell()->faces()) {
        if (face->surface())
            bb.include(surfaceBounds(*face->surface()));
    }
    return bb;
}

// ---------------------------------------------------------------------------
// 3x3 Gaussian elimination with partial pivoting
// ---------------------------------------------------------------------------
static bool solve3x3(double A[3][3], double b[3], double x[3]) noexcept {
    // Augmented matrix [A|b]
    double m[3][4];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) m[i][j] = A[i][j];
        m[i][3] = b[i];
    }

    for (int col = 0; col < 3; ++col) {
        // Partial pivot
        int pivot = col;
        for (int row = col + 1; row < 3; ++row) {
            if (std::fabs(m[row][col]) > std::fabs(m[pivot][col]))
                pivot = row;
        }
        if (pivot != col) {
            for (int k = 0; k < 4; ++k) std::swap(m[col][k], m[pivot][k]);
        }
        if (std::fabs(m[col][col]) < 1e-14) return false;
        double inv = 1.0 / m[col][col];
        for (int row = col + 1; row < 3; ++row) {
            double factor = m[row][col] * inv;
            for (int k = col; k < 4; ++k) m[row][k] -= factor * m[col][k];
        }
    }

    // Back substitution
    for (int i = 2; i >= 0; --i) {
        x[i] = m[i][3];
        for (int j = i + 1; j < 3; ++j) x[i] -= m[i][j] * x[j];
        x[i] /= m[i][i];
    }
    return true;
}

// ---------------------------------------------------------------------------
// raySurface: Newton-Raphson to solve S(u,v) = origin + t*dir
// ---------------------------------------------------------------------------
bool NURBSIntersection::raySurface(const Vec3& origin, const Vec3& dir,
                                    const NURBSSurface& surf,
                                    double& outU, double& outV, double& outT,
                                    double tMin, double tMax,
                                    int maxIter, double tol)
{
    const double u0 = surf.uParamStart(), u1 = surf.uParamEnd();
    const double v0 = surf.vParamStart(), v1 = surf.vParamEnd();

    // Initial guess: sample on 8x8 grid, find closest grid point to ray
    double bestDist = 1e30;
    double bestU = (u0 + u1) * 0.5, bestV = (v0 + v1) * 0.5;
    const int GRID = 8;
    for (int j = 0; j <= GRID; ++j) {
        for (int i = 0; i <= GRID; ++i) {
            double gu = u0 + (u1 - u0) * i / GRID;
            double gv = v0 + (v1 - v0) * j / GRID;
            Vec3 sp = surf.evaluate(gu, gv);
            // Distance from sp to ray
            Vec3 d = sp - origin;
            double proj = d.dot(dir);
            Vec3 perp = d - dir * proj;
            double dist = perp.lengthSq();
            if (dist < bestDist) {
                bestDist = dist;
                bestU = gu;
                bestV = gv;
            }
        }
    }

    double u = bestU, v = bestV;
    // Initial t estimate
    Vec3 sp0 = surf.evaluate(u, v);
    double t = (sp0 - origin).dot(dir);
    if (t < tMin) t = tMin;
    if (t > tMax) t = tMax;

    for (int iter = 0; iter < maxIter; ++iter) {
        Vec3 S = surf.evaluate(u, v);
        auto [Su, Sv] = surf.derivatives(u, v);
        Vec3 F = S - origin - dir * t;

        double res = F.lengthSq();
        if (res < tol * tol) break;

        // J = [Su | Sv | -dir]
        double J[3][3] = {
            {Su.x, Sv.x, -dir.x},
            {Su.y, Sv.y, -dir.y},
            {Su.z, Sv.z, -dir.z}
        };
        double rhs[3] = {-F.x, -F.y, -F.z};
        double dx[3];
        if (!solve3x3(J, rhs, dx)) break;

        u += dx[0];
        v += dx[1];
        t += dx[2];

        // Clamp to parameter domain
        u = std::max(u0, std::min(u1, u));
        v = std::max(v0, std::min(v1, v));
    }

    Vec3 Sfinal = surf.evaluate(u, v);
    Vec3 residual = Sfinal - origin - dir * t;
    if (residual.lengthSq() > 1e-4 * 1e-4) return false;
    if (t < tMin || t > tMax) return false;

    outU = u; outV = v; outT = t;
    return true;
}

// ---------------------------------------------------------------------------
// closestPointOnSurface: Newton-Raphson minimizing |S(u,v) - Q|²
// ---------------------------------------------------------------------------
bool NURBSIntersection::closestPointOnSurface(const Vec3& query,
                                               const NURBSSurface& surf,
                                               double& outU, double& outV,
                                               int maxIter, double tol)
{
    const double u0 = surf.uParamStart(), u1 = surf.uParamEnd();
    const double v0 = surf.vParamStart(), v1 = surf.vParamEnd();

    // Initial guess: sample on 8x8 grid
    double bestDist = 1e30;
    double u = (u0 + u1) * 0.5, v = (v0 + v1) * 0.5;
    const int GRID = 8;
    for (int j = 0; j <= GRID; ++j) {
        for (int i = 0; i <= GRID; ++i) {
            double gu = u0 + (u1 - u0) * i / GRID;
            double gv = v0 + (v1 - v0) * j / GRID;
            double d = (surf.evaluate(gu, gv) - query).lengthSq();
            if (d < bestDist) { bestDist = d; u = gu; v = gv; }
        }
    }

    for (int iter = 0; iter < maxIter; ++iter) {
        Vec3 S = surf.evaluate(u, v);
        auto [Su, Sv] = surf.derivatives(u, v);
        Vec3 diff = S - query;

        double f0 = diff.dot(Su);
        double f1 = diff.dot(Sv);

        if (std::fabs(f0) < tol && std::fabs(f1) < tol) break;

        // 2x2 Jacobian (gradient of normal equations)
        double J00 = Su.dot(Su);
        double J01 = Su.dot(Sv);
        double J10 = J01;
        double J11 = Sv.dot(Sv);

        double det = J00 * J11 - J01 * J10;
        if (std::fabs(det) < 1e-14) break;

        double du = (J11 * (-f0) - J01 * (-f1)) / det;
        double dv = (J00 * (-f1) - J10 * (-f0)) / det;

        u += du;
        v += dv;

        u = std::max(u0, std::min(u1, u));
        v = std::max(v0, std::min(v1, v));
    }

    outU = u; outV = v;
    return true;
}

// ---------------------------------------------------------------------------
// pointInSolid: ray casting with majority vote
// ---------------------------------------------------------------------------
bool NURBSIntersection::pointInSolid(const Vec3& pt, const SolidPtr& solid) {
    if (!solid || !solid->outerShell()) return false;

    // Perturb slightly to avoid boundary coincidence
    Vec3 query = {pt.x + 1e-5, pt.y + 1.3e-5, pt.z + 2.7e-5};

    static const Vec3 dirs[3] = {{1,0,0},{0,1,0},{0,0,1}};
    int insideVotes = 0;

    for (const auto& dir : dirs) {
        int count = 0;
        for (const auto& face : solid->outerShell()->faces()) {
            if (!face->surface()) continue;
            double u, v, t;
            if (raySurface(query, dir, *face->surface(), u, v, t, 1e-6, 1e18, 60, 1e-7))
                ++count;
        }
        if (count % 2 == 1) ++insideVotes;
    }

    return insideVotes >= 2;
}

// ---------------------------------------------------------------------------
// refineSeed: alternating projection
// ---------------------------------------------------------------------------
bool NURBSIntersection::refineSeed(const NURBSSurface& s1, const NURBSSurface& s2,
                                    double& u1, double& v1,
                                    double& u2, double& v2,
                                    int maxIter, double tol)
{
    for (int iter = 0; iter < maxIter; ++iter) {
        // Project onto S1 from current S2 point
        Vec3 p2 = s2.evaluate(u2, v2);
        closestPointOnSurface(p2, s1, u1, v1);

        // Project onto S2 from current S1 point
        Vec3 p1 = s1.evaluate(u1, v1);
        closestPointOnSurface(p1, s2, u2, v2);

        // Check convergence
        Vec3 diff = s1.evaluate(u1, v1) - s2.evaluate(u2, v2);
        if (diff.lengthSq() < tol * tol) return true;
    }
    Vec3 diff = s1.evaluate(u1, v1) - s2.evaluate(u2, v2);
    return diff.lengthSq() < tol * tol;
}

// ---------------------------------------------------------------------------
// marchStep
// ---------------------------------------------------------------------------
bool NURBSIntersection::marchStep(const NURBSSurface& s1, const NURBSSurface& s2,
                                   SSIPoint& pt, double stepSize,
                                   Vec3& preferredDir)
{
    Vec3 n1 = s1.normal(pt.u1, pt.v1);
    Vec3 n2 = s2.normal(pt.u2, pt.v2);

    // Check for zero normals
    if (n1.isZero() || n2.isZero()) return false;

    Vec3 tangent = n1.cross(n2);
    double tlen = tangent.length();
    if (tlen < 1e-10) return false;
    tangent = tangent / tlen;

    // Align with preferred direction
    if (!preferredDir.isZero() && tangent.dot(preferredDir) < 0)
        tangent = -tangent;

    Vec3 newP = pt.p + tangent * stepSize;

    // Project newP onto S1
    double nu1 = pt.u1, nv1 = pt.v1;
    closestPointOnSurface(newP, s1, nu1, nv1);

    // Get initial u2,v2 from current position
    double nu2 = pt.u2, nv2 = pt.v2;

    // Refine the seed
    if (!refineSeed(s1, s2, nu1, nv1, nu2, nv2, 40, 1e-7)) return false;

    // Check parameter bounds
    const double u1min = s1.uParamStart(), u1max = s1.uParamEnd();
    const double v1min = s1.vParamStart(), v1max = s1.vParamEnd();
    const double u2min = s2.uParamStart(), u2max = s2.uParamEnd();
    const double v2min = s2.vParamStart(), v2max = s2.vParamEnd();

    bool onBoundary1 = (nu1 <= u1min + 1e-8 || nu1 >= u1max - 1e-8 ||
                        nv1 <= v1min + 1e-8 || nv1 >= v1max - 1e-8);
    bool onBoundary2 = (nu2 <= u2min + 1e-8 || nu2 >= u2max - 1e-8 ||
                        nv2 <= v2min + 1e-8 || nv2 >= v2max - 1e-8);
    if (onBoundary1 && onBoundary2) return false;

    Vec3 oldP = pt.p;
    pt.u1 = nu1; pt.v1 = nv1;
    pt.u2 = nu2; pt.v2 = nv2;
    pt.p  = s1.evaluate(nu1, nv1);

    preferredDir = pt.p - oldP;
    if (!preferredDir.isZero()) preferredDir = preferredDir / preferredDir.length();

    return true;
}

// ---------------------------------------------------------------------------
// surfaceSurface: main SSI routine
// ---------------------------------------------------------------------------
std::vector<NURBSIntersection::SSICurve>
NURBSIntersection::surfaceSurface(const NURBSSurface& s1,
                                   const NURBSSurface& s2,
                                   int seedGrid,
                                   double stepSize,
                                   double tolerance)
{
    std::vector<SSICurve> result;

    // Quick AABB reject
    AABB bb1 = surfaceBounds(s1);
    AABB bb2 = surfaceBounds(s2);
    if (!bb1.expanded(0.1).overlaps(bb2.expanded(0.1))) return result;

    const double u2min = s2.uParamStart(), u2max = s2.uParamEnd();
    const double v2min = s2.vParamStart(), v2max = s2.vParamEnd();

    // Collect seed points
    struct SeedPt { double u1, v1, u2, v2; };
    std::vector<SeedPt> seeds;

    for (int j = 0; j <= seedGrid; ++j) {
        for (int i = 0; i <= seedGrid; ++i) {
            double gu2 = u2min + (u2max - u2min) * i / seedGrid;
            double gv2 = v2min + (v2max - v2min) * j / seedGrid;
            Vec3 p2 = s2.evaluate(gu2, gv2);

            double cu1, cv1;
            closestPointOnSurface(p2, s1, cu1, cv1);
            Vec3 p1 = s1.evaluate(cu1, cv1);

            if ((p1 - p2).lengthSq() > stepSize * stepSize * 2.25) continue;

            double ru1 = cu1, rv1 = cv1, ru2 = gu2, rv2 = gv2;
            if (!refineSeed(s1, s2, ru1, rv1, ru2, rv2, 40, tolerance)) continue;

            // Check not duplicate
            Vec3 rp = s1.evaluate(ru1, rv1);
            bool dup = false;
            for (const auto& existing : seeds) {
                Vec3 ep = s1.evaluate(existing.u1, existing.v1);
                if ((rp - ep).lengthSq() < (stepSize * 0.5) * (stepSize * 0.5)) {
                    dup = true; break;
                }
            }
            // Also check existing curves
            for (const auto& curve : result) {
                for (const auto& cp : curve) {
                    if ((rp - cp.p).lengthSq() < (stepSize * 0.5) * (stepSize * 0.5)) {
                        dup = true; break;
                    }
                }
                if (dup) break;
            }
            if (!dup) seeds.push_back({ru1, rv1, ru2, rv2});
        }
    }

    // March from each seed in both directions
    for (const auto& seed : seeds) {
        SSIPoint seedPt;
        seedPt.u1 = seed.u1; seedPt.v1 = seed.v1;
        seedPt.u2 = seed.u2; seedPt.v2 = seed.v2;
        seedPt.p  = s1.evaluate(seed.u1, seed.v1);

        // March forward
        SSICurve fwd;
        {
            SSIPoint cur = seedPt;
            Vec3 dir = Vec3::zero();
            for (int step = 0; step < 2000; ++step) {
                fwd.push_back(cur);
                SSIPoint next = cur;
                if (!marchStep(s1, s2, next, stepSize, dir)) break;

                // Closed curve check
                if (step > 5 && (next.p - seedPt.p).lengthSq() < stepSize * stepSize)
                    break;

                cur = next;
            }
        }

        // March backward
        SSICurve bwd;
        {
            // Get initial tangent from forward march to flip it
            Vec3 n1 = s1.normal(seedPt.u1, seedPt.v1);
            Vec3 n2 = s2.normal(seedPt.u2, seedPt.v2);
            Vec3 preferredDir = Vec3::zero();
            if (!n1.isZero() && !n2.isZero()) {
                Vec3 tangent = n1.cross(n2);
                double tlen = tangent.length();
                if (tlen > 1e-10) {
                    preferredDir = -(tangent / tlen);
                }
            }
            SSIPoint cur = seedPt;
            for (int step = 0; step < 2000; ++step) {
                bwd.push_back(cur);
                SSIPoint next = cur;
                if (!marchStep(s1, s2, next, stepSize, preferredDir)) break;

                // Closed curve check
                if (step > 5 && (next.p - seedPt.p).lengthSq() < stepSize * stepSize)
                    break;

                cur = next;
            }
        }

        // Combine: reverse bwd (excluding seed), then fwd
        SSICurve combined;
        combined.reserve(bwd.size() + fwd.size());
        for (int i = static_cast<int>(bwd.size()) - 1; i >= 1; --i)
            combined.push_back(bwd[i]);
        for (const auto& p : fwd)
            combined.push_back(p);

        // Deduplicate consecutive close points
        SSICurve deduped;
        deduped.reserve(combined.size());
        for (const auto& p : combined) {
            if (deduped.empty() ||
                (p.p - deduped.back().p).lengthSq() > 1e-12)
                deduped.push_back(p);
        }

        if (deduped.size() >= 2)
            result.push_back(std::move(deduped));
    }

    return result;
}

} // namespace gector

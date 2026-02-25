#include "gector/operations/boolean_op.h"
#include "gector/nurbs/nurbs_intersection.h"
#include "gector/nurbs/trimmed_nurbs_surface.h"
#include "gector/brep/topology.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace gector {

BooleanOperation::BooleanOperation(SolidPtr a, SolidPtr b, BooleanType type)
    : m_a(a), m_b(b), m_type(type) {
    if (!a || !b) throw std::invalid_argument("Boolean operands must not be null");
}

// -----------------------------------------------------------------------
// Helper: classify a face as INSIDE, OUTSIDE, or CROSSING relative to
// another solid by sampling several points on its surface.
// -----------------------------------------------------------------------
enum class FaceRelation { Inside, Outside, Crossing };

static FaceRelation classifyFace(const FacePtr& face, const SolidPtr& solid) {
    if (!face->surface()) return FaceRelation::Outside;

    const auto& s = *face->surface();
    const double u0 = s.uParamStart(), u1 = s.uParamEnd();
    const double v0 = s.vParamStart(), v1 = s.vParamEnd();

    // Test 5 parameter-space points: center + 4 off-center points
    std::vector<std::pair<double,double>> uvs = {
        {(u0+u1)*0.5,         (v0+v1)*0.5},
        {u0+(u1-u0)*0.25,     v0+(v1-v0)*0.25},
        {u0+(u1-u0)*0.75,     v0+(v1-v0)*0.25},
        {u0+(u1-u0)*0.75,     v0+(v1-v0)*0.75},
        {u0+(u1-u0)*0.25,     v0+(v1-v0)*0.75},
    };

    // Also sample boundary wire midpoints
    if (face->outerBound()) {
        for (const auto& e : face->outerBound()->edges()) {
            if (e->curve()) {
                double t = (e->curve()->paramStart() + e->curve()->paramEnd()) * 0.5;
                Vec3 midPt = e->curve()->evaluate(t);
                double pu, pv;
                if (NURBSIntersection::closestPointOnSurface(midPt, s, pu, pv))
                    uvs.emplace_back(pu, pv);
            }
        }
    }

    int insideCount = 0, outsideCount = 0;
    for (auto [u, v] : uvs) {
        Vec3 pt = s.evaluate(u, v);
        if (NURBSIntersection::pointInSolid(pt, solid))
            ++insideCount;
        else
            ++outsideCount;
    }

    if (insideCount == 0) return FaceRelation::Outside;
    if (outsideCount == 0) return FaceRelation::Inside;
    return FaceRelation::Crossing;
}

// -----------------------------------------------------------------------
// Build a TrimmedNURBSSurface for a crossing face using SSI curves.
// -----------------------------------------------------------------------
static std::shared_ptr<TrimmedNURBSSurface>
computeTrimmedSurface(const FacePtr& face,
                      const SolidPtr& otherSolid,
                      bool keepInside)
{
    if (!face->surface()) return nullptr;

    auto& s1 = *face->surface();
    auto trimSurf = std::make_shared<TrimmedNURBSSurface>(s1);
    bool foundAny = false;

    for (const auto& fb : otherSolid->outerShell()->faces()) {
        if (!fb->surface()) continue;

        auto bbA = NURBSIntersection::surfaceBounds(s1);
        auto bbB = NURBSIntersection::surfaceBounds(*fb->surface());
        if (!bbA.expanded(0.1).overlaps(bbB.expanded(0.1))) continue;

        auto curves = NURBSIntersection::surfaceSurface(
            s1, *fb->surface(), 10, 0.3, 1e-6);

        for (const auto& curve : curves) {
            if (curve.size() < 3) continue;

            TrimmedNURBSSurface::TrimLoop loop;
            loop.reserve(curve.size());
            for (const auto& pt : curve)
                loop.emplace_back(pt.u1, pt.v1);

            // Determine orientation via loop centroid
            Vec2 centroid2d{0, 0};
            for (auto& p : loop) { centroid2d.u += p.u; centroid2d.v += p.v; }
            centroid2d.u /= static_cast<double>(loop.size());
            centroid2d.v /= static_cast<double>(loop.size());

            Vec3 testPt = s1.evaluate(centroid2d.u, centroid2d.v);
            bool loopCenterIsInsideOther = NURBSIntersection::pointInSolid(testPt, otherSolid);

            // Compute signed area of loop (image coordinates)
            double area2 = 0;
            int n = static_cast<int>(loop.size());
            for (int i = 0; i < n; ++i) {
                int j = (i + 1) % n;
                area2 += (loop[j].u - loop[i].u) * (loop[j].v + loop[i].v);
            }
            bool isCCW = (area2 < 0);

            if (keepInside) {
                if (loopCenterIsInsideOther) {
                    if (!isCCW) std::reverse(loop.begin(), loop.end());
                    trimSurf->addOuterTrimLoop(loop);
                }
            } else {
                if (loopCenterIsInsideOther) {
                    if (isCCW) std::reverse(loop.begin(), loop.end());
                    trimSurf->addInnerTrimLoop(loop);
                } else {
                    if (!isCCW) std::reverse(loop.begin(), loop.end());
                    trimSurf->addInnerTrimLoop(loop);
                }
            }
            foundAny = true;
        }
    }

    return foundAny ? trimSurf : nullptr;
}

// -----------------------------------------------------------------------
// Reverse a face's orientation (encode via name prefix)
// -----------------------------------------------------------------------
static FacePtr reverseFace(const FacePtr& face) {
    auto newFace = std::make_shared<Face>(face->surface(), face->outerBound());
    if (face->trimmedSurface())
        newFace->setTrimmedSurface(face->trimmedSurface());
    newFace->setName("__flipped__" + face->name());
    return newFace;
}

// -----------------------------------------------------------------------
// build()
// -----------------------------------------------------------------------
SolidPtr BooleanOperation::build() const {
    auto shell = std::make_shared<Shell>();
    shell->setClosed(true);

    auto addFace = [&](FacePtr f, bool flip) {
        if (flip) f = reverseFace(f);
        shell->addFace(f);
    };

    // ---- Process faces of A ----
    if (m_a->outerShell()) {
        for (const auto& fa : m_a->outerShell()->faces()) {
            FaceRelation rel = classifyFace(fa, m_b);

            bool keep = false;
            if (m_type == BooleanType::Union)
                keep = (rel == FaceRelation::Outside || rel == FaceRelation::Crossing);
            else if (m_type == BooleanType::Difference)
                keep = (rel == FaceRelation::Outside || rel == FaceRelation::Crossing);
            else if (m_type == BooleanType::Intersection)
                keep = (rel == FaceRelation::Inside || rel == FaceRelation::Crossing);

            if (!keep) continue;

            if (rel == FaceRelation::Crossing) {
                bool keepInside = (m_type == BooleanType::Intersection);
                auto trimmed = computeTrimmedSurface(fa, m_b, keepInside);
                if (trimmed && trimmed->hasTrimLoops()) {
                    auto tf = std::make_shared<Face>(fa->surface(), fa->outerBound());
                    tf->setTrimmedSurface(trimmed);
                    addFace(tf, false);
                } else {
                    addFace(fa, false);
                }
            } else {
                addFace(fa, false);
            }
        }
    }

    // ---- Process faces of B ----
    if (m_b->outerShell()) {
        for (const auto& fb : m_b->outerShell()->faces()) {
            FaceRelation rel = classifyFace(fb, m_a);

            bool keep = false;
            bool flip = false;

            if (m_type == BooleanType::Union) {
                keep = (rel == FaceRelation::Outside || rel == FaceRelation::Crossing);
                flip = false;
            } else if (m_type == BooleanType::Difference) {
                keep = (rel == FaceRelation::Inside || rel == FaceRelation::Crossing);
                flip = true;
            } else if (m_type == BooleanType::Intersection) {
                keep = (rel == FaceRelation::Inside || rel == FaceRelation::Crossing);
                flip = false;
            }

            if (!keep) continue;

            if (rel == FaceRelation::Crossing) {
                bool keepInside = (m_type == BooleanType::Difference ||
                                   m_type == BooleanType::Intersection);
                auto trimmed = computeTrimmedSurface(fb, m_a, keepInside);
                if (trimmed && trimmed->hasTrimLoops()) {
                    auto tf = std::make_shared<Face>(fb->surface(), fb->outerBound());
                    tf->setTrimmedSurface(trimmed);
                    addFace(tf, flip);
                } else {
                    addFace(fb, flip);
                }
            } else {
                addFace(fb, flip);
            }
        }
    }

    return std::make_shared<Solid>(shell);
}

} // namespace gector

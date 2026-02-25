#include "gector/operations/boolean_op.h"
#include "gector/nurbs/nurbs_surface.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/math/vec3.h"
#include <stdexcept>
#include <cmath>
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>

namespace gector {

BooleanOperation::BooleanOperation(SolidPtr a, SolidPtr b, BooleanType type)
    : m_a(std::move(a)), m_b(std::move(b)), m_type(type)
{
    if (!m_a || !m_b)
        throw std::invalid_argument("Boolean operands must not be null");
}

// ===========================================================================
// Internal triangle mesh type
// ===========================================================================
struct BoolTri { Vec3 v[3]; Vec3 n[3]; };
using TriMesh = std::vector<BoolTri>;

// ===========================================================================
// 2-D ear-clipping triangulation
// ===========================================================================
static double cross2d(double ax, double ay, double bx, double by) noexcept {
    return ax * by - ay * bx;
}
static bool ptInTri2d(double ax, double ay, double bx, double by,
                      double cx, double cy, double px, double py) noexcept {
    return cross2d(bx-ax,by-ay,px-ax,py-ay) > -1e-10 &&
           cross2d(cx-bx,cy-by,px-bx,py-by) > -1e-10 &&
           cross2d(ax-cx,ay-cy,px-cx,py-cy) > -1e-10;
}
static void earClip(const std::vector<Vec3>& poly, const Vec3& normal,
                    std::vector<std::array<int,3>>& out)
{
    int n = static_cast<int>(poly.size());
    if (n < 3) return;
    if (n == 3) { out.push_back({0,1,2}); return; }

    Vec3 xA = buildXAxisFromNormal(normal);
    Vec3 yA = normal.cross(xA).normalized();
    std::vector<double> px(n), py(n);
    for (int i = 0; i < n; ++i) {
        px[i] = poly[i].dot(xA);
        py[i] = poly[i].dot(yA);
    }
    double area = 0;
    for (int i = 0, j = 1; j < n; ++i, ++j)
        area += (px[j]-px[i])*(py[j]+py[i]);
    area += (px[0]-px[n-1])*(py[0]+py[n-1]);

    std::vector<int> ring(n);
    if (area > 0) for (int i = 0; i < n; ++i) ring[i] = n-1-i;
    else          std::iota(ring.begin(), ring.end(), 0);

    for (int iter = 0; static_cast<int>(ring.size()) >= 3 && iter < n*n+1; ++iter) {
        int sz = static_cast<int>(ring.size());
        bool clipped = false;
        for (int i = 0; i < sz; ++i) {
            int a = ring[(i+sz-1)%sz], b = ring[i], c = ring[(i+1)%sz];
            if (cross2d(px[b]-px[a],py[b]-py[a],px[c]-px[b],py[c]-py[b]) <= 1e-10)
                continue;
            bool ear = true;
            for (int j = 0; j < sz && ear; ++j) {
                int vv = ring[j];
                if (vv==a||vv==b||vv==c) continue;
                if (ptInTri2d(px[a],py[a],px[b],py[b],px[c],py[c],px[vv],py[vv]))
                    ear = false;
            }
            if (ear) {
                out.push_back({a,b,c});
                ring.erase(ring.begin()+i);
                clipped = true; break;
            }
        }
        if (!clipped && static_cast<int>(ring.size()) >= 3) {
            out.push_back({ring[0],ring[1],ring[2]});
            ring.erase(ring.begin()+1);
        }
    }
    if (static_cast<int>(ring.size()) == 3)
        out.push_back({ring[0],ring[1],ring[2]});
}

// ===========================================================================
// Tessellation helpers
// ===========================================================================
static constexpr int BTESS  = 14;   // grid steps for curved NURBS
static constexpr int BCSAMP = 10;   // curve samples per edge for flat faces

static void tessSurface(const NURBSSurface& s, TriMesh& out, int steps) {
    const double u0=s.uParamStart(), u1=s.uParamEnd();
    const double v0=s.vParamStart(), v1=s.vParamEnd();
    for (int j = 0; j < steps; ++j) {
        for (int i = 0; i < steps; ++i) {
            double ua=u0+(u1-u0)*i/steps,  ub=u0+(u1-u0)*(i+1)/steps;
            double va=v0+(v1-v0)*j/steps,  vb=v0+(v1-v0)*(j+1)/steps;
            Vec3 p00=s.evaluate(ua,va), p10=s.evaluate(ub,va);
            Vec3 p01=s.evaluate(ua,vb), p11=s.evaluate(ub,vb);
            Vec3 n00=s.normal(ua,va),   n10=s.normal(ub,va);
            Vec3 n01=s.normal(ua,vb),   n11=s.normal(ub,vb);
            out.push_back({{p00,p10,p11},{n00,n10,n11}});
            out.push_back({{p00,p11,p01},{n00,n11,n01}});
        }
    }
}

static void tessFlatFace(const FacePtr& face, TriMesh& out, int cSamp) {
    if (!face->outerBound()) return;
    std::vector<Vec3> poly;
    for (const auto& e : face->outerBound()->edges()) {
        auto c = e->curve();
        if (!c) { poly.push_back(e->startVertex()->position()); continue; }
        int N = (c->degree() == 1) ? 1 : cSamp;
        double t0=c->paramStart(), t1=c->paramEnd();
        for (int i = 0; i < N; ++i)
            poly.push_back(c->evaluate(t0+(t1-t0)*i/N));
    }
    while (poly.size() > 3 && (poly.back()-poly.front()).lengthSq() < 1e-12)
        poly.pop_back();
    if (static_cast<int>(poly.size()) < 3) return;

    Vec3 fn;
    if (face->surface()) {
        const auto& sf = *face->surface();
        fn = sf.normal((sf.uParamStart()+sf.uParamEnd())*.5,
                       (sf.vParamStart()+sf.vParamEnd())*.5);
    }
    if (fn.isZero()) {
        Vec3 e1=poly[1]-poly[0], e2=poly.back()-poly[0];
        fn = e1.cross(e2);
        fn = fn.isZero() ? Vec3::unitZ() : fn.normalized();
    }

    std::vector<std::array<int,3>> tris;
    earClip(poly, fn, tris);
    for (const auto& t : tris)
        out.push_back({{poly[t[0]],poly[t[1]],poly[t[2]]},{fn,fn,fn}});
}

static TriMesh buildMesh(const SolidPtr& solid) {
    TriMesh m;
    if (!solid) return m;
    auto shell = solid->outerShell();
    if (!shell) return m;
    for (const auto& face : shell->faces()) {
        auto surf = face->surface();
        bool flat = !surf || (surf->uDegree()==1 && surf->vDegree()==1);
        if (flat) tessFlatFace(face, m, BCSAMP);
        else      tessSurface(*surf, m, BTESS);
    }
    return m;
}

// If solid already has a precomputed mesh (from a prior Boolean), reuse it
// as the reference mesh rather than re-tessellating the (now-empty) B-Rep.
static TriMesh toTriMesh(const SolidPtr& solid) {
    if (solid->hasComputedMesh()) {
        TriMesh m;
        m.reserve(solid->computedMesh().size());
        for (const auto& t : solid->computedMesh()) {
            BoolTri bt;
            bt.v[0]=t.v0; bt.v[1]=t.v1; bt.v[2]=t.v2;
            bt.n[0]=t.n0; bt.n[1]=t.n1; bt.n[2]=t.n2;
            m.push_back(bt);
        }
        return m;
    }
    return buildMesh(solid);
}

// ===========================================================================
// Moller-Trumbore ray-triangle intersection
// ===========================================================================
static bool rayTri(const Vec3& o, const Vec3& d,
                   const Vec3& v0, const Vec3& v1, const Vec3& v2) noexcept {
    Vec3 e1=v1-v0, e2=v2-v0;
    Vec3 h=d.cross(e2);
    double a=e1.dot(h);
    if (std::fabs(a) < 1e-9) return false;
    double f=1.0/a;
    Vec3 s=o-v0;
    double u=f*s.dot(h); if (u<0.0||u>1.0) return false;
    Vec3 q=s.cross(e1);
    double v=f*d.dot(q); if (v<0.0||u+v>1.0) return false;
    return f*e2.dot(q) > 1e-9;
}

/// Majority-vote over three axis-aligned rays for robustness.
static bool pointInMesh(const Vec3& p, const TriMesh& mesh) noexcept {
    static const Vec3 dirs[3] = {{1,0,0},{0,1,0},{0,0,1}};
    int inside = 0;
    for (const auto& d : dirs) {
        int cnt = 0;
        for (const auto& t : mesh)
            if (rayTri(p, d, t.v[0], t.v[1], t.v[2])) ++cnt;
        if (cnt % 2 == 1) ++inside;
    }
    return inside >= 2;
}

// ===========================================================================
// Build
// ===========================================================================
SolidPtr BooleanOperation::build() const {
    // ------------------------------------------------------------------
    // Tessellate both operands into triangle meshes for classification.
    // ------------------------------------------------------------------
    const TriMesh meshA = toTriMesh(m_a);
    const TriMesh meshB = toTriMesh(m_b);

    // ------------------------------------------------------------------
    // Classify each triangle by testing its centroid against the
    // opposite mesh, then select / flip per operation type.
    // ------------------------------------------------------------------
    auto addTri = [](std::vector<Triangle>& out, const BoolTri& t, bool flip) {
        Triangle tri;
        if (!flip) {
            tri.v0=t.v[0]; tri.v1=t.v[1]; tri.v2=t.v[2];
            tri.n0=t.n[0]; tri.n1=t.n[1]; tri.n2=t.n[2];
        } else {
            // Reverse winding + negate normals = inner cavity wall
            tri.v0=t.v[0]; tri.v1=t.v[2]; tri.v2=t.v[1];
            tri.n0=-t.n[0]; tri.n1=-t.n[2]; tri.n2=-t.n[1];
        }
        out.push_back(tri);
    };

    std::vector<Triangle> result;
    result.reserve(meshA.size() + meshB.size());

    // Compute a test point for a triangle: centroid shifted slightly in the
    // direction of the triangle's average normal.  This prevents boundary
    // ambiguity when the centroid lies exactly on a face of the other mesh
    // (e.g. two boxes sharing the z=0 plane).
    auto testPt = [](const BoolTri& t) -> Vec3 {
        Vec3 c = (t.v[0]+t.v[1]+t.v[2]) / 3.0;
        Vec3 avgN = (t.n[0]+t.n[1]+t.n[2]) / 3.0;
        double len = avgN.length();
        if (len > 1e-12) c += avgN * (1e-5 / len);
        return c;
    };

    switch (m_type) {

        case BooleanType::Union:
            // A-outside-B  +  B-outside-A
            for (const auto& t : meshA) {
                if (!pointInMesh(testPt(t), meshB)) addTri(result, t, false);
            }
            for (const auto& t : meshB) {
                if (!pointInMesh(testPt(t), meshA)) addTri(result, t, false);
            }
            break;

        case BooleanType::Difference:
            // A-outside-B  +  B-inside-A (inner wall, normals flipped)
            for (const auto& t : meshA) {
                if (!pointInMesh(testPt(t), meshB)) addTri(result, t, false);
            }
            for (const auto& t : meshB) {
                if (pointInMesh(testPt(t), meshA)) addTri(result, t, true);
            }
            break;

        case BooleanType::Intersection:
            // A-inside-B  +  B-inside-A
            for (const auto& t : meshA) {
                if (pointInMesh(testPt(t), meshB)) addTri(result, t, false);
            }
            for (const auto& t : meshB) {
                if (pointInMesh(testPt(t), meshA)) addTri(result, t, false);
            }
            break;
    }

    // ------------------------------------------------------------------
    // Package result.  The B-Rep outer shell is left empty (geometry is
    // fully captured in the precomputed mesh). ObjExporter detects the
    // non-empty computedMesh and writes it directly.
    // ------------------------------------------------------------------
    auto shell = std::make_shared<Shell>();
    shell->setClosed(true);
    auto solid = std::make_shared<Solid>(shell);
    solid->setComputedMesh(std::move(result));
    return solid;
}

} // namespace gector

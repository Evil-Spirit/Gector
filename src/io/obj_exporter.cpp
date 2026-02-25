#include "gector/io/obj_exporter.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"
#include "gector/math/vec3.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <array>
#include <numeric>
#include <cmath>
#include <algorithm>

namespace gector {

// ---------------------------------------------------------------------------
// Internal stream helpers
// ---------------------------------------------------------------------------
static void writeVertex(std::ostream& os, const Point3D& p) {
    os << "v " << std::fixed << std::setprecision(6)
       << p.x << ' ' << p.y << ' ' << p.z << '\n';
}

static void writeNormal(std::ostream& os, const Vec3& n) {
    os << "vn " << std::fixed << std::setprecision(6)
       << n.x << ' ' << n.y << ' ' << n.z << '\n';
}

// Write a triangle.  flip=true reverses winding (for void/inner shells).
// Format: f v1//n1 v2//n2 v3//n3
static void writeFaceTri(std::ostream& os, int v1, int v2, int v3, bool flip) {
    if (flip) std::swap(v2, v3);
    os << "f " << v1 << "//" << v1
       << ' '  << v2 << "//" << v2
       << ' '  << v3 << "//" << v3 << '\n';
}

// ---------------------------------------------------------------------------
// Ear-clipping polygon triangulation (§3.2 of "Computational Geometry")
// ---------------------------------------------------------------------------

// 2D signed cross-product (z-component).
static double cross2d(double ax, double ay, double bx, double by) {
    return ax * by - ay * bx;
}

// Returns true if point P strictly lies inside the CCW triangle ABC.
static bool pointInEar(double ax, double ay,
                        double bx, double by,
                        double cx, double cy,
                        double px, double py)
{
    double d1 = cross2d(bx-ax, by-ay, px-ax, py-ay);
    double d2 = cross2d(cx-bx, cy-by, px-bx, py-by);
    double d3 = cross2d(ax-cx, ay-cy, px-cx, py-cy);
    // Point is inside if all three half-plane tests agree (all ≥ 0)
    return d1 > -1e-10 && d2 > -1e-10 && d3 > -1e-10;
}

// Ear-clipping triangulation of a planar polygon.
// pts: 3D polygon vertices assumed to be planar.
// normal: face normal used to project to 2D.
// Returns a list of {a, b, c} index triplets (CCW relative to normal).
/*static*/ std::vector<std::array<int,3>>
ObjExporter::triangulatePolygon(const std::vector<Vec3>& pts,
                                  const Vec3& normal)
{
    std::vector<std::array<int,3>> tris;
    int n = static_cast<int>(pts.size());
    if (n < 3) return tris;
    if (n == 3) { tris.push_back({0, 1, 2}); return tris; }

    // Build an orthonormal frame in the face plane for 2D projection.
    Vec3 xAxis = buildXAxisFromNormal(normal);
    Vec3 yAxis = normal.cross(xAxis).normalized();

    std::vector<double> px(n), py(n);
    for (int i = 0; i < n; ++i) {
        px[i] = pts[i].dot(xAxis);
        py[i] = pts[i].dot(yAxis);
    }

    // Compute signed area in 2D (positive = CCW).
    double area2 = 0.0;
    for (int i = 0, j = 1; j < n; ++i, ++j)
        area2 += (px[j] - px[i]) * (py[j] + py[i]);
    area2 += (px[0] - px[n-1]) * (py[0] + py[n-1]);

    // Build vertex ring; reverse if the polygon is CW in 2D.
    std::vector<int> ring(n);
    if (area2 > 0.0) {
        for (int i = 0; i < n; ++i) ring[i] = n - 1 - i;
    } else {
        std::iota(ring.begin(), ring.end(), 0);
    }

    // Ear-clipping main loop.
    const int maxIter = n * n + n + 1;
    for (int iter = 0; static_cast<int>(ring.size()) >= 3 && iter < maxIter; ++iter) {
        int sz = static_cast<int>(ring.size());
        bool clipped = false;
        for (int i = 0; i < sz; ++i) {
            int a = ring[(i + sz - 1) % sz];
            int b = ring[i];
            int c = ring[(i + 1) % sz];

            // Check convex at b (CCW turn: cross > 0).
            double ax2 = px[b] - px[a], ay2 = py[b] - py[a];
            double bx2 = px[c] - px[b], by2 = py[c] - py[b];
            if (cross2d(ax2, ay2, bx2, by2) <= 1e-10) continue;

            // Check no other polygon vertex inside this ear triangle.
            bool ear = true;
            for (int j = 0; j < sz && ear; ++j) {
                int v = ring[j];
                if (v == a || v == b || v == c) continue;
                if (pointInEar(px[a], py[a], px[b], py[b], px[c], py[c],
                               px[v], py[v]))
                    ear = false;
            }
            if (ear) {
                tris.push_back({a, b, c});
                ring.erase(ring.begin() + i);
                clipped = true;
                break;
            }
        }
        // Degenerate case fallback (collinear / near-collinear polygon).
        if (!clipped && ring.size() >= 3) {
            tris.push_back({ring[0], ring[1], ring[2]});
            ring.erase(ring.begin() + 1);
        }
    }
    if (ring.size() == 3)
        tris.push_back({ring[0], ring[1], ring[2]});

    return tris;
}

// ---------------------------------------------------------------------------
// isFlatFace
// ---------------------------------------------------------------------------
bool ObjExporter::isFlatFace(const FacePtr& face) const noexcept {
    auto surf = face->surface();
    if (!surf) return true;                              // no surface → use boundary
    if (surf->uDegree() == 1 && surf->vDegree() == 1) return true; // bilinear plane
    return false;
}

// ---------------------------------------------------------------------------
// sampleBoundaryPolygon
// Degree-1 (line) curves: 1 sample (just the start vertex = corner).
// Higher-degree curves: curveSamples samples.
// The last sample of each edge is skipped to avoid duplication with the
// next edge's first sample.
// ---------------------------------------------------------------------------
std::vector<Vec3> ObjExporter::sampleBoundaryPolygon(const FacePtr& face) const {
    std::vector<Vec3> poly;
    if (!face->outerBound()) return poly;

    for (const auto& edge : face->outerBound()->edges()) {
        auto curve = edge->curve();
        if (!curve) {
            poly.push_back(edge->startVertex()->position());
            continue;
        }
        const double t0 = curve->paramStart();
        const double t1 = curve->paramEnd();
        // Degree-1 curves need only the start vertex.
        const int N = (curve->degree() == 1) ? 1 : m_curveSamples;
        for (int i = 0; i < N; ++i) {
            double t = t0 + (t1 - t0) * static_cast<double>(i) / N;
            poly.push_back(curve->evaluate(t));
        }
    }
    return poly;
}

// ---------------------------------------------------------------------------
// tessellatePlanarCap – ear-clipping triangulation of boundary polygon
// ---------------------------------------------------------------------------
int ObjExporter::tessellatePlanarCap(const FacePtr& face,
                                      std::ostream& os,
                                      int vOffset,
                                      bool flip) const {
    auto poly = sampleBoundaryPolygon(face);

    // Drop exact duplicate closing vertex (closed loop).
    while (poly.size() > 3 &&
           (poly.back() - poly.front()).lengthSq() < 1e-12)
        poly.pop_back();

    // Fallback: if the boundary polygon is degenerate (< 3 points), the
    // face boundary wire is incomplete (e.g. an extrusion side face whose
    // wire only has top and bottom edges but not the two lateral connectors).
    // In this case the underlying NURBS surface still covers the full patch,
    // so use the surface grid sampler instead.
    if (static_cast<int>(poly.size()) < 3) {
        if (face->surface())
            return tessellateNURBSFace(face, os, vOffset, flip);
        return 0;
    }

    // Compute face normal.
    Vec3 faceNormal;
    if (face->surface()) {
        auto surf = face->surface();
        Vec3 sn = surf->normal(
            (surf->uParamStart() + surf->uParamEnd()) * 0.5,
            (surf->vParamStart() + surf->vParamEnd()) * 0.5);
        if (!sn.isZero()) faceNormal = sn;
    }
    if (faceNormal.isZero()) {
        // Estimate from polygon cross-product.
        Vec3 e1 = poly[1] - poly[0];
        Vec3 e2 = poly.back() - poly[0];
        Vec3 c  = e1.cross(e2);
        faceNormal = c.isZero() ? Vec3::unitZ() : c.normalized();
    }
    if (flip) faceNormal = -faceNormal;

    // Write all polygon vertices and their shared normal.
    const int nPoly = static_cast<int>(poly.size());
    for (const auto& p : poly) {
        writeVertex(os, p);
        writeNormal(os, faceNormal);
    }

    // Ear-clipping → triangle list.
    // When flip=true the triangulatePolygon result is CCW in the plane
    // opposite to `flip`; writeFaceTri re-reverses winding to make them
    // face outward relative to the (negated) normal.
    auto tris = triangulatePolygon(poly, flip ? -faceNormal : faceNormal);
    for (const auto& t : tris) {
        writeFaceTri(os,
                     vOffset + t[0],
                     vOffset + t[1],
                     vOffset + t[2],
                     flip);
    }

    return nPoly;
}

// ---------------------------------------------------------------------------
// tessellateNURBSFace – uniform parameter-space grid
// ---------------------------------------------------------------------------
int ObjExporter::tessellateNURBSFace(const FacePtr& face,
                                      std::ostream& os,
                                      int vOffset,
                                      bool flip) const {
    auto surf = face->surface();
    if (!surf) return 0;

    const int nu = m_uSteps;
    const int nv = m_vSteps;
    const double u0 = surf->uParamStart();
    const double u1 = surf->uParamEnd();
    const double v0 = surf->vParamStart();
    const double v1 = surf->vParamEnd();

    // Sample (nu+1) × (nv+1) points.
    std::vector<Point3D> pts;
    std::vector<Vec3>    nrm;
    pts.reserve(static_cast<std::size_t>((nu+1)*(nv+1)));
    nrm.reserve(static_cast<std::size_t>((nu+1)*(nv+1)));

    for (int j = 0; j <= nv; ++j) {
        double v = v0 + (v1 - v0) * static_cast<double>(j) / nv;
        for (int i = 0; i <= nu; ++i) {
            double u = u0 + (u1 - u0) * static_cast<double>(i) / nu;
            pts.push_back(surf->evaluate(u, v));
            Vec3 n = surf->normal(u, v);
            nrm.push_back(flip ? -n : n);
        }
    }

    for (std::size_t k = 0; k < pts.size(); ++k) {
        writeVertex(os, pts[k]);
        writeNormal(os, nrm[k]);
    }

    // Triangulate quads (CCW winding for front face).
    auto idx = [&](int i, int j) -> int {
        return vOffset + j * (nu + 1) + i;  // 1-based
    };
    for (int j = 0; j < nv; ++j) {
        for (int i = 0; i < nu; ++i) {
            int v00 = idx(i,   j);
            int v10 = idx(i+1, j);
            int v01 = idx(i,   j+1);
            int v11 = idx(i+1, j+1);
            writeFaceTri(os, v00, v10, v11, flip);
            writeFaceTri(os, v00, v11, v01, flip);
        }
    }

    return static_cast<int>(pts.size());
}

// ---------------------------------------------------------------------------
// writeSolidToStream
// ---------------------------------------------------------------------------
void ObjExporter::writeSolidToStream(const SolidPtr& solid,
                                      const std::string& objectName,
                                      std::ostream& os,
                                      int& vertexOffset) const {
    if (!solid) return;

    os << "o " << objectName << "\n\n";

    auto tessellateShell = [&](const ShellPtr& shell,
                                const std::string& prefix,
                                bool flip) {
        if (!shell) return;
        int faceIdx = 0;
        for (const auto& face : shell->faces()) {
            ++faceIdx;
            os << "g " << prefix << "_face" << faceIdx << '\n';
            int added = 0;
            if (isFlatFace(face))
                added = tessellatePlanarCap(face, os, vertexOffset, flip);
            else
                added = tessellateNURBSFace(face, os, vertexOffset, flip);
            vertexOffset += added;
        }
    };

    tessellateShell(solid->outerShell(), objectName + "_outer", false);
    // Void shells (e.g. the B operand of Difference) are rendered flipped
    // so they appear as interior cavity surfaces.
    for (std::size_t i = 0; i < solid->voids().size(); ++i)
        tessellateShell(solid->voids()[i],
                        objectName + "_void" + std::to_string(i),
                        /*flip=*/true);

    os << '\n';
}

// ---------------------------------------------------------------------------
// writeSolid
// ---------------------------------------------------------------------------
bool ObjExporter::writeSolid(const SolidPtr& solid,
                               const std::string& filename) const {
    if (!solid) return false;

    std::ofstream ofs(filename);
    if (!ofs.is_open()) return false;

    ofs << "# Gector NURBS solid export\n";
    ofs << "# " << filename << "\n\n";

    int offset = 1;
    writeSolidToStream(solid, "solid", ofs, offset);
    return ofs.good();
}

// ---------------------------------------------------------------------------
// writeSolids
// ---------------------------------------------------------------------------
bool ObjExporter::writeSolids(const std::vector<SolidPtr>& solids,
                                const std::string& filename) const {
    if (solids.empty()) return false;

    std::ofstream ofs(filename);
    if (!ofs.is_open()) return false;

    ofs << "# Gector NURBS solid export\n";
    ofs << "# " << filename << "\n\n";

    int offset = 1;
    for (std::size_t i = 0; i < solids.size(); ++i) {
        const std::string name = "solid" + std::to_string(i + 1);
        writeSolidToStream(solids[i], name, ofs, offset);
    }
    return ofs.good();
}

} // namespace gector

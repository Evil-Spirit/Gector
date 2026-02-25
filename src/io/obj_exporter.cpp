#include "gector/io/obj_exporter.h"
#include "gector/nurbs/nurbs_surface.h"
#include "gector/math/vec3.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

namespace gector {

// ---------------------------------------------------------------------------
// Stream helpers
// ---------------------------------------------------------------------------
static void writeVertex(std::ostream& os, const Point3D& p) {
    os << "v " << std::fixed << std::setprecision(6)
       << p.x << ' ' << p.y << ' ' << p.z << '\n';
}

static void writeNormal(std::ostream& os, const Vec3& n) {
    os << "vn " << std::fixed << std::setprecision(6)
       << n.x << ' ' << n.y << ' ' << n.z << '\n';
}

// Write a triangular face with 1-based OBJ indices: f v1//n1 v2//n2 v3//n3
static void writeFaceTri(std::ostream& os, int v1, int v2, int v3) {
    os << "f " << v1 << "//" << v1
       << ' '  << v2 << "//" << v2
       << ' '  << v3 << "//" << v3 << '\n';
}

// ---------------------------------------------------------------------------
// isFlatFace – decide whether to use planar-cap tessellation
// ---------------------------------------------------------------------------
bool ObjExporter::isFlatFace(const FacePtr& face) const noexcept {
    auto surf = face->surface();
    // No surface: definitely use boundary curves
    if (!surf) return true;
    // Degree-1 in both directions → bilinear (flat) plane patch
    if (surf->uDegree() == 1 && surf->vDegree() == 1) return true;
    return false;
}

// ---------------------------------------------------------------------------
// sampleBoundaryPolygon – walk all boundary edges and sample their curves
// ---------------------------------------------------------------------------
std::vector<Vec3> ObjExporter::sampleBoundaryPolygon(const FacePtr& face) const {
    std::vector<Vec3> poly;
    if (!face->outerBound()) return poly;

    const int N = m_curveSamples;
    for (const auto& edge : face->outerBound()->edges()) {
        auto curve = edge->curve();
        if (!curve) {
            // Degenerate: just take the start vertex
            poly.push_back(edge->startVertex()->position());
            continue;
        }
        const double t0 = curve->paramStart();
        const double t1 = curve->paramEnd();
        // Sample N interior segments; skip the very last point to avoid
        // duplication with the next edge's first point.
        for (int i = 0; i < N; ++i) {
            double t = t0 + (t1 - t0) * static_cast<double>(i) / N;
            poly.push_back(curve->evaluate(t));
        }
    }
    return poly;
}

// ---------------------------------------------------------------------------
// tessellatePlanarCap – fan-triangulate from centroid
// ---------------------------------------------------------------------------
int ObjExporter::tessellatePlanarCap(const FacePtr& face,
                                      std::ostream& os,
                                      int vOffset) const {
    auto poly = sampleBoundaryPolygon(face);
    if (poly.size() < 3) return 0;

    // Compute centroid
    Vec3 centroid;
    for (const auto& p : poly) centroid += p;
    centroid /= static_cast<double>(poly.size());

    // Approximate face normal from first edge cross product
    Vec3 edge1 = poly[1] - poly[0];
    Vec3 edge2 = poly.back() - poly[0];
    Vec3 faceNormal = edge1.cross(edge2);
    if (!faceNormal.isZero()) faceNormal = faceNormal.normalized();
    else                       faceNormal = Vec3::unitZ();

    // If the face has a surface with a computable normal, use that instead
    if (face->surface()) {
        auto surf = face->surface();
        Vec3 sn = surf->normal(
            (surf->uParamStart() + surf->uParamEnd()) * 0.5,
            (surf->vParamStart() + surf->vParamEnd()) * 0.5);
        if (!sn.isZero()) faceNormal = sn;
    }

    // Write centroid + polygon vertices
    writeVertex(os, centroid);
    writeNormal(os, faceNormal);
    int centIdx = vOffset; // 1-based
    int nPoly   = static_cast<int>(poly.size());
    for (const auto& p : poly) {
        writeVertex(os, p);
        writeNormal(os, faceNormal);
    }

    // Fan triangles: centroid → poly[i] → poly[(i+1)%n]
    for (int i = 0; i < nPoly; ++i) {
        int a = centIdx;
        int b = vOffset + 1 + i;
        int c = vOffset + 1 + (i + 1) % nPoly;
        writeFaceTri(os, a, b, c);
    }

    return 1 + nPoly; // centroid + polygon vertices
}

// ---------------------------------------------------------------------------
// tessellateNURBSFace – sample curved surface as a uniform grid
// ---------------------------------------------------------------------------
int ObjExporter::tessellateNURBSFace(const FacePtr& face,
                                      std::ostream& os,
                                      int vOffset) const {
    auto surf = face->surface();
    if (!surf) return 0;

    const int nu = m_uSteps;
    const int nv = m_vSteps;
    const double u0 = surf->uParamStart();
    const double u1 = surf->uParamEnd();
    const double v0 = surf->vParamStart();
    const double v1 = surf->vParamEnd();

    // Sample (nu+1) × (nv+1) points
    std::vector<Point3D> pts;
    std::vector<Vec3>    nrm;
    pts.reserve((nu+1)*(nv+1));
    nrm.reserve((nu+1)*(nv+1));

    for (int j = 0; j <= nv; ++j) {
        double v = v0 + (v1 - v0) * static_cast<double>(j) / nv;
        for (int i = 0; i <= nu; ++i) {
            double u = u0 + (u1 - u0) * static_cast<double>(i) / nu;
            pts.push_back(surf->evaluate(u, v));
            nrm.push_back(surf->normal(u, v));
        }
    }

    // Emit vertices + normals
    for (std::size_t k = 0; k < pts.size(); ++k) {
        writeVertex(os, pts[k]);
        writeNormal(os, nrm[k]);
    }

    // Emit triangulated quads (CCW)
    auto idx = [&](int i, int j) -> int {
        return vOffset + j * (nu + 1) + i;  // 1-based
    };
    for (int j = 0; j < nv; ++j) {
        for (int i = 0; i < nu; ++i) {
            int v00 = idx(i,   j);
            int v10 = idx(i+1, j);
            int v01 = idx(i,   j+1);
            int v11 = idx(i+1, j+1);
            writeFaceTri(os, v00, v10, v11);
            writeFaceTri(os, v00, v11, v01);
        }
    }

    return static_cast<int>(pts.size());
}

// ---------------------------------------------------------------------------
// Write a single solid to a stream
// ---------------------------------------------------------------------------
void ObjExporter::writeSolidToStream(const SolidPtr& solid,
                                      const std::string& objectName,
                                      std::ostream& os,
                                      int& vertexOffset) const {
    if (!solid) return;

    os << "o " << objectName << "\n\n";

    auto tessellateShell = [&](const ShellPtr& shell, const std::string& prefix) {
        if (!shell) return;
        int faceIdx = 0;
        for (const auto& face : shell->faces()) {
            ++faceIdx;
            os << "g " << prefix << "_face" << faceIdx << '\n';

            int added = 0;
            if (isFlatFace(face)) {
                added = tessellatePlanarCap(face, os, vertexOffset);
            } else {
                added = tessellateNURBSFace(face, os, vertexOffset);
            }
            vertexOffset += added;
        }
    };

    tessellateShell(solid->outerShell(), objectName + "_outer");
    for (std::size_t i = 0; i < solid->voids().size(); ++i)
        tessellateShell(solid->voids()[i],
                        objectName + "_void" + std::to_string(i));

    os << '\n';
}

// ---------------------------------------------------------------------------
// Public: write single solid to file
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
// Public: write multiple solids to file
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

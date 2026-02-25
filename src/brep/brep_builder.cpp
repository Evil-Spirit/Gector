#include "gector/brep/brep_builder.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"
#include <cmath>
#include <stdexcept>

namespace gector {

// ---------------------------------------------------------------------------
// Low-level helpers
// ---------------------------------------------------------------------------
EdgePtr BRepBuilder::makeLinearEdge(const Point3D& p0, const Point3D& p1) {
    auto v0 = std::make_shared<Vertex>(p0);
    auto v1 = std::make_shared<Vertex>(p1);
    auto curve = std::make_shared<NURBSCurve>(NURBSCurve::makeLine(p0, p1));
    return std::make_shared<Edge>(v0, v1, curve);
}

EdgePtr BRepBuilder::makeEdge(VertexPtr v0, VertexPtr v1,
                               std::shared_ptr<NURBSCurve> curve) {
    return std::make_shared<Edge>(std::move(v0), std::move(v1), std::move(curve));
}

FacePtr BRepBuilder::makePlanarFace(WirePtr boundary, const Vec3& normal) {
    // Create a large enough plane patch through the wire's centroid.
    // The exact extent does not matter for topology; evaluation uses the surface.
    Vec3 centroid;
    auto verts = boundary->vertices();
    for (const auto& v : verts) centroid += v->position();
    if (!verts.empty()) centroid /= static_cast<double>(verts.size());

    // Build a coordinate frame in the plane
    Vec3 n = normal.isZero() ? Vec3::unitZ() : normal.normalized();
    Vec3 u;
    if (std::abs(n.dot(Vec3::unitX())) < 0.9)
        u = n.cross(Vec3::unitX()).normalized();
    else
        u = n.cross(Vec3::unitY()).normalized();
    Vec3 v = n.cross(u).normalized();

    auto surf = std::make_shared<NURBSSurface>(
        NURBSSurface::makePlane(centroid - u * 1e3 - v * 1e3, u, v, 2e3, 2e3));
    return std::make_shared<Face>(surf, std::move(boundary));
}

FacePtr BRepBuilder::makeRuledFace(EdgePtr bottomEdge, EdgePtr topEdge) {
    auto botCurve = bottomEdge->curve();
    auto topCurve = topEdge->curve();
    if (!botCurve || !topCurve)
        throw std::invalid_argument("Both edges must have a curve geometry");

    auto surf = std::make_shared<NURBSSurface>(
        NURBSSurface::makeRuledSurface(*botCurve, *topCurve));

    // Build a wire: bot + right + top(reversed) + left
    auto w = std::make_shared<Wire>();
    w->addEdge(bottomEdge);
    w->addEdge(topEdge);
    return std::make_shared<Face>(surf, w);
}

// ---------------------------------------------------------------------------
// Box
// ---------------------------------------------------------------------------
SolidPtr BRepBuilder::makeBox(double sx, double sy, double sz) {
    if (sx <= 0 || sy <= 0 || sz <= 0)
        throw std::invalid_argument("Box dimensions must be positive");

    // 8 vertices
    auto V = [](double x, double y, double z){ return std::make_shared<Vertex>(Point3D{x,y,z}); };
    auto v000 = V(0,0,0);  auto v100 = V(sx,0,0);
    auto v110 = V(sx,sy,0);auto v010 = V(0,sy,0);
    auto v001 = V(0,0,sz); auto v101 = V(sx,0,sz);
    auto v111 = V(sx,sy,sz);auto v011 = V(0,sy,sz);

    // Helper: make a closed rectangular wire from 4 vertices (in order)
    auto makeRect = [](VertexPtr a, VertexPtr b, VertexPtr c, VertexPtr d) {
        auto w = std::make_shared<Wire>();
        auto addEdge = [&](VertexPtr p, VertexPtr q) {
            auto curve = std::make_shared<NURBSCurve>(
                NURBSCurve::makeLine(p->position(), q->position()));
            w->addEdge(std::make_shared<Edge>(p, q, curve));
        };
        addEdge(a,b); addEdge(b,c); addEdge(c,d); addEdge(d,a);
        return w;
    };

    // Helper: make a planar face
    auto makeFace = [&](WirePtr wire, Vec3 normal) {
        return makePlanarFace(wire, normal);
    };

    auto shell = std::make_shared<Shell>();

    // Bottom (z=0, normal -Z)
    shell->addFace(makeFace(makeRect(v000,v010,v110,v100), -Vec3::unitZ()));
    // Top (z=sz, normal +Z)
    shell->addFace(makeFace(makeRect(v001,v101,v111,v011),  Vec3::unitZ()));
    // Front (y=0, normal -Y)
    shell->addFace(makeFace(makeRect(v000,v100,v101,v001), -Vec3::unitY()));
    // Back (y=sy, normal +Y)
    shell->addFace(makeFace(makeRect(v010,v011,v111,v110),  Vec3::unitY()));
    // Left (x=0, normal -X)
    shell->addFace(makeFace(makeRect(v000,v001,v011,v010), -Vec3::unitX()));
    // Right (x=sx, normal +X)
    shell->addFace(makeFace(makeRect(v100,v110,v111,v101),  Vec3::unitX()));

    shell->setClosed(true);
    return std::make_shared<Solid>(shell);
}

// ---------------------------------------------------------------------------
// Cylinder
// ---------------------------------------------------------------------------
SolidPtr BRepBuilder::makeCylinder(const Point3D& base,
                                    const Vec3& axisDir,
                                    double radius, double height) {
    if (radius <= 0 || height <= 0)
        throw std::invalid_argument("Cylinder radius and height must be positive");

    Vec3 axis = axisDir.normalized();

    // Side surface
    auto sideSurf = std::make_shared<NURBSSurface>(
        NURBSSurface::makeCylinder(base, axis, radius, height));

    // Build circle wires for top and bottom faces
    auto makeCircleWire = [&](const Point3D& center, double r, const Vec3& n) {
        auto curve = std::make_shared<NURBSCurve>(NURBSCurve::makeCircle(center, r, n));
        auto v = std::make_shared<Vertex>(center + Vec3{r, 0, 0});
        auto w = std::make_shared<Wire>();
        w->addEdge(std::make_shared<Edge>(v, v, curve));
        return w;
    };

    Vec3 normal = axis;
    auto botWire = makeCircleWire(base, radius, -normal);
    auto topPt   = base + axis * height;
    auto topWire = makeCircleWire(topPt, radius,  normal);

    auto botFace  = makePlanarFace(botWire, -normal);
    auto topFace  = makePlanarFace(topWire,  normal);

    // Side face (we create a boundary wire for it)
    auto sideWire = std::make_shared<Wire>();
    auto shell = std::make_shared<Shell>();
    shell->addFace(botFace);
    shell->addFace(topFace);
    // Side face with the cylinder surface
    shell->addFace(std::make_shared<Face>(sideSurf, sideWire));
    shell->setClosed(true);
    return std::make_shared<Solid>(shell);
}

// ---------------------------------------------------------------------------
// Sphere
// ---------------------------------------------------------------------------
SolidPtr BRepBuilder::makeSphere(const Point3D& center, double radius) {
    if (radius <= 0)
        throw std::invalid_argument("Sphere radius must be positive");

    auto surf = std::make_shared<NURBSSurface>(NURBSSurface::makeSphere(center, radius));
    auto wire = std::make_shared<Wire>(); // seam / degenerate boundary
    auto face  = std::make_shared<Face>(surf, wire);
    auto shell = std::make_shared<Shell>();
    shell->addFace(face);
    shell->setClosed(true);
    return std::make_shared<Solid>(shell);
}

// ---------------------------------------------------------------------------
// Cone
// ---------------------------------------------------------------------------
SolidPtr BRepBuilder::makeCone(const Point3D& base,
                                const Vec3& axisDir,
                                double bottomRadius, double topRadius,
                                double height) {
    if (height <= 0 || bottomRadius < 0 || topRadius < 0)
        throw std::invalid_argument("Invalid cone parameters");

    Vec3 axis = axisDir.normalized();
    auto coneSurf = std::make_shared<NURBSSurface>(
        NURBSSurface::makeCone(base, axis, bottomRadius, topRadius, height));

    auto makeCircleWire = [&](const Point3D& center, double r, const Vec3& n) {
        auto curve = std::make_shared<NURBSCurve>(NURBSCurve::makeCircle(center, r, n));
        auto v = std::make_shared<Vertex>(center + Vec3{r, 0, 0});
        auto w = std::make_shared<Wire>();
        w->addEdge(std::make_shared<Edge>(v, v, curve));
        return w;
    };

    auto shell = std::make_shared<Shell>();

    if (bottomRadius > 1e-10) {
        auto botWire = makeCircleWire(base, bottomRadius, -axis);
        shell->addFace(makePlanarFace(botWire, -axis));
    }
    if (topRadius > 1e-10) {
        auto topWire = makeCircleWire(base + axis*height, topRadius, axis);
        shell->addFace(makePlanarFace(topWire, axis));
    }

    auto sideWire = std::make_shared<Wire>();
    shell->addFace(std::make_shared<Face>(coneSurf, sideWire));
    shell->setClosed(true);
    return std::make_shared<Solid>(shell);
}

} // namespace gector

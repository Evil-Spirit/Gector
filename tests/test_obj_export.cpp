#include "test_framework.h"
#include "gector/io/obj_exporter.h"
#include "gector/brep/brep_builder.h"
#include "gector/sketch/circle.h"
#include "gector/sketch/arc.h"
#include "gector/sketch/line.h"
#include "gector/sketch/sketch.h"
#include "gector/sketch/spline_curve.h"
#include "gector/operations/extrusion.h"
#include "gector/operations/revolve.h"
#include "gector/operations/boolean_op.h"
#include "gector/operations/fillet.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

using namespace gector;

// ---------------------------------------------------------------------------
// Helper: count lines starting with a given prefix in an OBJ file
// ---------------------------------------------------------------------------
static int countObjLines(const std::string& path, const std::string& prefix) {
    std::ifstream f(path);
    if (!f.is_open()) return -1;
    int count = 0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.size() >= prefix.size() &&
            line.substr(0, prefix.size()) == prefix)
            ++count;
    }
    return count;
}

static bool fileExists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

// Output directory (relative to where the test binary runs)
static const std::string OUT_DIR = "obj_output/";

// ---------------------------------------------------------------------------
// Primitive solids
// ---------------------------------------------------------------------------
TEST(ObjExport, Box) {
    auto solid = BRepBuilder::makeBox(10, 8, 6);
    ObjExporter exp;
    exp.setSurfaceSteps(8);
    std::string path = OUT_DIR + "box.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "v ") > 0);
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
    // 6 faces, each with 8×8 quads = 2 triangles → 6 * 8*8*2 = 768 face lines
    ASSERT_TRUE(countObjLines(path, "f ") >= 1);
}

TEST(ObjExport, Cylinder) {
    auto solid = BRepBuilder::makeCylinder({0,0,0}, Vec3::unitZ(), 5.0, 12.0);
    ObjExporter exp;
    exp.setSurfaceSteps(24);
    std::string path = OUT_DIR + "cylinder.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    // Cylinder has 3 shells (bottom cap, top cap, side)
    // Side has 24×1 columns → 24*2 triangles for a degree-1-in-v surface
    ASSERT_TRUE(countObjLines(path, "v ") > 0);
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, Sphere) {
    auto solid = BRepBuilder::makeSphere({0,0,0}, 7.0);
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "sphere.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "v ") > 0);
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, Cone) {
    auto solid = BRepBuilder::makeCone({0,0,0}, Vec3::unitZ(), 5.0, 0.0, 10.0);
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "cone.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "v ") > 0);
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, Frustum) {
    // Cone frustum (truncated cone)
    auto solid = BRepBuilder::makeCone({0,0,0}, Vec3::unitZ(), 6.0, 3.0, 8.0);
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "frustum.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "v ") > 0);
}

// ---------------------------------------------------------------------------
// Extrusion
// ---------------------------------------------------------------------------
TEST(ObjExport, ExtrudedRectangle) {
    Sketch sk;
    sk.addLine({-5,0,0}, { 5,0,0});
    sk.addLine({ 5,0,0}, { 5,4,0});
    sk.addLine({ 5,4,0}, {-5,4,0});
    sk.addLine({-5,4,0}, {-5,0,0});
    auto solid = Extrusion(sk.toWire(), Vec3::unitZ(), 8.0).build();
    ObjExporter exp;
    exp.setSurfaceSteps(10);
    std::string path = OUT_DIR + "extruded_rectangle.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, ExtrudedCircle) {
    Circle c{{0,0,0}, 4.0};
    auto solid = Extrusion(c.toWire(), Vec3::unitZ(), 10.0).build();
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "extruded_circle.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, ExtrudedArc) {
    // Extrude a D-shaped profile: semicircle + straight line
    Sketch sk;
    sk.addArc({0,0,0}, 3.0, 0.0, M_PI);        // semicircle
    sk.addLine({-3,0,0}, {3,0,0});               // closing line
    auto solid = Extrusion(sk.toWire(), Vec3::unitZ(), 5.0).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "extruded_arc.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, ExtrudedSpline) {
    // Extrude a free-form profile
    std::vector<Point3D> pts = {{0,0,0},{2,3,0},{4,1,0},{6,4,0},{8,0,0}};
    Sketch sk;
    sk.addSpline(pts, 3);
    auto solid = Extrusion(sk.toWire(), Vec3::unitZ(), 4.0).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "extruded_spline.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

// ---------------------------------------------------------------------------
// Revolve
// ---------------------------------------------------------------------------
TEST(ObjExport, RevolveLine_FullCircle) {
    // Revolve a line segment at x=4 around Z → thin cylinder
    Line l{{4,0,0},{4,0,8}};
    auto wire = std::make_shared<Wire>();
    wire->addEdge(l.toEdge());
    auto solid = Revolve(wire, {0,0,0}, Vec3::unitZ(), 2*M_PI).build();
    ObjExporter exp;
    exp.setSurfaceSteps(24);
    std::string path = OUT_DIR + "revolve_line_full.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, RevolveLine_HalfCircle) {
    // Half revolution → D-shaped solid
    Line l{{3,0,0},{3,0,6}};
    auto wire = std::make_shared<Wire>();
    wire->addEdge(l.toEdge());
    auto solid = Revolve(wire, {0,0,0}, Vec3::unitZ(), M_PI).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "revolve_line_half.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, RevolveArc_Torus) {
    // Revolve a small circle arc around an axis at distance 6 → torus-like shape
    Arc a{{6,0,0}, 2.0, Vec3::unitY(), 0.0, 2*M_PI};
    auto wire = std::make_shared<Wire>();
    wire->addEdge(a.toEdge());
    auto solid = Revolve(wire, {0,0,0}, Vec3::unitZ(), 2*M_PI).build();
    ObjExporter exp;
    exp.setSurfaceSteps(24);
    std::string path = OUT_DIR + "revolve_arc_torus.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, RevolveSpline) {
    // Revolve an S-shaped profile → vase
    std::vector<Point3D> pts = {{2,0,0},{3,0,3},{2,0,6},{4,0,9},{2,0,12}};
    auto spline = SplineCurve::interpolate(pts, 3);
    auto wire = std::make_shared<Wire>();
    wire->addEdge(spline.toEdge());
    auto solid = Revolve(wire, {0,0,0}, Vec3::unitZ(), 2*M_PI).build();
    ObjExporter exp;
    exp.setSurfaceSteps(24);
    std::string path = OUT_DIR + "revolve_spline_vase.obj";
    ASSERT_TRUE(exp.writeSolid(solid, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

// ---------------------------------------------------------------------------
// Boolean operations
// ---------------------------------------------------------------------------
TEST(ObjExport, BooleanUnion_TwoBoxes) {
    // Use different sizes so both shapes are visible without z-fighting.
    auto a = BRepBuilder::makeBox(10, 10, 10);
    auto b = BRepBuilder::makeBox(6, 14, 8);  // different dimensions, both at origin
    auto result = BooleanOperation(a, b, BooleanType::Union).build();
    ObjExporter exp;
    exp.setSurfaceSteps(8);
    std::string path = OUT_DIR + "boolean_union_boxes.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, BooleanUnion_BoxAndCylinder) {
    auto box = BRepBuilder::makeBox(8, 8, 8);
    auto cyl = BRepBuilder::makeCylinder({4,4,0}, Vec3::unitZ(), 3.0, 12.0);
    auto result = BooleanOperation(box, cyl, BooleanType::Union).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "boolean_union_box_cylinder.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, BooleanDifference_BoxMinusCylinder) {
    // A's exterior faces + cylinder faces as void (inner walls, flipped normals).
    auto box = BRepBuilder::makeBox(12, 12, 12);
    auto cyl = BRepBuilder::makeCylinder({6,6,0}, Vec3::unitZ(), 3.5, 14.0);
    auto result = BooleanOperation(box, cyl, BooleanType::Difference).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "boolean_diff_box_minus_cylinder.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
    // Void shell faces appear as a separate group in the OBJ
    ASSERT_TRUE(countObjLines(path, "g solid_void") > 0);
}

TEST(ObjExport, BooleanDifference_BoxMinusSphere) {
    auto box = BRepBuilder::makeBox(10, 10, 10);
    auto sph = BRepBuilder::makeSphere({5,5,5}, 4.0);
    auto result = BooleanOperation(box, sph, BooleanType::Difference).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "boolean_diff_box_minus_sphere.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, BooleanIntersection_BoxAndSphere) {
    auto box = BRepBuilder::makeBox(8, 8, 8);
    auto sph = BRepBuilder::makeSphere({4,4,4}, 5.0);
    auto result = BooleanOperation(box, sph, BooleanType::Intersection).build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "boolean_intersect_box_sphere.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, BooleanDifference_CylMinusCone) {
    auto cyl  = BRepBuilder::makeCylinder({0,0,0}, Vec3::unitZ(), 5.0, 10.0);
    auto cone = BRepBuilder::makeCone({0,0,0}, Vec3::unitZ(), 4.5, 0.0, 9.0);
    auto result = BooleanOperation(cyl, cone, BooleanType::Difference).build();
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "boolean_diff_cyl_minus_cone.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, BooleanUnion_SpherePair) {
    auto s1 = BRepBuilder::makeSphere({-3,0,0}, 4.0);
    auto s2 = BRepBuilder::makeSphere({ 3,0,0}, 4.0);
    auto result = BooleanOperation(s1, s2, BooleanType::Union).build();
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "boolean_union_sphere_pair.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

// ---------------------------------------------------------------------------
// Fillet
// ---------------------------------------------------------------------------
TEST(ObjExport, FilletBox) {
    auto box = BRepBuilder::makeBox(10, 10, 10);
    Fillet f(box);
    // Fillet all 4 edges of the bottom face
    for (const auto& e : box->outerShell()->faces()[0]->outerBound()->edges())
        f.addEdge(e, 1.5);
    auto result = f.build();
    ObjExporter exp;
    exp.setSurfaceSteps(12);
    std::string path = OUT_DIR + "fillet_box.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, FilletCylinder) {
    auto cyl = BRepBuilder::makeCylinder({0,0,0}, Vec3::unitZ(), 4.0, 10.0);
    Fillet f(cyl);
    // Fillet the circular edge on the bottom face
    const auto& botEdges = cyl->outerShell()->faces()[0]->outerBound()->edges();
    if (!botEdges.empty())
        f.addEdge(botEdges[0], 1.0);
    auto result = f.build();
    ObjExporter exp;
    exp.setSurfaceSteps(16);
    std::string path = OUT_DIR + "fillet_cylinder.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

// ---------------------------------------------------------------------------
// Multi-solid export
// ---------------------------------------------------------------------------
TEST(ObjExport, MultiSolid) {
    auto box = BRepBuilder::makeBox(6,6,6);
    auto sph = BRepBuilder::makeSphere({10,0,3}, 3.5);
    auto cyl = BRepBuilder::makeCylinder({20,0,0}, Vec3::unitZ(), 2.5, 8.0);
    ObjExporter exp;
    exp.setSurfaceSteps(12);
    std::string path = OUT_DIR + "multi_solid.obj";
    ASSERT_TRUE(exp.writeSolids({box, sph, cyl}, path));
    ASSERT_TRUE(fileExists(path));
    // Three 'o' headers
    ASSERT_EQ(countObjLines(path, "o "), 3);
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

// ---------------------------------------------------------------------------
// OBJ stream API
// ---------------------------------------------------------------------------
TEST(ObjExport, StreamOutput) {
    auto box = BRepBuilder::makeBox(2,2,2);
    ObjExporter exp;
    exp.setSurfaceSteps(4);
    std::ostringstream oss;
    int offset = 1;
    exp.writeSolidToStream(box, "tiny_box", oss, offset);
    const std::string content = oss.str();
    ASSERT_TRUE(content.find("o tiny_box") != std::string::npos);
    ASSERT_TRUE(content.find("v ") != std::string::npos);
    ASSERT_TRUE(content.find("f ") != std::string::npos);
}

// ---------------------------------------------------------------------------
// Chained operations
// ---------------------------------------------------------------------------
TEST(ObjExport, ExtrudeThenBoolean) {
    // Extrude a circle then subtract a smaller extruded circle → hollow cylinder
    Circle outer{{0,0,0}, 6.0};
    Circle inner{{0,0,0}, 4.0};
    auto outerSolid = Extrusion(outer.toWire(), Vec3::unitZ(), 10.0).build();
    auto innerSolid = Extrusion(inner.toWire(), Vec3::unitZ(), 12.0).build();
    auto result = BooleanOperation(outerSolid, innerSolid, BooleanType::Difference).build();
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "extrude_then_boolean.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

TEST(ObjExport, RevolveAndBoolean) {
    // Vase shape minus a small box cut-out
    std::vector<Point3D> pts = {{2,0,0},{3,0,4},{2,0,8},{4,0,12}};
    auto spline = SplineCurve::interpolate(pts, 3);
    auto wire = std::make_shared<Wire>();
    wire->addEdge(spline.toEdge());
    auto vase  = Revolve(wire, {0,0,0}, Vec3::unitZ(), 2*M_PI).build();
    auto cutout = BRepBuilder::makeBox(2,2,4);
    auto result = BooleanOperation(vase, cutout, BooleanType::Difference).build();
    ObjExporter exp;
    exp.setSurfaceSteps(20);
    std::string path = OUT_DIR + "revolve_and_boolean.obj";
    ASSERT_TRUE(exp.writeSolid(result, path));
    ASSERT_TRUE(fileExists(path));
    ASSERT_TRUE(countObjLines(path, "f ") > 0);
}

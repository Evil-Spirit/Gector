#include "test_framework.h"
#include "gector/sketch/circle.h"
#include "gector/sketch/line.h"
#include "gector/sketch/sketch.h"
#include "gector/brep/brep_builder.h"
#include "gector/operations/extrusion.h"
#include "gector/operations/revolve.h"
#include "gector/operations/boolean_op.h"
#include "gector/operations/fillet.h"
#include <cmath>

using namespace gector;

// ---------------------------------------------------------------------------
// Extrusion
// ---------------------------------------------------------------------------
TEST(Operations, ExtrudeCircle) {
    // Extrude a circle to produce a cylinder-like solid
    Circle c{{0,0,0}, 5.0};
    auto wire = c.toWire();
    auto solid = Extrusion(wire, Vec3::unitZ(), 10.0).build();
    ASSERT_TRUE(solid != nullptr);
    ASSERT_TRUE(solid->outerShell()->isClosed());
    // 1 bottom face + 1 top face + 1 side face per profile edge
    ASSERT_TRUE(solid->outerShell()->faceCount() >= 2);
}

TEST(Operations, ExtrudeRectangle) {
    // Extrude a rectangle (3 lines make open; let's use a 4-line sketch)
    Sketch sk;
    sk.addLine({0,0,0},{4,0,0});
    sk.addLine({4,0,0},{4,3,0});
    sk.addLine({4,3,0},{0,3,0});
    sk.addLine({0,3,0},{0,0,0});
    auto wire = sk.toWire();
    auto solid = Extrusion(wire, Vec3::unitZ(), 5.0).build();
    ASSERT_TRUE(solid != nullptr);
    // 2 caps + 4 side faces = 6 faces
    ASSERT_EQ(solid->outerShell()->faceCount(), std::size_t{6});
}

TEST(Operations, ExtrudeZeroDistanceThrows) {
    Circle c{{0,0,0}, 1.0};
    ASSERT_THROWS(Extrusion(c.toWire(), Vec3::unitZ(), 0.0));
}

TEST(Operations, ExtrudeNullWireThrows) {
    ASSERT_THROWS(Extrusion(nullptr, Vec3::unitZ(), 1.0));
}

TEST(Operations, ExtrudeByVector) {
    Circle c{{0,0,0}, 2.0};
    auto solid = Extrusion(c.toWire(), Vec3{0,0,7}).build();
    ASSERT_TRUE(solid != nullptr);
}

// ---------------------------------------------------------------------------
// Revolve
// ---------------------------------------------------------------------------
TEST(Operations, RevolveFullCircle) {
    // Revolve a vertical line at x=3 around the Z axis → cylinder
    Line l{{3,0,0},{3,0,10}};
    auto wire = std::make_shared<Wire>();
    wire->addEdge(l.toEdge());
    auto solid = Revolve(wire, {0,0,0}, Vec3::unitZ(), 2*M_PI).build();
    ASSERT_TRUE(solid != nullptr);
    ASSERT_TRUE(solid->outerShell()->isClosed());
}

TEST(Operations, RevolvePartialAngle) {
    Line l{{5,0,0},{5,0,3}};
    auto wire = std::make_shared<Wire>();
    wire->addEdge(l.toEdge());
    auto solid = Revolve(wire, {0,0,0}, Vec3::unitZ(), M_PI).build();
    ASSERT_TRUE(solid != nullptr);
    // Partial revolution produces start cap + end cap + surface
    ASSERT_TRUE(solid->outerShell()->faceCount() >= 3);
}

TEST(Operations, RevolveZeroAngleThrows) {
    Line l{{3,0,0},{3,0,5}};
    auto wire = std::make_shared<Wire>();
    wire->addEdge(l.toEdge());
    ASSERT_THROWS(Revolve(wire, {0,0,0}, Vec3::unitZ(), 0.0));
}

TEST(Operations, RevolveNullWireThrows) {
    ASSERT_THROWS(Revolve(nullptr, {0,0,0}, Vec3::unitZ(), M_PI));
}

// ---------------------------------------------------------------------------
// Boolean
// ---------------------------------------------------------------------------
TEST(Operations, BooleanUnion) {
    auto a = BRepBuilder::makeBox(5,5,5);
    auto b = BRepBuilder::makeBox(3,3,3);
    auto result = BooleanOperation(a, b, BooleanType::Union).build();
    ASSERT_TRUE(result != nullptr);
    ASSERT_TRUE(result->outerShell()->isClosed());
}

TEST(Operations, BooleanIntersection) {
    auto a = BRepBuilder::makeBox(5,5,5);
    auto b = BRepBuilder::makeBox(3,3,3);
    auto result = BooleanOperation(a, b, BooleanType::Intersection).build();
    ASSERT_TRUE(result != nullptr);
}

TEST(Operations, BooleanDifference) {
    auto a = BRepBuilder::makeBox(10,10,10);
    auto b = BRepBuilder::makeBox(5,5,5);
    auto result = BooleanOperation(a, b, BooleanType::Difference).build();
    ASSERT_TRUE(result != nullptr);
}

TEST(Operations, BooleanNullOperandThrows) {
    auto a = BRepBuilder::makeBox(1,1,1);
    ASSERT_THROWS(BooleanOperation(a, nullptr, BooleanType::Union));
    ASSERT_THROWS(BooleanOperation(nullptr, a, BooleanType::Union));
}

TEST(Operations, BooleanPreservesType) {
    auto a = BRepBuilder::makeBox(1,1,1);
    auto b = BRepBuilder::makeBox(1,1,1);
    BooleanOperation op(a, b, BooleanType::Difference);
    ASSERT_TRUE(op.type() == BooleanType::Difference);
    ASSERT_TRUE(op.operandA() == a);
    ASSERT_TRUE(op.operandB() == b);
}

// ---------------------------------------------------------------------------
// Fillet
// ---------------------------------------------------------------------------
TEST(Operations, FilletAddEdge) {
    auto box = BRepBuilder::makeBox(5,5,5);
    Fillet f(box);
    auto edge = box->outerShell()->faces()[0]->outerBound()->edges()[0];
    f.addEdge(edge, 1.0);
    ASSERT_EQ(f.specs().size(), std::size_t{1});
}

TEST(Operations, FilletBuild) {
    auto box = BRepBuilder::makeBox(5,5,5);
    Fillet f(box);
    auto edge = box->outerShell()->faces()[0]->outerBound()->edges()[0];
    f.addEdge(edge, 0.5);
    auto result = f.build();
    ASSERT_TRUE(result != nullptr);
    // Result has at least as many faces as original + one blend
    ASSERT_TRUE(result->outerShell()->faceCount() >= box->outerShell()->faceCount());
}

TEST(Operations, FilletNullSolidThrows) {
    ASSERT_THROWS(Fillet(nullptr));
}

TEST(Operations, FilletZeroRadiusThrows) {
    auto box = BRepBuilder::makeBox(1,1,1);
    Fillet f(box);
    auto edge = box->outerShell()->faces()[0]->outerBound()->edges()[0];
    ASSERT_THROWS(f.addEdge(edge, 0.0));
}

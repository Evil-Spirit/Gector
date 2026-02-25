#include "test_framework.h"
#include "gector/sketch/line.h"
#include "gector/sketch/arc.h"
#include "gector/sketch/circle.h"
#include "gector/sketch/spline_curve.h"
#include "gector/sketch/sketch.h"
#include <cmath>

using namespace gector;

// ---------------------------------------------------------------------------
// Line
// ---------------------------------------------------------------------------
TEST(Sketch, LineLength) {
    Line l{{0,0,0}, {3,4,0}};
    ASSERT_NEAR(l.length(), 5.0, 1e-10);
}

TEST(Sketch, LineDirection) {
    Line l{{0,0,0}, {1,0,0}};
    Vec3 d = l.direction();
    ASSERT_NEAR(d.x, 1.0, 1e-10);
}

TEST(Sketch, LineZeroLengthThrows) {
    ASSERT_THROWS(Line({1,2,3},{1,2,3}));
}

TEST(Sketch, LineToCurve) {
    Line l{{0,0,0},{10,0,0}};
    auto curve = l.toCurve();
    ASSERT_TRUE(curve != nullptr);
    ASSERT_NEAR(curve->evaluate(0.0).x, 0.0,  1e-8);
    ASSERT_NEAR(curve->evaluate(1.0).x, 10.0, 1e-8);
}

// ---------------------------------------------------------------------------
// Arc
// ---------------------------------------------------------------------------
TEST(Sketch, ArcLength) {
    // Quarter arc, radius 1 → length = π/2
    Arc a{{0,0,0}, 1.0, Vec3::unitZ(), 0.0, M_PI/2};
    ASSERT_NEAR(a.length(), M_PI/2, 1e-10);
}

TEST(Sketch, ArcStartEndPoints) {
    Arc a{{0,0,0}, 2.0, Vec3::unitZ(), 0.0, M_PI/2};
    ASSERT_NEAR(a.startPoint().x, 2.0, 1e-8);
    ASSERT_NEAR(a.startPoint().y, 0.0, 1e-8);
    ASSERT_NEAR(a.endPoint().x,   0.0, 1e-8);
    ASSERT_NEAR(a.endPoint().y,   2.0, 1e-8);
}

TEST(Sketch, ArcInvalidRadius) {
    ASSERT_THROWS(Arc({0,0,0}, -1.0, Vec3::unitZ(), 0, M_PI));
}

TEST(Sketch, ArcToCurveOnCircle) {
    Arc a{{0,0,0}, 3.0, Vec3::unitZ(), 0.0, M_PI};
    auto c = a.toCurve();
    ASSERT_TRUE(c != nullptr);
    // Mid-point should be at ~(0,3,0)
    Point3D mid = c->evaluate(0.5);
    double r = std::sqrt(mid.x*mid.x + mid.y*mid.y);
    ASSERT_NEAR(r, 3.0, 1e-6);
}

// ---------------------------------------------------------------------------
// Circle
// ---------------------------------------------------------------------------
TEST(Sketch, CircleLength) {
    Circle c{{0,0,0}, 5.0};
    ASSERT_NEAR(c.length(), 2*M_PI*5.0, 1e-8);
}

TEST(Sketch, CircleStartEqualsEnd) {
    Circle c{{0,0,0}, 1.0};
    ASSERT_NEAR(c.startPoint().distanceTo(c.endPoint()), 0.0, 1e-10);
}

TEST(Sketch, CircleInvalidRadius) {
    ASSERT_THROWS(Circle({0,0,0}, 0.0));
}

TEST(Sketch, CircleToWire) {
    Circle c{{0,0,0}, 2.0};
    auto w = c.toWire();
    ASSERT_TRUE(w != nullptr);
    ASSERT_EQ(w->edgeCount(), std::size_t{1});
    ASSERT_TRUE(w->isClosed());
}

// ---------------------------------------------------------------------------
// SplineCurve
// ---------------------------------------------------------------------------
TEST(Sketch, SplineInterpolate) {
    std::vector<Point3D> pts = {{0,0,0},{1,1,0},{2,0,0},{3,1,0}};
    auto s = SplineCurve::interpolate(pts, 3);
    ASSERT_NEAR(s.startPoint().x, 0.0, 1e-7);
    ASSERT_NEAR(s.endPoint().x,   3.0, 1e-7);
}

TEST(Sketch, SplineLength) {
    // Straight spline should have length ≈ straight-line distance
    std::vector<Point3D> pts = {{0,0,0},{5,0,0}};
    auto s = SplineCurve::interpolate(pts, 1);
    ASSERT_NEAR(s.length(), 5.0, 1e-4);
}

// ---------------------------------------------------------------------------
// Sketch
// ---------------------------------------------------------------------------
TEST(Sketch, SketchAddEntities) {
    Sketch sk;
    sk.addLine({0,0,0},{1,0,0});
    sk.addLine({1,0,0},{1,1,0});
    sk.addLine({1,1,0},{0,0,0});
    ASSERT_EQ(sk.entities().size(), std::size_t{3});
}

TEST(Sketch, SketchToWire) {
    Sketch sk;
    sk.addLine({0,0,0},{1,0,0});
    sk.addLine({1,0,0},{1,1,0});
    sk.addLine({1,1,0},{0,0,0});
    auto w = sk.toWire();
    ASSERT_EQ(w->edgeCount(), std::size_t{3});
}

TEST(Sketch, SketchToClosedWire) {
    Sketch sk;
    sk.addLine({0,0,0},{1,0,0});
    sk.addLine({1,0,0},{0,0,0});
    auto w = sk.toClosedWire();
    ASSERT_TRUE(w->isClosed());
}

TEST(Sketch, SketchCircle) {
    Sketch sk;
    sk.addCircle({0,0,0}, 5.0);
    ASSERT_EQ(sk.entities().size(), std::size_t{1});
    auto w = sk.toWire();
    ASSERT_EQ(w->edgeCount(), std::size_t{1});
}

TEST(Sketch, SketchNullEntityThrows) {
    Sketch sk;
    ASSERT_THROWS(sk.addEntity(nullptr));
}

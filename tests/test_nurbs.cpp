#include "test_framework.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"
#include <cmath>

using namespace gector;

// ---------------------------------------------------------------------------
// NURBSCurve: Line
// ---------------------------------------------------------------------------
TEST(NURBS, LineEvaluateEndpoints) {
    auto c = NURBSCurve::makeLine({0,0,0}, {10,0,0});
    Point3D p0 = c.evaluate(c.paramStart());
    Point3D p1 = c.evaluate(c.paramEnd());
    ASSERT_NEAR(p0.x, 0.0,  1e-10);
    ASSERT_NEAR(p1.x, 10.0, 1e-10);
}

TEST(NURBS, LineMidpoint) {
    auto c = NURBSCurve::makeLine({0,0,0}, {10,0,0});
    Point3D mid = c.evaluate(0.5);
    ASSERT_NEAR(mid.x, 5.0, 1e-10);
    ASSERT_NEAR(mid.y, 0.0, 1e-10);
}

TEST(NURBS, LineDerivative) {
    auto c = NURBSCurve::makeLine({0,0,0}, {10,0,0});
    Vec3 d = c.derivative(0.5);
    // Tangent should point in +X direction (un-normalised magnitude = 10 for parametric speed)
    ASSERT_TRUE(d.x > 0);
    ASSERT_NEAR(d.y, 0.0, 1e-10);
}

// ---------------------------------------------------------------------------
// NURBSCurve: Arc
// ---------------------------------------------------------------------------
TEST(NURBS, ArcQuarterEndpoints) {
    // Quarter arc in XY plane: 0 → π/2
    auto c = NURBSCurve::makeArc({0,0,0}, 1.0,
                                   Vec3::unitX(), Vec3::unitY(),
                                   0.0, M_PI / 2.0);
    Point3D p0 = c.evaluate(c.paramStart());
    Point3D p1 = c.evaluate(c.paramEnd());
    ASSERT_NEAR(p0.x, 1.0, 1e-8);
    ASSERT_NEAR(p0.y, 0.0, 1e-8);
    ASSERT_NEAR(p1.x, 0.0, 1e-8);
    ASSERT_NEAR(p1.y, 1.0, 1e-8);
}

TEST(NURBS, ArcMidpointOnCircle) {
    // Point at π/4 should lie on the unit circle
    auto c = NURBSCurve::makeArc({0,0,0}, 1.0,
                                   Vec3::unitX(), Vec3::unitY(),
                                   0.0, M_PI / 2.0);
    Point3D mid = c.evaluate(0.5);
    double r = std::sqrt(mid.x*mid.x + mid.y*mid.y);
    ASSERT_NEAR(r, 1.0, 1e-8);
}

TEST(NURBS, ArcSemicircle) {
    // Semicircle: 0 → π
    auto c = NURBSCurve::makeArc({0,0,0}, 5.0,
                                   Vec3::unitX(), Vec3::unitY(),
                                   0.0, M_PI);
    Point3D p0 = c.evaluate(c.paramStart());
    Point3D p1 = c.evaluate(c.paramEnd());
    ASSERT_NEAR(p0.x,  5.0, 1e-7);
    ASSERT_NEAR(p1.x, -5.0, 1e-7);
    // Midpoint should be at (0,5,0)
    Point3D mid = c.evaluate(0.5);
    ASSERT_NEAR(mid.x, 0.0, 1e-7);
    ASSERT_NEAR(mid.y, 5.0, 1e-7);
}

// ---------------------------------------------------------------------------
// NURBSCurve: Circle
// ---------------------------------------------------------------------------
TEST(NURBS, CircleRadius) {
    auto c = NURBSCurve::makeCircle({1,2,0}, 3.0, Vec3::unitZ());
    // Sample several points; all should be distance 3 from centre
    for (int i = 0; i <= 8; ++i) {
        double t = c.paramStart() +
                   (c.paramEnd() - c.paramStart()) * i / 8.0;
        // Avoid exact endpoint at full 2π boundary
        if (i == 8) t -= 1e-9;
        Point3D p = c.evaluate(t);
        double r = std::sqrt((p.x-1)*(p.x-1) + (p.y-2)*(p.y-2));
        ASSERT_NEAR(r, 3.0, 1e-6);
    }
}

// ---------------------------------------------------------------------------
// NURBSCurve: Interpolated spline
// ---------------------------------------------------------------------------
TEST(NURBS, InterpolatedPassesThrough) {
    std::vector<Point3D> pts = {{0,0,0}, {1,1,0}, {2,0,0}, {3,1,0}};
    auto c = NURBSCurve::makeInterpolated(pts, 3);
    // Endpoints must match exactly
    ASSERT_NEAR(c.evaluate(c.paramStart()).x, 0.0, 1e-7);
    ASSERT_NEAR(c.evaluate(c.paramEnd()).x,   3.0, 1e-7);
}

// ---------------------------------------------------------------------------
// NURBSCurve: Split
// ---------------------------------------------------------------------------
TEST(NURBS, SplitContinuity) {
    auto c = NURBSCurve::makeLine({0,0,0}, {10,0,0});
    auto [left, right] = c.split(0.5);
    Point3D leftEnd  = left.evaluate(left.paramEnd());
    Point3D rightStart = right.evaluate(right.paramStart());
    ASSERT_NEAR(leftEnd.x,    5.0, 1e-8);
    ASSERT_NEAR(rightStart.x, 5.0, 1e-8);
}

// ---------------------------------------------------------------------------
// NURBSSurface: Plane
// ---------------------------------------------------------------------------
TEST(NURBS, PlaneCorners) {
    auto s = NURBSSurface::makePlane({0,0,0}, Vec3::unitX(), Vec3::unitY(), 1.0, 1.0);
    Point3D p00 = s.evaluate(0,0);
    Point3D p10 = s.evaluate(1,0);
    Point3D p01 = s.evaluate(0,1);
    Point3D p11 = s.evaluate(1,1);
    ASSERT_NEAR(p00.x, 0.0, 1e-10); ASSERT_NEAR(p00.y, 0.0, 1e-10);
    ASSERT_NEAR(p10.x, 1.0, 1e-10); ASSERT_NEAR(p10.y, 0.0, 1e-10);
    ASSERT_NEAR(p01.x, 0.0, 1e-10); ASSERT_NEAR(p01.y, 1.0, 1e-10);
    ASSERT_NEAR(p11.x, 1.0, 1e-10); ASSERT_NEAR(p11.y, 1.0, 1e-10);
}

TEST(NURBS, PlaneNormal) {
    auto s = NURBSSurface::makePlane({0,0,0}, Vec3::unitX(), Vec3::unitY(), 2.0, 2.0);
    Vec3 n = s.normal(0.5, 0.5);
    ASSERT_NEAR(std::abs(n.z), 1.0, 1e-8);
}

// ---------------------------------------------------------------------------
// NURBSSurface: Cylinder
// ---------------------------------------------------------------------------
TEST(NURBS, CylinderRadius) {
    auto s = NURBSSurface::makeCylinder({0,0,0}, Vec3::unitZ(), 4.0, 10.0);
    // Sample several u values at v=0 (bottom ring)
    for (int i = 0; i <= 8; ++i) {
        double u = s.uParamStart() +
                   (s.uParamEnd() - s.uParamStart()) * i / 8.0;
        if (i == 8) u = s.uParamEnd() - 1e-9;
        Point3D p = s.evaluate(u, 0.0);
        double r = std::sqrt(p.x*p.x + p.y*p.y);
        ASSERT_NEAR(r, 4.0, 1e-6);
        ASSERT_NEAR(p.z, 0.0, 1e-6);
    }
}

// ---------------------------------------------------------------------------
// NURBSSurface: Revolution surface
// ---------------------------------------------------------------------------
TEST(NURBS, RevolutionSurfaceRadius) {
    // Profile: vertical line at x=3 from z=0 to z=5
    auto profile = NURBSCurve::makeLine({3,0,0}, {3,0,5});
    auto s = NURBSSurface::makeRevolutionSurface(
        profile, {0,0,0}, Vec3::unitZ(), 2*M_PI);
    // All points should be at radius 3
    for (int i = 0; i <= 4; ++i) {
        double u = s.uParamStart() +
                   (s.uParamEnd() - s.uParamStart()) * i / 4.0;
        if (i == 4) u = s.uParamEnd() - 1e-9;
        Point3D p = s.evaluate(u, 0.5);
        double r = std::sqrt(p.x*p.x + p.y*p.y);
        ASSERT_NEAR(r, 3.0, 1e-5);
    }
}

// ---------------------------------------------------------------------------
// NURBSSurface: Sphere (approximate test)
// ---------------------------------------------------------------------------
TEST(NURBS, SphereRadius) {
    auto s = NURBSSurface::makeSphere({0,0,0}, 2.0);
    // Interior sample
    Point3D c = s.evaluate(0.5, 0.5);
    double r = c.length();
    ASSERT_NEAR(r, 2.0, 1e-5);
}

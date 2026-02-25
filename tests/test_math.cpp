#include "test_framework.h"
#include "gector/math/vec3.h"
#include "gector/math/mat4.h"
#include <cmath>

using namespace gector;

// ---------------------------------------------------------------------------
// Vec3 tests
// ---------------------------------------------------------------------------
TEST(Math, Vec3DefaultConstruct) {
    Vec3 v;
    ASSERT_NEAR(v.x, 0.0, 1e-15);
    ASSERT_NEAR(v.y, 0.0, 1e-15);
    ASSERT_NEAR(v.z, 0.0, 1e-15);
}

TEST(Math, Vec3ArithmeticAdd) {
    Vec3 a{1,2,3}, b{4,5,6};
    Vec3 c = a + b;
    ASSERT_NEAR(c.x, 5.0, 1e-15);
    ASSERT_NEAR(c.y, 7.0, 1e-15);
    ASSERT_NEAR(c.z, 9.0, 1e-15);
}

TEST(Math, Vec3Scale) {
    Vec3 v{1,2,3};
    Vec3 s = v * 2.0;
    ASSERT_NEAR(s.x, 2.0, 1e-15);
    ASSERT_NEAR(s.y, 4.0, 1e-15);
    ASSERT_NEAR(s.z, 6.0, 1e-15);
}

TEST(Math, Vec3Dot) {
    Vec3 a{1,0,0}, b{0,1,0};
    ASSERT_NEAR(a.dot(b), 0.0, 1e-15);
    ASSERT_NEAR(a.dot(a), 1.0, 1e-15);
}

TEST(Math, Vec3Cross) {
    Vec3 x = Vec3::unitX();
    Vec3 y = Vec3::unitY();
    Vec3 z = x.cross(y);
    ASSERT_NEAR(z.x, 0.0, 1e-15);
    ASSERT_NEAR(z.y, 0.0, 1e-15);
    ASSERT_NEAR(z.z, 1.0, 1e-15);
}

TEST(Math, Vec3Length) {
    Vec3 v{3,4,0};
    ASSERT_NEAR(v.length(), 5.0, 1e-12);
}

TEST(Math, Vec3Normalized) {
    Vec3 v{3,4,0};
    Vec3 n = v.normalized();
    ASSERT_NEAR(n.length(), 1.0, 1e-12);
    ASSERT_NEAR(n.x, 0.6, 1e-12);
    ASSERT_NEAR(n.y, 0.8, 1e-12);
}

TEST(Math, Vec3NormalizeZeroThrows) {
    ASSERT_THROWS(Vec3{}.normalized());
}

TEST(Math, Vec3Lerp) {
    Vec3 a{0,0,0}, b{10,0,0};
    Vec3 m = a.lerp(b, 0.5);
    ASSERT_NEAR(m.x, 5.0, 1e-12);
}

TEST(Math, Vec3Distance) {
    Vec3 a{0,0,0}, b{3,4,0};
    ASSERT_NEAR(a.distanceTo(b), 5.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Mat4 tests
// ---------------------------------------------------------------------------
TEST(Math, Mat4Identity) {
    Mat4 m = Mat4::identity();
    Point3D p{1,2,3};
    Point3D t = m.transformPoint(p);
    ASSERT_NEAR(t.x, 1.0, 1e-12);
    ASSERT_NEAR(t.y, 2.0, 1e-12);
    ASSERT_NEAR(t.z, 3.0, 1e-12);
}

TEST(Math, Mat4Translation) {
    Mat4 m = Mat4::translation({5,6,7});
    Point3D t = m.transformPoint({1,2,3});
    ASSERT_NEAR(t.x, 6.0, 1e-12);
    ASSERT_NEAR(t.y, 8.0, 1e-12);
    ASSERT_NEAR(t.z, 10.0, 1e-12);
}

TEST(Math, Mat4Scaling) {
    Mat4 m = Mat4::scale(2,3,4);
    Point3D t = m.transformPoint({1,1,1});
    ASSERT_NEAR(t.x, 2.0, 1e-12);
    ASSERT_NEAR(t.y, 3.0, 1e-12);
    ASSERT_NEAR(t.z, 4.0, 1e-12);
}

TEST(Math, Mat4RotationZ90) {
    // Rotate (1,0,0) by 90° around Z → (0,1,0)
    Mat4 m = Mat4::rotation(Vec3::unitZ(), M_PI / 2.0);
    Vec3 v = m.transformVector({1,0,0});
    ASSERT_NEAR(v.x, 0.0, 1e-12);
    ASSERT_NEAR(v.y, 1.0, 1e-12);
    ASSERT_NEAR(v.z, 0.0, 1e-12);
}

TEST(Math, Mat4Multiply) {
    Mat4 t = Mat4::translation({1,0,0});
    Mat4 s = Mat4::scale(2,2,2);
    Mat4 ts = t * s; // scale first, then translate
    Point3D p = ts.transformPoint({1,0,0});
    // s(1,0,0) = (2,0,0), t(2,0,0) = (3,0,0)
    ASSERT_NEAR(p.x, 3.0, 1e-12);
}

TEST(Math, Mat4Inverse) {
    Mat4 t = Mat4::translation({3,4,5});
    Mat4 inv = t.inverse();
    Point3D p = inv.transformPoint({3,4,5});
    ASSERT_NEAR(p.x, 0.0, 1e-10);
    ASSERT_NEAR(p.y, 0.0, 1e-10);
    ASSERT_NEAR(p.z, 0.0, 1e-10);
}

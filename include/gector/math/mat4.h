#pragma once
#include "vec3.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace gector {

/// 4x4 column-major double matrix used for homogeneous transforms.
/// Stored as m[col][row] so that matrix * column-vector is natural.
struct Mat4 {
    // m[col][row]
    double m[4][4];

    /// Initialise to identity by default.
    Mat4() noexcept {
        for (int c = 0; c < 4; ++c)
            for (int r = 0; r < 4; ++r)
                m[c][r] = (c == r) ? 1.0 : 0.0;
    }

    /// Access element at row r, column c.
    double  at(int r, int c) const { return m[c][r]; }
    double& at(int r, int c)       { return m[c][r]; }

    // ---------------------------------------------------------------
    // Algebra
    // ---------------------------------------------------------------
    Mat4 operator*(const Mat4& b) const noexcept {
        Mat4 res;
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c) {
                double s = 0;
                for (int k = 0; k < 4; ++k)
                    s += at(r,k) * b.at(k,c);
                res.at(r,c) = s;
            }
        return res;
    }

    /// Transform a homogeneous 4-vector.
    std::array<double,4> operator*(const std::array<double,4>& v) const noexcept {
        std::array<double,4> res{0,0,0,0};
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                res[r] += at(r,c) * v[c];
        return res;
    }

    /// Transform a Point3D (w=1).
    Point3D transformPoint(const Point3D& p) const noexcept {
        double x = at(0,0)*p.x + at(0,1)*p.y + at(0,2)*p.z + at(0,3);
        double y = at(1,0)*p.x + at(1,1)*p.y + at(1,2)*p.z + at(1,3);
        double z = at(2,0)*p.x + at(2,1)*p.y + at(2,2)*p.z + at(2,3);
        double w = at(3,0)*p.x + at(3,1)*p.y + at(3,2)*p.z + at(3,3);
        return {x/w, y/w, z/w};
    }

    /// Transform a direction vector (w=0, no translation).
    Vec3 transformVector(const Vec3& v) const noexcept {
        return {
            at(0,0)*v.x + at(0,1)*v.y + at(0,2)*v.z,
            at(1,0)*v.x + at(1,1)*v.y + at(1,2)*v.z,
            at(2,0)*v.x + at(2,1)*v.y + at(2,2)*v.z
        };
    }

    Mat4 transposed() const noexcept {
        Mat4 t;
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                t.at(r,c) = at(c,r);
        return t;
    }

    /// Inverse for a general 4×4 matrix (Gauss-Jordan).
    Mat4 inverse() const {
        double a[4][8];
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) a[r][c]   = at(r,c);
            for (int c = 0; c < 4; ++c) a[r][c+4] = (r==c) ? 1.0 : 0.0;
        }
        for (int col = 0; col < 4; ++col) {
            // Pivot
            int piv = col;
            for (int r = col+1; r < 4; ++r)
                if (std::abs(a[r][col]) > std::abs(a[piv][col])) piv = r;
            if (std::abs(a[piv][col]) < 1e-15)
                throw std::runtime_error("Matrix is singular");
            std::swap(a[col], a[piv]);
            double s = 1.0 / a[col][col];
            for (int c = 0; c < 8; ++c) a[col][c] *= s;
            for (int r = 0; r < 4; ++r) {
                if (r == col) continue;
                double f = a[r][col];
                for (int c = 0; c < 8; ++c) a[r][c] -= f * a[col][c];
            }
        }
        Mat4 inv;
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                inv.at(r,c) = a[r][c+4];
        return inv;
    }

    // ---------------------------------------------------------------
    // Static factories
    // ---------------------------------------------------------------
    static Mat4 identity() noexcept { return {}; }

    static Mat4 translation(const Vec3& t) noexcept {
        Mat4 m;
        m.at(0,3) = t.x;
        m.at(1,3) = t.y;
        m.at(2,3) = t.z;
        return m;
    }

    static Mat4 scale(double sx, double sy, double sz) noexcept {
        Mat4 m;
        m.at(0,0) = sx;
        m.at(1,1) = sy;
        m.at(2,2) = sz;
        return m;
    }
    static Mat4 scale(double s) noexcept { return scale(s,s,s); }

    /// Rotation by `angle` radians around unit vector `axis`.
    static Mat4 rotation(const Vec3& axis, double angle) noexcept {
        double c = std::cos(angle);
        double s = std::sin(angle);
        double t = 1.0 - c;
        const Vec3& a = axis; // assumed unit
        Mat4 m;
        m.at(0,0) = t*a.x*a.x + c;
        m.at(0,1) = t*a.x*a.y - s*a.z;
        m.at(0,2) = t*a.x*a.z + s*a.y;
        m.at(1,0) = t*a.x*a.y + s*a.z;
        m.at(1,1) = t*a.y*a.y + c;
        m.at(1,2) = t*a.y*a.z - s*a.x;
        m.at(2,0) = t*a.x*a.z - s*a.y;
        m.at(2,1) = t*a.y*a.z + s*a.x;
        m.at(2,2) = t*a.z*a.z + c;
        return m;
    }

    /// Build a rotation that maps fromDir → toDir (both should be unit vectors).
    static Mat4 rotationFromTo(const Vec3& from, const Vec3& to) {
        Vec3 axis = from.cross(to);
        if (axis.isZero()) return identity();
        double angle = std::acos(std::max(-1.0, std::min(1.0, from.dot(to))));
        return rotation(axis.normalized(), angle);
    }
};

} // namespace gector

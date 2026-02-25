#pragma once
#include <cmath>
#include <stdexcept>
#include <ostream>

// Portable math constant: M_PI is not part of the C++ standard.
// _USE_MATH_DEFINES (set in CMakeLists.txt) exposes it on MSVC; the guard
// below provides a fallback for any other build configuration.
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace gector {

/// 3-component double-precision vector / point.
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vec3() noexcept = default;
    constexpr Vec3(double x, double y, double z) noexcept : x(x), y(y), z(z) {}

    // Arithmetic
    constexpr Vec3  operator+ (const Vec3& o) const noexcept { return {x+o.x, y+o.y, z+o.z}; }
    constexpr Vec3  operator- (const Vec3& o) const noexcept { return {x-o.x, y-o.y, z-o.z}; }
    constexpr Vec3  operator- ()              const noexcept { return {-x, -y, -z}; }
    constexpr Vec3  operator* (double s)      const noexcept { return {x*s, y*s, z*s}; }
    constexpr Vec3  operator/ (double s)      const noexcept { return {x/s, y/s, z/s}; }
    constexpr Vec3& operator+=(const Vec3& o)       noexcept { x+=o.x; y+=o.y; z+=o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o)       noexcept { x-=o.x; y-=o.y; z-=o.z; return *this; }
    constexpr Vec3& operator*=(double s)             noexcept { x*=s;   y*=s;   z*=s;   return *this; }
    constexpr Vec3& operator/=(double s)             noexcept { x/=s;   y/=s;   z/=s;   return *this; }

    // Comparison
    constexpr bool operator==(const Vec3& o) const noexcept { return x==o.x && y==o.y && z==o.z; }
    constexpr bool operator!=(const Vec3& o) const noexcept { return !(*this==o); }

    // Vector operations
    constexpr double dot(const Vec3& o) const noexcept {
        return x*o.x + y*o.y + z*o.z;
    }
    constexpr Vec3 cross(const Vec3& o) const noexcept {
        return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
    }
    double length()    const noexcept { return std::sqrt(x*x + y*y + z*z); }
    constexpr double lengthSq() const noexcept { return x*x + y*y + z*z; }

    Vec3 normalized() const {
        const double len = length();
        if (len < 1e-15)
            throw std::runtime_error("Cannot normalize a zero-length vector");
        return *this / len;
    }
    bool isZero(double tol = 1e-10) const noexcept { return lengthSq() < tol*tol; }

    /// Linear interpolation: (1-t)*this + t*other
    constexpr Vec3 lerp(const Vec3& other, double t) const noexcept {
        return *this * (1.0 - t) + other * t;
    }
    double distanceTo(const Vec3& o) const noexcept { return (*this - o).length(); }

    // Index access
    double  operator[](int i) const { return (&x)[i]; }
    double& operator[](int i)       { return (&x)[i]; }

    // Static factories
    static constexpr Vec3 zero()  noexcept { return {0,0,0}; }
    static constexpr Vec3 unitX() noexcept { return {1,0,0}; }
    static constexpr Vec3 unitY() noexcept { return {0,1,0}; }
    static constexpr Vec3 unitZ() noexcept { return {0,0,1}; }
};

constexpr inline Vec3 operator*(double s, const Vec3& v) noexcept { return v * s; }

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

using Point3D = Vec3;

/// Build a unit vector that lies in the plane with the given unit normal,
/// preferring alignment with the global X axis.
inline Vec3 buildXAxisFromNormal(const Vec3& normal) noexcept {
    // Project the global X (or Y when nearly parallel) onto the plane
    Vec3 ref = (std::abs(normal.x) < 0.9) ? Vec3{1,0,0} : Vec3{0,1,0};
    Vec3 x = ref - normal * normal.dot(ref);
    double len = x.length();
    return (len < 1e-15) ? Vec3{0,1,0} : x / len;
}

} // namespace gector

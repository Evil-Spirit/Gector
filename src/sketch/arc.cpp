#include "gector/sketch/arc.h"
#include <stdexcept>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gector {

Arc::Arc(const Point3D& center, double radius, const Vec3& normal,
         double startAngle, double endAngle)
    : m_center(center)
    , m_radius(radius)
    , m_normal(normal.normalized())
    , m_startAngle(startAngle)
    , m_endAngle(endAngle)
{
    if (radius <= 0)
        throw std::invalid_argument("Arc radius must be positive");
    if (endAngle <= startAngle)
        throw std::invalid_argument("Arc endAngle must be greater than startAngle");
    if (endAngle - startAngle > 2.0 * M_PI + 1e-10)
        throw std::invalid_argument("Arc sweep must not exceed 2π");
}

Vec3 Arc::xAxis() const {
    return buildXAxisFromNormal(m_normal);
}

Vec3 Arc::yAxis() const {
    return m_normal.cross(xAxis()).normalized();
}

Point3D Arc::startPoint() const {
    Vec3 x = xAxis(), y = yAxis();
    return m_center
        + x * (m_radius * std::cos(m_startAngle))
        + y * (m_radius * std::sin(m_startAngle));
}

Point3D Arc::endPoint() const {
    Vec3 x = xAxis(), y = yAxis();
    return m_center
        + x * (m_radius * std::cos(m_endAngle))
        + y * (m_radius * std::sin(m_endAngle));
}

double Arc::length() const {
    return m_radius * (m_endAngle - m_startAngle);
}

std::shared_ptr<NURBSCurve> Arc::toCurve() const {
    return std::make_shared<NURBSCurve>(
        NURBSCurve::makeArc(m_center, m_radius, xAxis(), yAxis(),
                             m_startAngle, m_endAngle));
}

} // namespace gector

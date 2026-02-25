#include "gector/sketch/circle.h"
#include <cmath>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gector {

Circle::Circle(const Point3D& center, double radius, const Vec3& normal)
    : m_center(center), m_radius(radius), m_normal(normal.normalized())
{
    if (radius <= 0)
        throw std::invalid_argument("Circle radius must be positive");
}

Point3D Circle::startPoint() const {
    Vec3 x = buildXAxisFromNormal(m_normal);
    return m_center + x * m_radius;
}

Point3D Circle::endPoint() const {
    return startPoint(); // closed
}

double Circle::length() const {
    return 2.0 * M_PI * m_radius;
}

std::shared_ptr<NURBSCurve> Circle::toCurve() const {
    return std::make_shared<NURBSCurve>(
        NURBSCurve::makeCircle(m_center, m_radius, m_normal));
}

WirePtr Circle::toWire() const {
    auto curve = toCurve();
    auto v = std::make_shared<Vertex>(startPoint());
    auto w = std::make_shared<Wire>();
    w->addEdge(std::make_shared<Edge>(v, v, curve));
    return w;
}

} // namespace gector

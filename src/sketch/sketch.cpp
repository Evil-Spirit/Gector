#include "gector/sketch/sketch.h"
#include "gector/sketch/line.h"
#include "gector/sketch/arc.h"
#include "gector/sketch/circle.h"
#include "gector/sketch/spline_curve.h"
#include <stdexcept>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gector {

Sketch::Sketch(const Point3D& origin, const Vec3& normal)
    : m_origin(origin), m_normal(normal.normalized()) {}

Vec3 Sketch::xAxis() const {
    Vec3 n = m_normal;
    Vec3 x;
    if (std::abs(n.dot(Vec3::unitX())) < 0.9)
        x = n.cross(Vec3::unitX()).normalized();
    else
        x = n.cross(Vec3::unitY()).normalized();
    return x;
}

Vec3 Sketch::yAxis() const {
    return m_normal.cross(xAxis()).normalized();
}

void Sketch::addEntity(std::shared_ptr<SketchEntity> entity) {
    if (!entity)
        throw std::invalid_argument("Cannot add null entity to sketch");
    m_entities.push_back(std::move(entity));
}

void Sketch::addLine(const Point3D& p0, const Point3D& p1) {
    addEntity(std::make_shared<Line>(p0, p1));
}

void Sketch::addArc(const Point3D& center, double radius,
                    double startAngle, double endAngle) {
    addEntity(std::make_shared<Arc>(center, radius, m_normal, startAngle, endAngle));
}

void Sketch::addCircle(const Point3D& center, double radius) {
    addEntity(std::make_shared<Circle>(center, radius, m_normal));
}

void Sketch::addSpline(const std::vector<Point3D>& points, int degree) {
    addEntity(std::make_shared<SplineCurve>(SplineCurve::interpolate(points, degree)));
}

WirePtr Sketch::toWire() const {
    auto wire = std::make_shared<Wire>();
    for (const auto& entity : m_entities)
        wire->addEdge(entity->toEdge());
    return wire;
}

WirePtr Sketch::toClosedWire() const {
    auto wire = toWire();
    if (!wire->isClosed())
        throw std::runtime_error("Sketch entities do not form a closed wire");
    return wire;
}

} // namespace gector

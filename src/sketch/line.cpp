#include "gector/sketch/line.h"
#include <stdexcept>

namespace gector {

Line::Line(const Point3D& start, const Point3D& end)
    : m_start(start), m_end(end)
{
    if (start.distanceTo(end) < 1e-15)
        throw std::invalid_argument("Line start and end points must be distinct");
}

std::shared_ptr<NURBSCurve> Line::toCurve() const {
    return std::make_shared<NURBSCurve>(NURBSCurve::makeLine(m_start, m_end));
}

Vec3 Line::direction() const {
    return (m_end - m_start).normalized();
}

} // namespace gector

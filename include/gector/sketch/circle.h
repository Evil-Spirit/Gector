#pragma once
#include "sketch_entity.h"

namespace gector {

/**
 * @brief Full circle (closed arc from 0 to 2π).
 *
 * A circle's start and end points are identical.
 */
class Circle : public SketchEntity {
public:
    Circle(const Point3D& center, double radius,
           const Vec3& normal = Vec3::unitZ());

    const char* typeName() const noexcept override { return "Circle"; }

    Point3D startPoint() const override;
    Point3D endPoint()   const override;
    double  length()     const override;

    std::shared_ptr<NURBSCurve> toCurve() const override;

    const Point3D& center() const noexcept { return m_center; }
    double         radius() const noexcept { return m_radius; }
    const Vec3&    normal() const noexcept { return m_normal; }

    /// Closed wire wrapping the circle edge (convenient for operations).
    WirePtr toWire() const;

private:
    Point3D m_center;
    double  m_radius;
    Vec3    m_normal;
};

} // namespace gector

#pragma once
#include "sketch_entity.h"

namespace gector {

/**
 * @brief Circular arc defined by centre, radius, plane normal and angular range.
 *
 * Angles are measured from xAxis (the reference axis in the arc plane)
 * counter-clockwise when viewed along the normal direction.
 */
class Arc : public SketchEntity {
public:
    /**
     * @param center      Centre of the circle
     * @param radius      Radius (> 0)
     * @param normal      Unit normal to the arc plane
     * @param startAngle  Start angle in radians
     * @param endAngle    End angle in radians (> startAngle, ≤ startAngle + 2π)
     */
    Arc(const Point3D& center, double radius, const Vec3& normal,
        double startAngle, double endAngle);

    const char* typeName() const noexcept override { return "Arc"; }

    Point3D startPoint() const override;
    Point3D endPoint()   const override;
    double  length()     const override;

    std::shared_ptr<NURBSCurve> toCurve() const override;

    const Point3D& center()     const noexcept { return m_center; }
    double         radius()     const noexcept { return m_radius; }
    const Vec3&    normal()     const noexcept { return m_normal; }
    double         startAngle() const noexcept { return m_startAngle; }
    double         endAngle()   const noexcept { return m_endAngle; }

    /// Reference X axis in the arc plane.
    Vec3 xAxis() const;
    /// Reference Y axis in the arc plane.
    Vec3 yAxis() const;

private:
    Point3D m_center;
    double  m_radius;
    Vec3    m_normal;
    double  m_startAngle;
    double  m_endAngle;
};

} // namespace gector

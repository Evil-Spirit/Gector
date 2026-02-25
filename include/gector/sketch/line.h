#pragma once
#include "sketch_entity.h"

namespace gector {

/**
 * @brief Straight line segment between two 3-D points.
 */
class Line : public SketchEntity {
public:
    Line(const Point3D& start, const Point3D& end);

    const char* typeName() const noexcept override { return "Line"; }

    Point3D startPoint() const override { return m_start; }
    Point3D endPoint()   const override { return m_end; }
    double  length()     const override { return m_start.distanceTo(m_end); }

    std::shared_ptr<NURBSCurve> toCurve() const override;

    /// Direction (unit) vector from start to end.
    Vec3 direction() const;

private:
    Point3D m_start;
    Point3D m_end;
};

} // namespace gector

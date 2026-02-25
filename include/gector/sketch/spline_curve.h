#pragma once
#include "sketch_entity.h"
#include <vector>

namespace gector {

/**
 * @brief Free-form NURBS spline curve used in a sketch.
 *
 * Can be created by specifying control points + weights + knots directly,
 * or through the convenient interpolation factory that fits a curve through
 * a set of pass-through points.
 */
class SplineCurve : public SketchEntity {
public:
    /// Direct NURBS specification.
    explicit SplineCurve(NURBSCurve curve);

    /// Interpolating spline through the given points (degree 3 by default).
    static SplineCurve interpolate(const std::vector<Point3D>& points,
                                    int degree = 3);

    const char* typeName() const noexcept override { return "SplineCurve"; }

    Point3D startPoint() const override;
    Point3D endPoint()   const override;
    double  length()     const override;

    std::shared_ptr<NURBSCurve> toCurve() const override;

    const NURBSCurve& curve() const noexcept { return m_curve; }

private:
    NURBSCurve m_curve;
};

} // namespace gector

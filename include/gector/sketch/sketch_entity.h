#pragma once
#include "gector/math/vec3.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/brep/topology.h"
#include <memory>

namespace gector {

/**
 * @brief Abstract base for all sketch entities.
 *
 * A sketch entity lives in a plane and can produce:
 *   - A NURBS curve (geometric description),
 *   - A B-Rep edge (topological description).
 */
class SketchEntity {
public:
    virtual ~SketchEntity() = default;

    /// Unique entity type name for identification / serialisation.
    virtual const char* typeName() const noexcept = 0;

    /// 3-D start point.
    virtual Point3D startPoint() const = 0;

    /// 3-D end point.
    virtual Point3D endPoint() const = 0;

    /// Approximate arc length.
    virtual double length() const = 0;

    /// Export as a NURBS curve.
    virtual std::shared_ptr<NURBSCurve> toCurve() const = 0;

    /// Export as a B-Rep edge (creates start/end vertices automatically).
    virtual EdgePtr toEdge() const;
};

} // namespace gector

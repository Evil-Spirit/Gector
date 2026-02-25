#pragma once
#include "sketch_entity.h"
#include "gector/math/mat4.h"
#include <vector>
#include <memory>

namespace gector {

/**
 * @brief A 2-D sketch living in a plane defined by an origin and normal.
 *
 * A Sketch contains an ordered collection of SketchEntity objects.
 * It can produce a closed or open Wire for use in operations.
 */
class Sketch {
public:
    /**
     * @param origin  A point on the sketch plane.
     * @param normal  Unit normal to the sketch plane.
     */
    Sketch(const Point3D& origin = Point3D{},
           const Vec3&    normal = Vec3::unitZ());

    // -----------------------------------------------------------------------
    // Plane
    // -----------------------------------------------------------------------
    const Point3D& origin() const noexcept { return m_origin; }
    const Vec3&    normal() const noexcept { return m_normal; }
    /// Local X axis in the sketch plane.
    Vec3 xAxis() const;
    /// Local Y axis in the sketch plane.
    Vec3 yAxis() const;

    // -----------------------------------------------------------------------
    // Entity management
    // -----------------------------------------------------------------------
    void addEntity(std::shared_ptr<SketchEntity> entity);

    /// Convenience: add a line between two 3-D points.
    void addLine(const Point3D& p0, const Point3D& p1);

    /// Convenience: add a circular arc.
    void addArc(const Point3D& center, double radius,
                double startAngle, double endAngle);

    /// Convenience: add a full circle.
    void addCircle(const Point3D& center, double radius);

    /// Convenience: add an interpolating spline.
    void addSpline(const std::vector<Point3D>& points, int degree = 3);

    const std::vector<std::shared_ptr<SketchEntity>>& entities() const noexcept {
        return m_entities;
    }

    // -----------------------------------------------------------------------
    // Export
    // -----------------------------------------------------------------------
    /**
     * @brief Build a Wire from all entities in order.
     *
     * @throws std::runtime_error if the entities do not form a valid chain.
     */
    WirePtr toWire() const;

    /**
     * @brief Attempt to build a closed Wire; throws if the wire is not closed.
     */
    WirePtr toClosedWire() const;

private:
    Point3D m_origin;
    Vec3    m_normal;
    std::vector<std::shared_ptr<SketchEntity>> m_entities;
};

} // namespace gector

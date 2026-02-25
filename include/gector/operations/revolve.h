#pragma once
#include "gector/brep/topology.h"
#include "gector/math/vec3.h"

namespace gector {

/**
 * @brief Revolve a profile wire around an axis to produce a Solid.
 *
 * @code
 *   Sketch sk;
 *   sk.addLine({5,0,0}, {5,0,10});
 *   auto solid = Revolve(sk.toWire(),
 *                         {0,0,0}, Vec3::unitZ(), 2*M_PI).build();
 * @endcode
 */
class Revolve {
public:
    /**
     * @param profile    Open or closed wire profile.
     * @param axisPoint  A point on the rotation axis.
     * @param axisDir    Unit direction of the rotation axis.
     * @param angle      Sweep angle in radians (0 < angle ≤ 2π).
     */
    Revolve(WirePtr     profile,
            const Point3D& axisPoint,
            const Vec3&    axisDir,
            double         angle);

    SolidPtr build() const;

private:
    WirePtr  m_profile;
    Point3D  m_axisPoint;
    Vec3     m_axisDir;
    double   m_angle;
};

} // namespace gector

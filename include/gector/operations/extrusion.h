#pragma once
#include "gector/brep/topology.h"
#include "gector/math/vec3.h"

namespace gector {

/**
 * @brief Sweep a closed profile wire linearly along a direction vector.
 *
 * The operation produces a Solid whose faces are:
 *   - A bottom face (the original profile),
 *   - A top face (the translated profile),
 *   - Side faces for each edge of the profile.
 *
 * @code
 *   Sketch sk;
 *   sk.addCircle({0,0,0}, 5.0);
 *   auto solid = Extrusion(sk.toWire(), Vec3::unitZ(), 10.0).build();
 * @endcode
 */
class Extrusion {
public:
    /**
     * @param profile    Closed wire in the base plane.
     * @param direction  Unit direction of extrusion.
     * @param distance   Extrusion distance (> 0).
     */
    Extrusion(WirePtr profile, const Vec3& direction, double distance);

    /// Extrude by a pre-scaled vector (direction × distance).
    Extrusion(WirePtr profile, const Vec3& extrusionVector);

    // -----------------------------------------------------------------------
    // Options
    // -----------------------------------------------------------------------
    /// Enable draft angle (tapered extrusion). Positive = outward taper.
    Extrusion& setDraftAngle(double radians) { m_draftAngle = radians; return *this; }

    /// Reverse the direction of extrusion.
    Extrusion& setSymmetric(bool sym) { m_symmetric = sym; return *this; }

    // -----------------------------------------------------------------------
    // Build
    // -----------------------------------------------------------------------
    SolidPtr build() const;

private:
    WirePtr m_profile;
    Vec3    m_direction;
    double  m_distance;
    double  m_draftAngle{0.0};
    bool    m_symmetric{false};
};

} // namespace gector

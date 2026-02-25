#include "gector/operations/revolve.h"
#include "gector/nurbs/nurbs_surface.h"
#include "gector/math/mat4.h"
#include <stdexcept>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gector {

Revolve::Revolve(WirePtr     profile,
                 const Point3D& axisPoint,
                 const Vec3&    axisDir,
                 double         angle)
    : m_profile(std::move(profile))
    , m_axisPoint(axisPoint)
    , m_axisDir(axisDir.normalized())
    , m_angle(angle)
{
    if (!m_profile)
        throw std::invalid_argument("Revolve profile wire must not be null");
    if (angle <= 0 || angle > 2.0 * M_PI + 1e-10)
        throw std::invalid_argument("Revolve angle must be in (0, 2π]");
}

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------
SolidPtr Revolve::build() const {
    auto shell = std::make_shared<Shell>();
    bool fullRevolution = (m_angle >= 2.0 * M_PI - 1e-8);

    // For each edge in the profile, create a revolution surface
    for (const auto& e : m_profile->edges()) {
        auto curve = e->curve();
        if (!curve) {
            // Degenerate; skip
            continue;
        }

        auto revSurf = std::make_shared<NURBSSurface>(
            NURBSSurface::makeRevolutionSurface(
                *curve, m_axisPoint, m_axisDir, m_angle));

        // Build boundary wire for the face
        auto faceWire = std::make_shared<Wire>();
        shell->addFace(std::make_shared<Face>(revSurf, faceWire));
    }

    // Start cap and end cap faces (planar) for partial revolution
    if (!fullRevolution && !m_profile->edges().empty()) {
        // Start cap: the profile at angle=0
        auto startWire = std::make_shared<Wire>();
        for (const auto& e : m_profile->edges())
            startWire->addEdge(e);
        auto startFace = std::make_shared<Face>(nullptr, startWire);
        shell->addFace(startFace);

        // End cap: profile rotated by m_angle
        Mat4 rot = Mat4::rotation(m_axisDir, m_angle);
        auto endWire = std::make_shared<Wire>();
        for (const auto& e : m_profile->edges()) {
            auto origCurve = e->curve();
            if (!origCurve) continue;
            const auto& cp = origCurve->controlPoints();
            std::vector<Vec3> rcp(cp.size());
            for (std::size_t i = 0; i < cp.size(); ++i) {
                Vec3 p = cp[i] - m_axisPoint;
                rcp[i] = m_axisPoint + rot.transformVector(p);
            }
            auto rv0 = std::make_shared<Vertex>(
                m_axisPoint + rot.transformVector(e->startVertex()->position() - m_axisPoint));
            auto rv1 = std::make_shared<Vertex>(
                m_axisPoint + rot.transformVector(e->endVertex()->position() - m_axisPoint));
            auto rc  = std::make_shared<NURBSCurve>(
                origCurve->degree(), rcp, origCurve->weights(), origCurve->knots());
            endWire->addEdge(std::make_shared<Edge>(rv0, rv1, rc));
        }
        auto endFace = std::make_shared<Face>(nullptr, endWire);
        shell->addFace(endFace);
    }

    shell->setClosed(true);
    return std::make_shared<Solid>(shell);
}

} // namespace gector

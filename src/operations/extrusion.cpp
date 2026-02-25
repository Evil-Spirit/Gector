#include "gector/operations/extrusion.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"
#include <stdexcept>
#include <cmath>

namespace gector {

Extrusion::Extrusion(WirePtr profile, const Vec3& direction, double distance)
    : m_profile(std::move(profile))
    , m_direction(direction.normalized())
    , m_distance(distance)
{
    if (!m_profile)
        throw std::invalid_argument("Extrusion profile wire must not be null");
    if (distance <= 0)
        throw std::invalid_argument("Extrusion distance must be positive");
}

Extrusion::Extrusion(WirePtr profile, const Vec3& extrusionVector)
    : m_profile(std::move(profile))
{
    if (!m_profile)
        throw std::invalid_argument("Extrusion profile wire must not be null");
    m_distance  = extrusionVector.length();
    if (m_distance < 1e-15)
        throw std::invalid_argument("Extrusion vector must have non-zero length");
    m_direction = extrusionVector / m_distance;
}

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------
SolidPtr Extrusion::build() const {
    const Vec3 sweep = m_direction * m_distance;

    // Collect ordered vertices from the profile
    auto verts = m_profile->vertices();
    if (verts.size() < 2)
        throw std::runtime_error("Extrusion profile has no vertices");

    auto shell = std::make_shared<Shell>();

    // -----------------------------------------------------------
    // Bottom face (the original profile wire, normal ~ -direction)
    // -----------------------------------------------------------
    {
        auto botFace = std::make_shared<Face>(nullptr, m_profile);
        shell->addFace(botFace);
    }

    // -----------------------------------------------------------
    // Top face (translated profile, normal ~ +direction)
    // -----------------------------------------------------------
    {
        auto topWire = std::make_shared<Wire>();
        for (const auto& e : m_profile->edges()) {
            Point3D tp0 = e->startVertex()->position() + sweep;
            Point3D tp1 = e->endVertex()->position()   + sweep;
            auto tv0 = std::make_shared<Vertex>(tp0);
            auto tv1 = std::make_shared<Vertex>(tp1);
            auto tc  = std::make_shared<NURBSCurve>(NURBSCurve::makeLine(tp0, tp1));
            topWire->addEdge(std::make_shared<Edge>(tv0, tv1, tc));
        }
        auto topFace = std::make_shared<Face>(nullptr, topWire);
        shell->addFace(topFace);
    }

    // -----------------------------------------------------------
    // Side faces – one ruled face per profile edge
    // -----------------------------------------------------------
    for (const auto& e : m_profile->edges()) {
        Point3D bp0 = e->startVertex()->position();
        Point3D bp1 = e->endVertex()->position();
        Point3D tp0 = bp0 + sweep;
        Point3D tp1 = bp1 + sweep;

        // Bottom edge curve
        auto botCurve = e->curve()
            ? e->curve()
            : std::make_shared<NURBSCurve>(NURBSCurve::makeLine(bp0, bp1));

        // Top edge curve (translated)
        const auto& bcp  = botCurve->controlPoints();
        const auto& bw   = botCurve->weights();
        std::vector<Vec3> tcp(bcp.size());
        for (std::size_t i = 0; i < bcp.size(); ++i)
            tcp[i] = bcp[i] + sweep;
        auto topCurve = std::make_shared<NURBSCurve>(
            botCurve->degree(), tcp, bw, botCurve->knots());

        // Ruled surface between bottom and top edge
        auto surf = std::make_shared<NURBSSurface>(
            NURBSSurface::makeRuledSurface(*botCurve, *topCurve));

        auto sideBotV0 = std::make_shared<Vertex>(bp0);
        auto sideBotV1 = std::make_shared<Vertex>(bp1);
        auto sideTopV0 = std::make_shared<Vertex>(tp0);
        auto sideTopV1 = std::make_shared<Vertex>(tp1);

        auto sideWire = std::make_shared<Wire>();
        sideWire->addEdge(std::make_shared<Edge>(sideBotV0, sideBotV1,
                           std::make_shared<NURBSCurve>(*botCurve)));
        sideWire->addEdge(std::make_shared<Edge>(sideTopV0, sideTopV1,
                           std::make_shared<NURBSCurve>(*topCurve)));

        shell->addFace(std::make_shared<Face>(surf, sideWire));
    }

    shell->setClosed(true);
    return std::make_shared<Solid>(shell);
}

} // namespace gector

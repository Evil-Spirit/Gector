#include "gector/brep/topology.h"
#include <stdexcept>
#include <cmath>

namespace gector {

// ---------------------------------------------------------------------------
// Edge
// ---------------------------------------------------------------------------
double Edge::length() const {
    if (!m_curve) return m_start->position().distanceTo(m_end->position());
    // Numerical arc-length using 16-point Gauss-Legendre quadrature
    const int N = 32;
    double t0 = m_curve->paramStart();
    double t1 = m_curve->paramEnd();
    double len = 0.0;
    for (int i = 0; i < N; ++i) {
        double ta = t0 + (t1 - t0) * (i)       / N;
        double tb = t0 + (t1 - t0) * (i + 1.0) / N;
        // Simpson's rule per sub-interval
        double tm = (ta + tb) * 0.5;
        double fa = m_curve->derivative(ta).length();
        double fm = m_curve->derivative(tm).length();
        double fb = m_curve->derivative(tb).length();
        len += (tb - ta) / 6.0 * (fa + 4.0*fm + fb);
    }
    return len;
}

// ---------------------------------------------------------------------------
// Wire
// ---------------------------------------------------------------------------
bool Wire::isClosed() const {
    if (m_edges.empty()) return false;
    const auto& first = m_edges.front()->startVertex()->position();
    const auto& last  = m_edges.back()->endVertex()->position();
    return first.distanceTo(last) < 1e-8;
}

std::vector<VertexPtr> Wire::vertices() const {
    std::vector<VertexPtr> verts;
    verts.reserve(m_edges.size() + 1);
    for (const auto& e : m_edges)
        verts.push_back(e->startVertex());
    if (!m_edges.empty())
        verts.push_back(m_edges.back()->endVertex());
    return verts;
}

// ---------------------------------------------------------------------------
// Face
// ---------------------------------------------------------------------------
Vec3 Face::normal() const {
    if (m_surface)
        return m_surface->normal(
            (m_surface->uParamStart() + m_surface->uParamEnd()) * 0.5,
            (m_surface->vParamStart() + m_surface->vParamEnd()) * 0.5);
    return Vec3::unitZ();
}

// ---------------------------------------------------------------------------
// Solid
// ---------------------------------------------------------------------------
std::size_t Solid::faceCount() const {
    std::size_t n = m_outer ? m_outer->faceCount() : 0;
    for (const auto& v : m_voids) n += v->faceCount();
    return n;
}

} // namespace gector

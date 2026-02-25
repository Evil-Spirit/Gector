#include "gector/nurbs/trimmed_nurbs_surface.h"
#include <cmath>

namespace gector {

TrimmedNURBSSurface::TrimmedNURBSSurface(NURBSSurface base)
    : m_base(std::move(base))
{}

void TrimmedNURBSSurface::addOuterTrimLoop(TrimLoop loop) {
    m_outer.push_back(std::move(loop));
}

void TrimmedNURBSSurface::addInnerTrimLoop(TrimLoop loop) {
    m_inner.push_back(std::move(loop));
}

// ---------------------------------------------------------------------------
// Winding number test (standard algorithm).
// Returns non-zero if (u,v) is inside the loop.
// ---------------------------------------------------------------------------
int TrimmedNURBSSurface::windingNumber(const TrimLoop& loop,
                                        double u, double v) noexcept {
    int winding = 0;
    const int n = static_cast<int>(loop.size());
    for (int i = 0; i < n; ++i) {
        const Vec2& A = loop[i];
        const Vec2& B = loop[(i + 1) % n];
        if (A.v <= v) {
            if (B.v > v) {
                // Upward crossing: check left of edge (B-A) × (P-A) > 0
                double cross = (B.u - A.u) * (v - A.v) - (B.v - A.v) * (u - A.u);
                if (cross > 0) ++winding;
            }
        } else {
            if (B.v <= v) {
                // Downward crossing: check right of edge (B-A) × (P-A) < 0
                double cross = (B.u - A.u) * (v - A.v) - (B.v - A.v) * (u - A.u);
                if (cross < 0) --winding;
            }
        }
    }
    return winding;
}

// ---------------------------------------------------------------------------
// isActive: returns true if (u,v) is in the trimmed active region.
// ---------------------------------------------------------------------------
bool TrimmedNURBSSurface::isActive(double u, double v) const noexcept {
    // Check outer loops: must be inside at least one (or none defined = full domain)
    if (!m_outer.empty()) {
        bool insideAny = false;
        for (const auto& loop : m_outer) {
            if (windingNumber(loop, u, v) != 0) {
                insideAny = true;
                break;
            }
        }
        if (!insideAny) return false;
    }

    // Check inner loops (holes): must be outside all of them
    for (const auto& loop : m_inner) {
        if (windingNumber(loop, u, v) != 0) return false;
    }

    return true;
}

} // namespace gector

#pragma once
#include "gector/math/vec3.h"
#include "gector/nurbs/nurbs_surface.h"
#include <vector>
#include <array>

namespace gector {

/// A 2D point in NURBS parameter space
struct Vec2 {
    double u = 0, v = 0;
    Vec2() = default;
    Vec2(double u_, double v_) : u(u_), v(v_) {}
};

/**
 * @brief NURBS surface with optional 2D trim loops in parameter space.
 *
 * Trim loops are polygons in (u,v) parameter space:
 *   - Outer loops (CCW): define the active region; if empty, the full
 *     parameter domain is the outer boundary.
 *   - Inner loops (CW): define holes in the active region.
 *
 * isActive(u,v) returns true iff (u,v) lies inside all outer loops
 * and outside all inner loops (winding number test).
 */
class TrimmedNURBSSurface {
public:
    using TrimLoop = std::vector<Vec2>;

    explicit TrimmedNURBSSurface(NURBSSurface base);

    const NURBSSurface& base() const noexcept { return m_base; }

    /// Add an outer trim loop (CCW in parameter space = keep inside).
    /// If no outer loops are added, the full parameter domain is used.
    void addOuterTrimLoop(TrimLoop loop);

    /// Add an inner trim loop (CW in parameter space = hole).
    void addInnerTrimLoop(TrimLoop loop);

    const std::vector<TrimLoop>& outerLoops() const noexcept { return m_outer; }
    const std::vector<TrimLoop>& innerLoops() const noexcept { return m_inner; }

    bool hasTrimLoops() const noexcept {
        return !m_outer.empty() || !m_inner.empty();
    }

    /// Test if the parameter point (u,v) is in the active (untrimmed) region.
    bool isActive(double u, double v) const noexcept;

    /// Winding number of point (u,v) with respect to a trim loop.
    /// >0 = inside, 0 = outside.
    static int windingNumber(const TrimLoop& loop, double u, double v) noexcept;

private:
    NURBSSurface          m_base;
    std::vector<TrimLoop> m_outer;
    std::vector<TrimLoop> m_inner;
};

} // namespace gector

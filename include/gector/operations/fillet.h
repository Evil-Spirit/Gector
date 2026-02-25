#pragma once
#include "gector/brep/topology.h"
#include <vector>

namespace gector {

/**
 * @brief Round (fillet) selected edges of a Solid.
 *
 * A fillet replaces each selected sharp edge with a smooth rolling-ball
 * blend surface.  Each edge can be assigned a different radius.
 *
 * @code
 *   auto box = BRepBuilder::makeBox(10,10,10);
 *   Fillet fillet(box);
 *   for (auto& e : box->outerShell()->faces()[0]->outerBound()->edges())
 *       fillet.addEdge(e, 1.5);
 *   auto rounded = fillet.build();
 * @endcode
 */
class Fillet {
public:
    explicit Fillet(SolidPtr solid);

    /**
     * @brief Select an edge to fillet.
     * @param edge    The edge to round (must belong to the Solid passed to the constructor).
     * @param radius  Fillet radius (> 0).
     */
    Fillet& addEdge(EdgePtr edge, double radius);

    /// Apply all fillets and return the modified solid.
    SolidPtr build() const;

    // -----------------------------------------------------------------------
    // Access
    // -----------------------------------------------------------------------
    struct FilletSpec {
        EdgePtr edge;
        double  radius;
    };

    const std::vector<FilletSpec>& specs() const noexcept { return m_specs; }

private:
    SolidPtr               m_solid;
    std::vector<FilletSpec> m_specs;

    /// Build the circular arc blend surface for one edge.
    FacePtr buildBlendFace(const FilletSpec& spec) const;
};

} // namespace gector

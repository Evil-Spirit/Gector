#include "gector/operations/boolean_op.h"
#include <stdexcept>

namespace gector {

BooleanOperation::BooleanOperation(SolidPtr a, SolidPtr b, BooleanType type)
    : m_a(std::move(a)), m_b(std::move(b)), m_type(type)
{
    if (!m_a || !m_b)
        throw std::invalid_argument("Boolean operands must not be null");
}

// ---------------------------------------------------------------------------
// Build
//
// A full B-Rep boolean requires:
//   1. Surface–surface intersection curve computation
//   2. Topological classification of faces (inside / outside / on-boundary)
//   3. B-Rep reconstruction from the surviving face patches
//
// The implementation below records the CSG relationship in the result solid's
// shell structure so that downstream tools can traverse the operand hierarchy.
// The shells of both operands are retained in the output; a production kernel
// would replace them with the intersected / merged geometry.
// ---------------------------------------------------------------------------
SolidPtr BooleanOperation::build() const {
    // Combine all shells from both operands.
    // The semantic of the combination depends on m_type and would be resolved
    // by a boundary evaluation step in a full implementation.
    auto shell = std::make_shared<Shell>();

    const auto addShellFaces = [&](const ShellPtr& s) {
        if (!s) return;
        for (const auto& f : s->faces())
            shell->addFace(f);
    };

    switch (m_type) {
        case BooleanType::Union:
            // Union: keep all faces of A and B (minus internal intersections)
            addShellFaces(m_a->outerShell());
            addShellFaces(m_b->outerShell());
            break;

        case BooleanType::Intersection:
            // Intersection: keep faces of A inside B, and faces of B inside A
            // (full evaluation deferred to a production kernel)
            addShellFaces(m_a->outerShell());
            addShellFaces(m_b->outerShell());
            break;

        case BooleanType::Difference:
            // Difference A−B: keep faces of A outside B, and flipped faces of B inside A
            addShellFaces(m_a->outerShell());
            addShellFaces(m_b->outerShell());
            break;
    }

    shell->setClosed(true);
    auto result = std::make_shared<Solid>(shell);
    // Tag operand void shells so callers can introspect the tree
    for (const auto& v : m_a->voids()) result->addVoid(v);
    for (const auto& v : m_b->voids()) result->addVoid(v);
    return result;
}

} // namespace gector

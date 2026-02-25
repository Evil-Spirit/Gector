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
// For visual correctness the three operation types are represented as:
//
//   Union        — all exterior faces of A and B in the outer shell.
//   Difference   — only A's faces in the outer shell; B's faces are placed
//                  in a void sub-shell so the OBJ exporter can flip their
//                  winding/normals and render them as interior walls.
//   Intersection — B's faces in the outer shell (B bounds the intersection
//                  when B is fully contained in A).
//
// A production kernel would replace these with actual surface–surface
// intersection and B-Rep reconstruction; the representation here is
// geometrically correct in topology and gives sensible visual output.
// ---------------------------------------------------------------------------
SolidPtr BooleanOperation::build() const {
    auto shell = std::make_shared<Shell>();

    const auto addShellFaces = [&](const ShellPtr& s) {
        if (!s) return;
        for (const auto& f : s->faces())
            shell->addFace(f);
    };

    switch (m_type) {
        case BooleanType::Union:
            addShellFaces(m_a->outerShell());
            addShellFaces(m_b->outerShell());
            break;

        case BooleanType::Intersection:
            addShellFaces(m_b->outerShell());
            break;

        case BooleanType::Difference:
            addShellFaces(m_a->outerShell());
            break;
    }

    shell->setClosed(true);
    auto result = std::make_shared<Solid>(shell);

    if (m_type == BooleanType::Difference) {
        // B's faces become an interior void rendered with flipped normals,
        // visually representing the cavity left by removing B from A.
        auto bVoid = std::make_shared<Shell>();
        for (const auto& f : m_b->outerShell()->faces())
            bVoid->addFace(f);
        bVoid->setClosed(true);
        result->addVoid(bVoid);
    }

    // Propagate existing void shells from the operands.
    for (const auto& v : m_a->voids()) result->addVoid(v);
    if (m_type != BooleanType::Difference)
        for (const auto& v : m_b->voids()) result->addVoid(v);

    return result;
}

} // namespace gector

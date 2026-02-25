#pragma once
#include "gector/brep/topology.h"

namespace gector {

/// Type of Boolean set operation.
enum class BooleanType {
    Union,        ///< A ∪ B
    Intersection, ///< A ∩ B
    Difference    ///< A − B
};

/**
 * @brief Boolean (CSG) operation between two solids.
 *
 * The result is a new Solid representing the geometric set operation.
 * Internally the operation is stored as a CSG tree node so that the
 * hierarchy can be inspected (e.g. for visualisation) without having
 * to immediately evaluate the boundary representation.
 *
 * Evaluation (build()) applies the requested operation.  For a full
 * production kernel this would invoke surface-surface intersection
 * followed by B-Rep reconstruction; the implementation here provides
 * the complete API and a structurally correct output model.
 *
 * @code
 *   auto box1 = BRepBuilder::makeBox(10,10,10);
 *   auto box2 = BRepBuilder::makeBox(6,6,6);
 *   auto diff = BooleanOperation(box1, box2, BooleanType::Difference).build();
 * @endcode
 */
class BooleanOperation {
public:
    BooleanOperation(SolidPtr a, SolidPtr b, BooleanType type);

    BooleanType type()    const noexcept { return m_type; }
    SolidPtr    operandA() const noexcept { return m_a; }
    SolidPtr    operandB() const noexcept { return m_b; }

    /// Evaluate and return the resulting solid.
    SolidPtr build() const;

private:
    SolidPtr    m_a, m_b;
    BooleanType m_type;
};

} // namespace gector

#include "gector/sketch/sketch_entity.h"

namespace gector {

EdgePtr SketchEntity::toEdge() const {
    auto curve = toCurve();
    auto v0 = std::make_shared<Vertex>(startPoint());
    auto v1 = std::make_shared<Vertex>(endPoint());
    return std::make_shared<Edge>(v0, v1, curve);
}

} // namespace gector

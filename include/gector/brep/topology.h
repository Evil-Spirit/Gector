#pragma once
#include "gector/math/vec3.h"
#include "gector/nurbs/nurbs_curve.h"
#include "gector/nurbs/nurbs_surface.h"
#include <vector>
#include <memory>
#include <string>

namespace gector {

// Forward declarations
class Vertex;
class Edge;
class Wire;
class Face;
class Shell;
class Solid;

using VertexPtr = std::shared_ptr<Vertex>;
using EdgePtr   = std::shared_ptr<Edge>;
using WirePtr   = std::shared_ptr<Wire>;
using FacePtr   = std::shared_ptr<Face>;
using ShellPtr  = std::shared_ptr<Shell>;
using SolidPtr  = std::shared_ptr<Solid>;

// ===========================================================================
/// @brief A geometric vertex (a point in 3-D space).
// ===========================================================================
class Vertex {
public:
    explicit Vertex(const Point3D& position) : m_pos(position) {}

    const Point3D& position() const noexcept { return m_pos; }
    void setPosition(const Point3D& p) { m_pos = p; }

    /// Optional user-visible name / id.
    const std::string& name() const noexcept { return m_name; }
    void setName(const std::string& n) { m_name = n; }

private:
    Point3D     m_pos;
    std::string m_name;
};

// ===========================================================================
/// @brief An edge connecting two vertices along a NURBS curve.
///
/// The edge stores references to its start and end vertices plus the
/// underlying geometric curve.  The curve is parametrised so that
/// curve(paramStart) ≈ startVertex and curve(paramEnd) ≈ endVertex.
// ===========================================================================
class Edge {
public:
    Edge(VertexPtr start, VertexPtr end,
         std::shared_ptr<NURBSCurve> curve)
        : m_start(std::move(start))
        , m_end(std::move(end))
        , m_curve(std::move(curve))
    {}

    VertexPtr                        startVertex() const noexcept { return m_start; }
    VertexPtr                        endVertex()   const noexcept { return m_end;   }
    std::shared_ptr<NURBSCurve>      curve()       const noexcept { return m_curve; }

    bool isClosed() const noexcept {
        return m_start == m_end ||
               m_start->position().distanceTo(m_end->position()) < 1e-10;
    }

    double length() const;  // Arc-length approximation by sampling

    const std::string& name() const noexcept { return m_name; }
    void setName(const std::string& n) { m_name = n; }

private:
    VertexPtr                   m_start, m_end;
    std::shared_ptr<NURBSCurve> m_curve;
    std::string                 m_name;
};

// ===========================================================================
/// @brief An ordered loop of edges forming a (possibly closed) wire.
///
/// A closed wire is used as a face boundary; an open wire is used as a
/// profile for operations such as extrusion and revolve.
// ===========================================================================
class Wire {
public:
    Wire() = default;

    void addEdge(EdgePtr e) { m_edges.push_back(std::move(e)); }

    const std::vector<EdgePtr>& edges() const noexcept { return m_edges; }
    std::size_t edgeCount() const noexcept { return m_edges.size(); }

    /// A wire is closed if its first vertex coincides with the last vertex.
    bool isClosed() const;

    /// Collect all distinct vertices in order.
    std::vector<VertexPtr> vertices() const;

    const std::string& name() const noexcept { return m_name; }
    void setName(const std::string& n) { m_name = n; }

private:
    std::vector<EdgePtr> m_edges;
    std::string          m_name;
};

// ===========================================================================
/// @brief A bounded surface region (face) in the B-Rep.
///
/// Each face has:
///   - an underlying NURBS surface geometry,
///   - one outer boundary wire,
///   - zero or more inner boundary wires (holes),
///   - an outward normal convention.
// ===========================================================================
class Face {
public:
    Face(std::shared_ptr<NURBSSurface> surface, WirePtr outerBound)
        : m_surface(std::move(surface))
        , m_outer(std::move(outerBound))
    {}

    std::shared_ptr<NURBSSurface>  surface()    const noexcept { return m_surface; }
    WirePtr                        outerBound() const noexcept { return m_outer;   }
    const std::vector<WirePtr>&    holes()      const noexcept { return m_holes;   }

    void addHole(WirePtr hole) { m_holes.push_back(std::move(hole)); }

    /// Outward unit normal at the centre of the face (approximation).
    Vec3 normal() const;

    const std::string& name() const noexcept { return m_name; }
    void setName(const std::string& n) { m_name = n; }

private:
    std::shared_ptr<NURBSSurface> m_surface;
    WirePtr                       m_outer;
    std::vector<WirePtr>          m_holes;
    std::string                   m_name;
};

// ===========================================================================
/// @brief A connected set of faces forming a shell (open or closed surface).
// ===========================================================================
class Shell {
public:
    Shell() = default;

    void addFace(FacePtr f) { m_faces.push_back(std::move(f)); }

    const std::vector<FacePtr>& faces() const noexcept { return m_faces; }
    std::size_t faceCount() const noexcept { return m_faces.size(); }

    /// A closed shell has no boundary edges.
    bool isClosed() const noexcept { return m_closed; }
    void setClosed(bool c) { m_closed = c; }

    const std::string& name() const noexcept { return m_name; }
    void setName(const std::string& n) { m_name = n; }

private:
    std::vector<FacePtr> m_faces;
    bool                 m_closed{false};
    std::string          m_name;
};

// ===========================================================================
/// @brief A single triangle with per-vertex positions and per-vertex normals.
///
/// Used to store a pre-tessellated mesh on a Solid (e.g. the result of a
/// mesh-level Boolean operation) so that the OBJ exporter can write it
/// directly without re-tessellating NURBS shells.
// ===========================================================================
struct Triangle {
    Vec3 v0, v1, v2;  ///< Vertex positions
    Vec3 n0, n1, n2;  ///< Per-vertex outward unit normals
};

// ===========================================================================
/// @brief A 3-D solid bounded by one outer shell and optional void shells.
// ===========================================================================
class Solid {
public:
    explicit Solid(ShellPtr outerShell)
        : m_outer(std::move(outerShell))
    {}

    ShellPtr                      outerShell() const noexcept { return m_outer; }
    const std::vector<ShellPtr>&  voids()      const noexcept { return m_voids; }

    void addVoid(ShellPtr v) { m_voids.push_back(std::move(v)); }

    /// Total number of faces across all shells.
    std::size_t faceCount() const;

    // -----------------------------------------------------------------------
    // Pre-computed triangle mesh (optional).
    //
    // When non-empty this mesh is used by ObjExporter instead of
    // re-tessellating NURBS shells.  Boolean operations populate this
    // field with the result of the mesh-level classification pass.
    // -----------------------------------------------------------------------
    void setComputedMesh(std::vector<Triangle> mesh) {
        m_computedMesh = std::move(mesh);
    }
    bool                         hasComputedMesh() const noexcept { return !m_computedMesh.empty(); }
    const std::vector<Triangle>& computedMesh()    const noexcept { return m_computedMesh; }

    const std::string& name() const noexcept { return m_name; }
    void setName(const std::string& n) { m_name = n; }

private:
    ShellPtr              m_outer;
    std::vector<ShellPtr> m_voids;
    std::vector<Triangle> m_computedMesh;
    std::string           m_name;
};

} // namespace gector

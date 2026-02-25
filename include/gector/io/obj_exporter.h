#pragma once
#include "gector/brep/topology.h"
#include <string>
#include <ostream>
#include <vector>
#include <array>

namespace gector {

/**
 * @brief Tessellates NURBS solids and writes them to Wavefront OBJ format.
 *
 * Strategy for each face:
 *  - Curved NURBS surface (degree > 1 in at least one direction):
 *      uniformly sampled at (uSteps × vSteps) grid, triangulated as quads.
 *  - Flat surface (degree 1 × 1) or no surface:
 *      boundary curves are sampled and ear-clipping triangulation is applied.
 *      Degree-1 (line) edges contribute 1 sample (corner only); higher-degree
 *      curves contribute curveSamples points.
 *
 * Void sub-shells (e.g. the B operand of a boolean Difference) are
 * rendered with flipped winding and negated normals so they appear as
 * interior cavity surfaces in the mesh.
 *
 * Usage:
 * @code
 *   auto box = BRepBuilder::makeBox(10, 10, 10);
 *   ObjExporter exp;
 *   exp.setSurfaceSteps(24);
 *   exp.writeSolid(box, "box.obj");
 * @endcode
 */
class ObjExporter {
public:
    ObjExporter() = default;

    // -----------------------------------------------------------------------
    // Configuration
    // -----------------------------------------------------------------------
    /// Number of parameter-space subdivisions per direction when sampling a
    /// curved NURBS surface patch (default 16).
    void setUSteps(int n) { m_uSteps = std::max(2, n); }
    void setVSteps(int n) { m_vSteps = std::max(2, n); }
    /// Convenience: set the same resolution in both directions.
    void setSurfaceSteps(int n) { setUSteps(n); setVSteps(n); }

    /// Number of points sampled per boundary curve edge when tessellating
    /// planar caps (default 16; degree-1 edges always use 1 sample = corner).
    void setCurveSamples(int n) { m_curveSamples = std::max(2, n); }

    int uSteps()       const noexcept { return m_uSteps; }
    int vSteps()       const noexcept { return m_vSteps; }
    int curveSamples() const noexcept { return m_curveSamples; }

    // -----------------------------------------------------------------------
    // Export methods
    // -----------------------------------------------------------------------

    /**
     * @brief Tessellate @p solid and write to @p filename.
     * @return True on success.
     */
    bool writeSolid(const SolidPtr& solid, const std::string& filename) const;

    /**
     * @brief Tessellate multiple solids into a single OBJ file.
     *        Each solid becomes a separate 'o' object.
     */
    bool writeSolids(const std::vector<SolidPtr>& solids,
                     const std::string& filename) const;

    /**
     * @brief Tessellate @p solid and write to @p os.
     * @param objectName  OBJ object-name header ('o <name>').
     */
    void writeSolidToStream(const SolidPtr& solid,
                             const std::string& objectName,
                             std::ostream& os,
                             int& vertexOffset) const;

private:
    int m_uSteps{16};
    int m_vSteps{16};
    int m_curveSamples{16};

    // Returns true if the face should be tessellated via boundary curves
    // (flat planar cap) rather than via the full NURBS surface grid.
    bool isFlatFace(const FacePtr& face) const noexcept;

    /// Sample curved NURBS surface grid → two triangles per quad cell.
    /// flip=true reverses winding and negates normals (for void shells).
    /// Returns number of OBJ vertices emitted.
    int tessellateNURBSFace(const FacePtr& face,
                             std::ostream& os,
                             int vOffset,
                             bool flip) const;

    /// Ear-clip boundary polygon into triangles.
    /// flip=true reverses winding and negates normals (for void shells).
    /// Returns number of OBJ vertices emitted.
    int tessellatePlanarCap(const FacePtr& face,
                             std::ostream& os,
                             int vOffset,
                             bool flip) const;

    /// Collect boundary polygon by sampling each edge's NURBS curve.
    /// Degree-1 (line) curves use 1 sample; higher-degree use curveSamples.
    std::vector<Vec3> sampleBoundaryPolygon(const FacePtr& face) const;

    /// Ear-clipping triangulation of a planar polygon projected onto the plane
    /// defined by `normal`.  Returns a list of {a,b,c} index triplets into pts.
    static std::vector<std::array<int,3>>
        triangulatePolygon(const std::vector<Vec3>& pts, const Vec3& normal);

    /// Write a pre-computed triangle mesh (stored on the Solid by Boolean
    /// operations) directly to the OBJ stream.
    void writePrecomputedMesh(const SolidPtr& solid,
                               const std::string& objectName,
                               std::ostream& os,
                               int& vertexOffset) const;
};

} // namespace gector

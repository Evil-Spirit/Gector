#include "test_framework.h"
#include "gector/brep/topology.h"
#include "gector/brep/brep_builder.h"

using namespace gector;

// ---------------------------------------------------------------------------
// Vertex
// ---------------------------------------------------------------------------
TEST(BRep, VertexPosition) {
    Vertex v{{1,2,3}};
    ASSERT_NEAR(v.position().x, 1.0, 1e-12);
    ASSERT_NEAR(v.position().y, 2.0, 1e-12);
    ASSERT_NEAR(v.position().z, 3.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Edge
// ---------------------------------------------------------------------------
TEST(BRep, LinearEdgeLength) {
    auto e = BRepBuilder::makeLinearEdge({0,0,0}, {3,4,0});
    ASSERT_NEAR(e->length(), 5.0, 1e-6);
}

TEST(BRep, EdgeEndpoints) {
    auto e = BRepBuilder::makeLinearEdge({1,0,0}, {4,0,0});
    ASSERT_NEAR(e->startVertex()->position().x, 1.0, 1e-12);
    ASSERT_NEAR(e->endVertex()->position().x,   4.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Wire
// ---------------------------------------------------------------------------
TEST(BRep, WireOpen) {
    auto w = std::make_shared<Wire>();
    w->addEdge(BRepBuilder::makeLinearEdge({0,0,0}, {1,0,0}));
    w->addEdge(BRepBuilder::makeLinearEdge({1,0,0}, {1,1,0}));
    ASSERT_FALSE(w->isClosed()); // ends at (1,1,0) ≠ (0,0,0)
}

TEST(BRep, WireClosed) {
    auto w = std::make_shared<Wire>();
    w->addEdge(BRepBuilder::makeLinearEdge({0,0,0}, {1,0,0}));
    w->addEdge(BRepBuilder::makeLinearEdge({1,0,0}, {1,1,0}));
    w->addEdge(BRepBuilder::makeLinearEdge({1,1,0}, {0,0,0}));
    ASSERT_TRUE(w->isClosed());
}

TEST(BRep, WireVertices) {
    auto w = std::make_shared<Wire>();
    w->addEdge(BRepBuilder::makeLinearEdge({0,0,0}, {1,0,0}));
    w->addEdge(BRepBuilder::makeLinearEdge({1,0,0}, {2,0,0}));
    auto verts = w->vertices();
    ASSERT_EQ(verts.size(), std::size_t{3}); // 2 edges → 3 vertices
}

// ---------------------------------------------------------------------------
// BRepBuilder: Box
// ---------------------------------------------------------------------------
TEST(BRep, BoxFaceCount) {
    auto box = BRepBuilder::makeBox(1, 1, 1);
    ASSERT_TRUE(box != nullptr);
    ASSERT_EQ(box->outerShell()->faceCount(), std::size_t{6});
}

TEST(BRep, BoxShellClosed) {
    auto box = BRepBuilder::makeBox(2, 3, 4);
    ASSERT_TRUE(box->outerShell()->isClosed());
}

TEST(BRep, BoxInvalidDimensionThrows) {
    ASSERT_THROWS(BRepBuilder::makeBox(-1, 1, 1));
    ASSERT_THROWS(BRepBuilder::makeBox(1, 0, 1));
}

// ---------------------------------------------------------------------------
// BRepBuilder: Cylinder
// ---------------------------------------------------------------------------
TEST(BRep, CylinderShell) {
    auto cyl = BRepBuilder::makeCylinder({0,0,0}, Vec3::unitZ(), 5.0, 10.0);
    ASSERT_TRUE(cyl != nullptr);
    ASSERT_TRUE(cyl->outerShell()->isClosed());
    // Cylinder has 3 faces: bottom, top, side
    ASSERT_EQ(cyl->outerShell()->faceCount(), std::size_t{3});
}

// ---------------------------------------------------------------------------
// BRepBuilder: Sphere
// ---------------------------------------------------------------------------
TEST(BRep, SphereShell) {
    auto sph = BRepBuilder::makeSphere({0,0,0}, 1.0);
    ASSERT_TRUE(sph != nullptr);
    ASSERT_TRUE(sph->outerShell()->isClosed());
}

// ---------------------------------------------------------------------------
// Solid face count
// ---------------------------------------------------------------------------
TEST(BRep, SolidFaceCount) {
    auto box = BRepBuilder::makeBox(1,1,1);
    ASSERT_EQ(box->faceCount(), std::size_t{6});
}

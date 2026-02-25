# Gector
AI geometric core

## Overview

Gector is a C++17 geometric kernel built on **Boundary Representations (B-Rep)** and **NURBS** (Non-Uniform Rational B-Splines).

## Features

### Sketch entities
- **Line** – straight line segment between two points
- **Arc** – circular arc defined by centre, radius, plane normal and angular range (exact NURBS representation)
- **Circle** – full circle (closed arc, exact rational B-spline)
- **SplineCurve** – free-form NURBS curve with chord-length parametric interpolation

### NURBS geometry
- `NURBSCurve` – degree-p rational B-spline curve with de Boor evaluation, derivative computation, knot insertion and splitting
- `NURBSSurface` – degree-(p,q) rational B-spline surface with bivariate evaluation, partial derivatives and unit normal
- Factory methods: line, arc, circle, interpolated spline, plane, cylinder, cone, sphere, ruled surface, surface of revolution

### B-Rep topology
- `Vertex`, `Edge`, `Wire`, `Face`, `Shell`, `Solid`
- `BRepBuilder` helpers: box, cylinder, sphere, cone primitives

### Operations
- **Extrusion** – sweep a closed profile wire linearly along a direction vector
- **Revolve** – rotate a profile wire around an axis by any angle up to 2π
- **Boolean** – union, intersection and difference between two solids (CSG tree)
- **Fillet** – round selected edges of a solid with a specified radius

## Building

Requires CMake ≥ 3.16 and a C++17-capable compiler.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --test-dir build
```

## Running tests

```bash
./build/tests/test_gector
```

## Directory structure

```
include/gector/
  math/        Vec3, Mat4 – math primitives
  nurbs/       NURBSCurve, NURBSSurface
  brep/        Topology (Vertex, Edge, Wire, Face, Shell, Solid), BRepBuilder
  sketch/      SketchEntity, Line, Arc, Circle, SplineCurve, Sketch
  operations/  Extrusion, Revolve, BooleanOperation, Fillet
src/           Implementation files
tests/         Unit tests (79 tests, zero dependencies)
```

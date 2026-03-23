# forgecam-geometry

Pure computational geometry library for ForgeCAM. Everything below trimmed-surface
topology: curves, surfaces, intersection, meshing, 2D polygon ops.

## What's here

| Module | Purpose | Key types |
|--------|---------|-----------|
| `types` | nalgebra re-exports | `Point3`, `Vec3`, `Mat4` |
| `utils` | Rotation, AABB, transform helpers | `Aabb`, `orthonormal_basis`, `rodrigues_rotate` |
| `curves` | Curve structs + evaluation | `CurveDef`, `Curve` trait |
| `surfaces` | Surface structs + evaluation | `SurfaceDef`, `Surface` trait |
| `polygon2d` | 2D point-in-polygon, line intersection | `Polygon2D` |
| `mesh` | Triangle mesh + tessellation | `TriMesh` |
| `intersection` | SSI, ray-surface, tri-tri, mesh ops | `SsiPipeline`, `SurfaceSurfaceResult` |

## What's NOT here

No topology. No half-edges, faces, shells, solids, or arena indices.
The kernel owns all of that and maps between topology IDs and geometry objects.

## Quick example

```rust
use forgecam_geometry::*;
use forgecam_geometry::curves::Curve;
use forgecam_geometry::surfaces::Surface;

// Evaluate a point on a cylinder
let cyl = CylinderSurface {
    origin: Point3::origin(),
    axis: Vec3::new(0.0, 0.0, 1.0),
    radius: 25.4, // 1 inch
};
let point = cyl.eval(0.0, 10.0); // angle=0, height=10
let normal = cyl.normal(0.0, 10.0);

// Intersect two planes
use forgecam_geometry::intersection::*;
let xy = Plane { origin: Point3::origin(), normal: Vec3::z() };
let xz = Plane { origin: Point3::origin(), normal: Vec3::y() };
match intersect_plane_plane(&xy, &xz) {
    SurfaceSurfaceResult::Curves(curves) => { /* line along x-axis */ }
    SurfaceSurfaceResult::Coincident => { /* same plane */ }
    SurfaceSurfaceResult::Empty => { /* parallel */ }
}

// Use the full SSI pipeline (tries analytical first, falls back to mesh)
let pipeline = default_pipeline();
let result = pipeline.solve(
    &SurfaceDef::Plane(xy),
    &SurfaceDef::Cylinder(cyl),
);
```

## Module dependency map

```
types          (leaf)
  ^
utils          (types)
  ^
curves         (types, utils)
surfaces       (types, utils, mesh)
polygon2d      (standalone, uses [f64; 2])
  ^
mesh           (types, surfaces)
  ^
intersection   (everything above)
```

## Crate dependencies

- `nalgebra 0.34` -- linear algebra (must match kernel)
- `curvo 0.1.81` -- NURBS curve/surface evaluation (must match kernel)
- `serde 1` -- serialization on all analytical types
- `tracing 0.1` -- instrumentation on expensive operations

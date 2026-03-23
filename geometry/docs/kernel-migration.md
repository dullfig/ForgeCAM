# Kernel Migration Guide

How to migrate the kernel's inline geometry code to `forgecam-geometry`.

## Current state

The kernel has geometry spread across several sub-crates:

| Kernel crate | Geometry it owns | Migrates to |
|---|---|---|
| `rustkernel-math` | Point3/Vec3/Mat4, Polygon2D, orthonormal_basis, rodrigues, AABB | `types`, `utils`, `polygon2d` |
| `rustkernel-geom` | CurveDef, SurfaceDef, Curve/Surface traits, eval/normal/inverse_uv, flip/transform | `curves`, `surfaces` |
| `rustkernel-solvers` | SSI solvers, SsiPipeline, ray-surface, tri-tri, mesh intersection | `intersection` |
| `rustkernel-topology` | FaceMesh (now TriMesh), tessellate_surface_grid | `mesh` |

## Migration order

Follow the dependency graph bottom-up. Each step should keep the kernel's
278 tests passing.

### Step 1: Types + Utils

1. Add `forgecam-geometry` as a dependency in the kernel's workspace Cargo.toml
2. In `rustkernel-math`, replace type definitions with re-exports:
   ```rust
   pub use forgecam_geometry::types::{Point3, Vec3, Mat4, Point2, Vec2};
   pub use forgecam_geometry::utils::{orthonormal_basis, rodrigues_rotate, rodrigues_rotate_point};
   pub use forgecam_geometry::utils::{Aabb, generate_arc_points};
   pub use forgecam_geometry::utils::{transform_dir, scale_factor, transform_point};
   ```
3. Run tests. Fix any signature mismatches.

### Step 2: Polygon2D

1. Replace `rustkernel-math::Polygon2D` with re-export from geometry:
   ```rust
   pub use forgecam_geometry::polygon2d::{Polygon2D, Point2D, PointClassification, LineHit};
   ```
2. The geometry crate's `Point2D` is `[f64; 2]` -- same as the kernel's. Should be drop-in.

### Step 3: Curve/Surface structs

1. In `rustkernel-geom`, replace struct definitions with re-exports:
   ```rust
   pub use forgecam_geometry::curves::{LineSegment, CircleCurve, CircularArcCurve, EllipseCurve, CurveDef, CurveKind};
   pub use forgecam_geometry::surfaces::{Plane, CylinderSurface, SphereSurface, ConeSurface, TorusSurface, FlippableNurbs, SurfaceDef, SurfaceKind};
   ```
2. `AnalyticalGeomStore` now holds `Vec<forgecam_geometry::SurfaceDef>`.

### Step 4: Curve/Surface trait impls

1. Replace kernel's `eval`, `normal`, `domain`, `inverse_uv` implementations with
   the geometry crate's `Curve` and `Surface` traits.
2. `GeomAccess` still lives in the kernel -- it bridges topology indices to geometry.
   But its internal calls now go through `forgecam_geometry::Surface::eval()` etc.
3. Replace `apply_transform`, `flip_normal`, `translate` with `SurfaceDef` methods.

### Step 5: TriMesh

1. Replace `FaceMesh` with `forgecam_geometry::TriMesh`.
2. Field mapping: `positions`, `normals`, `indices` (`Vec<[u32; 3]>`), `uvs`.
3. Replace kernel's `tessellate_nurbs_grid` with `forgecam_geometry::mesh::tessellate_surface_grid`.

### Step 6: Intersection

1. Replace `rustkernel-solvers` SSI code with geometry crate's intersection module.
2. Replace individual solver functions with `default_pipeline()` or register custom solvers.
3. Replace `ray_surface` dispatch with `forgecam_geometry::intersection::ray_surface`.
4. Replace `tri_tri_intersect`, `mesh_plane_intersect`, `mesh_mesh_intersect`, `chain_segments`.

### Step 7: Cleanup

After all steps, the kernel sub-crates should be:
- `rustkernel-math` -- thin re-export (or deleted, consumers import geometry directly)
- `rustkernel-geom` -- just `AnalyticalGeomStore` + `GeomAccess` impl
- `rustkernel-solvers` -- deleted (absorbed into geometry::intersection)

## Things to watch for

### Sign convention
The geometry crate uses the same negative-radius convention as the kernel.
No conversion needed. But verify with tests that normals point the right way
after migration.

### NURBS inverse UV
The geometry crate implements 8x8 grid + 5 Gauss-Newton exactly as the kernel
does. The kernel's test suite validates this. If a test fails post-migration,
check that the knot domain is being passed through correctly.

### Serde
Analytical types (Plane, CylinderSurface, etc.) implement Serialize/Deserialize.
`CurveDef` and `SurfaceDef` do NOT derive serde yet because `NurbsCurve3D` from
curvo needs custom serialization. The kernel will need to keep its existing NURBS
serde handling until this is resolved.

### Transform
`CurveDef::apply_transform` and `SurfaceDef::apply_transform` take `&Mat4`.
For NURBS, they pass the matrix directly to curvo's `Transformable::transformed`.
For analytical types, they use `transform_point`/`transform_dir`/`scale_factor`.
This matches the kernel's current behavior.

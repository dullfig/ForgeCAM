# forgecam-geometry

Pure computational geometry library for the ForgeCAM CAM system.

## What this crate is

This is the mathematical geometry foundation — curves, surfaces, intersection algorithms,
meshing, 2D polygon operations, and geometric utilities. It is analogous to SMLib's
NLib + GSNLib layers: everything below trimmed-surface topology.

## What this crate is NOT

No topology. No half-edges, faces, shells, solids, or index-based stores. The B-Rep
kernel (`/ForgeCAM/kernel`) imports this crate and maps between topology indices and
geometry objects. That mapping lives in the kernel, not here.

## Key files

- `API.md` — **Read this first.** Complete specification of every type, trait, and function
  needed. Sections marked "KERNEL REQUIRES" are the immediate priority. Sections marked
  "CAM WILL NEED" can be deferred.
- `src/lib.rs` — Module structure
- `src/types.rs` — Point3, Vec3, Mat4 re-exports
- `src/curves.rs` — Curve structs, CurveDef enum, Curve trait
- `src/surfaces.rs` — Surface structs, SurfaceDef enum, Surface trait
- `src/intersection.rs` — SSI, ray-surface, tri-tri, mesh intersection
- `src/mesh.rs` — TriMesh type
- `src/polygon2d.rs` — 2D polygon point classification and line intersection
- `src/utils.rs` — orthonormal_basis, rodrigues_rotate, AABB, transform helpers

## Dependencies

- `nalgebra 0.34` — linear algebra (must match kernel's version)
- `curvo 0.1.81` — NURBS curve/surface evaluation (must match kernel's version)
- `serde 1` — all public types must be Serialize + Deserialize
- `tracing 0.1` — instrumentation on expensive operations

## Critical conventions

1. **Negative radius = flipped normal.** Cylinder, Sphere, Torus use negative radius
   to indicate inward-facing normals. Cone uses negative half_angle. The kernel depends
   on this sign convention — do not change it.

2. **All angles in radians.**

3. **Tolerances as parameters, not globals.** No global tolerance state.

4. **NURBS inverse UV: 8x8 grid search + 5-iteration Gauss-Newton.** This specific
   algorithm is validated by the kernel's test suite. Use finite-difference Jacobian
   with step 1e-6. Warn (tracing::warn!) if final error > 1e-4.

5. **Parametrizations must match these exactly** (kernel depends on them):
   - Cylinder: u = angle [0, 2pi], v = axial distance
   - Sphere: u = azimuth [0, 2pi], v = elevation [-pi/2, pi/2]
   - Cone: u = angle [0, 2pi], v = distance along axis from apex
   - Torus: u = major angle [0, 2pi], v = minor angle [0, 2pi]
   - LineSegment: t in [0, 1], linear interpolation start→end
   - Circle: t in [0, 1] maps to [0, 2pi]
   - CircularArcCurve: t in [0, 1] maps to [start_angle, end_angle]
   - Ellipse: t in [0, 1] maps to [0, 2pi]

## Implementation priority

Start with "KERNEL REQUIRES" items in API.md. That's everything the kernel uses today.
The "CAM WILL NEED" items are future work for toolpath computation, surface analysis, etc.

## Testing

Every function should have unit tests. For numerical algorithms (inverse UV, intersection,
refinement), test with known analytical solutions and verify tolerances.
The kernel has 278 tests that exercise this geometry indirectly — once migrated, those
tests validate correctness.

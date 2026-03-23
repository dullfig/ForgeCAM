# rustkernel-solvers — Sub-Agent

You own the surface-surface intersection solvers.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `plane_plane.rs` | — | Plane-Plane intersection (line) |
| `plane_sphere.rs` | — | Plane-Sphere intersection (circle) |
| `plane_cylinder.rs` | — | Plane-Cylinder intersection (line or ellipse) |
| `plane_cone.rs` | — | Plane-Cone intersection (conic section) |
| `plane_torus.rs` | — | Plane-Torus intersection (quartic curve, mesh-sampled) |
| `nurbs_solver.rs` | — | PlaneNurbs + NurbsNurbs intersection (mesh-guided) |
| `mesh_intersect.rs` | — | Mesh-mesh intersection utilities (triangle-level) |
| `refine.rs` | — | Newton refinement of intersection points to exact surfaces |
| `lib.rs` | ~39 | Module decls + `default_pipeline()` factory |

## Architecture

Each solver implements the `IntersectionSolver` trait (from `rustkernel_topology::intersection`):
- `can_handle(a: &SurfaceDef, b: &SurfaceDef) -> bool`
- `intersect(a, b, geom, topo) -> Vec<IntersectionCurve>`

`default_pipeline()` registers all solvers in priority order:
1. Analytical solvers first (fast, exact)
2. NURBS solvers as fallback (mesh-guided, slower)

### Analytical solvers (5)
Closed-form solutions for plane vs quadric surfaces. Return exact curve types
(LineSegment, Circle, Ellipse, CircularArc).

### NURBS solvers (2)
Mesh-guided approach: tessellate both surfaces, find mesh-mesh intersections via
`mesh_intersect`, then refine sample points to exact surfaces via `refine.rs`.
Return NURBS curves fit through refined points.

## Conventions

- `tracing::debug_span!` on each solver
- Tolerances as function parameters (no globals)
- Return empty vec (not error) when surfaces don't intersect

## Dependencies

- `rustkernel_topology` — IntersectionPipeline trait, GeomAccess
- `rustkernel_geom` — SurfaceDef, CurveDef
- `rustkernel_math` — Point3, Vec3, tri_tri
- `curvo` — NURBS curve fitting

## Future work

- New solver pairs as new surface types are added (cylinder-cylinder, sphere-cylinder, etc.)
- Performance: BVH for mesh-mesh intersection (Phase 12a)
- These solvers will eventually migrate to `forgecam-geometry::intersection`

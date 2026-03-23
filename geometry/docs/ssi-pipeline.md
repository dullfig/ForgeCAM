# Surface-Surface Intersection Pipeline

The SSI system uses a **solver registry** pattern. Analytical solvers handle
common cases with exact results; a mesh-based fallback handles everything else.

## Architecture

```
SsiPipeline
  |-- PlanePlaneSolver       (line, coincident, or empty)
  |-- PlaneSphereSolver      (circle or empty)
  |-- PlaneCylinderSolver    (circle, ellipse, or empty)
  |-- PlaneConeSolver        (circle, polyline, or empty)
  |-- PlaneTorusSolver       (circle pair or empty)
  |-- NurbsFallbackSolver    (tessellate + mesh intersect + chain + refine)
```

`pipeline.solve(a, b)` iterates solvers in registration order. Each solver
checks `can_solve(kind_a, kind_b)`. The first solver that returns a non-Empty
result wins. If all return Empty, the pipeline returns Empty.

## Usage

```rust
use forgecam_geometry::intersection::*;
use forgecam_geometry::surfaces::*;

// Use the default pipeline (all solvers registered)
let pipeline = default_pipeline();
let result = pipeline.solve(&surface_a, &surface_b);

// Or call analytical functions directly
let result = intersect_plane_sphere(&plane, &sphere);
```

## Result types

```rust
enum SurfaceSurfaceResult {
    Curves(Vec<IntersectionCurve>),  // one or more intersection curves
    Coincident,                       // surfaces are the same
    Empty,                            // no intersection
}

enum IntersectionCurve {
    Line(IntersectionLine),           // infinite line (plane-plane)
    Circle(IntersectionCircle),       // circle (plane-sphere, plane-cyl perp)
    Ellipse(IntersectionEllipse),     // ellipse (plane-cyl oblique)
    Polyline(IntersectionPolyline),   // sampled curve (NURBS, cone oblique)
}
```

## Analytical solvers

### Plane-Plane
Cross normals for direction. Solve 2x2 for a point on the line. Detects
parallel (empty) and coincident cases.

### Plane-Sphere
Signed distance from sphere center to plane. Circle radius from Pythagorean theorem.

### Plane-Cylinder
- **Perpendicular** (normal || axis): circle at the axis-plane intersection
- **Oblique**: ellipse with `semi_minor = |radius|`, `semi_major = |radius| / |dot(normal, axis)|`
- **Parallel** (normal perp axis): returns Empty (line intersections not yet exposed)

### Plane-Cone
- **Perpendicular**: circle at `v * tan(half_angle)` radius
- **Oblique ellipse**: parametric sampling (64 points), returned as Polyline
- **Parabola/hyperbola**: returns Empty (handled by mesh fallback)

### Plane-Torus
- **Axis-perpendicular**: inner and outer circles from meridional cross-section
- **Oblique**: returns Empty (handled by mesh fallback)

## Mesh-based fallback (NurbsFallbackSolver)

Handles any surface pair by brute force:

1. **Tessellate** both surfaces on a uniform grid (default 16x16)
2. **mesh_mesh_intersect**: O(n*m) tri-tri tests, collects raw segments
3. **chain_segments**: greedy nearest-endpoint chaining within tolerance
4. **refine_intersection_points**: alternating projection onto both surfaces

Parameters (stored in NurbsFallbackSolver):
| Parameter | Default | Purpose |
|-----------|---------|---------|
| `tess_divs` | 16 | Grid divisions per surface |
| `chain_tol` | 0.05 | Max gap between chained segments |
| `refine_iters` | 5 | Gauss-Newton iterations per point |
| `refine_tol` | 1e-8 | Convergence threshold |

For plane-vs-NURBS specifically, use `intersect_plane_nurbs()` which uses
`mesh_plane_intersect` (faster than full mesh-mesh).

## Ray-surface intersection

Separate from the SSI pipeline. Used by the kernel's face classifier for
point-in-solid ray casting.

```rust
use forgecam_geometry::intersection::*;

let hit = ray_surface(&ray_origin, &ray_dir, &surface);
if let Some(hit) = hit {
    println!("Hit at t={}, point={:?}, normal={:?}", hit.t, hit.point, hit.normal);
}
```

Supported: Plane, Sphere, Cylinder, Cone. Torus and NURBS return `None`
(torus requires quartic solver; NURBS requires ray-mesh proxy).

## Adding a new solver

```rust
struct MyCylinderCylinderSolver;

impl SsiSolver for MyCylinderCylinderSolver {
    fn can_solve(&self, a: SurfaceKind, b: SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Cylinder, SurfaceKind::Cylinder)
        )
    }

    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        // Your analytical cylinder-cylinder intersection here
        todo!()
    }
}

// Register before the fallback solver
let mut pipeline = SsiPipeline::new();
pipeline.register(Box::new(MyCylinderCylinderSolver));
pipeline.register(Box::new(NurbsFallbackSolver::default()));
```

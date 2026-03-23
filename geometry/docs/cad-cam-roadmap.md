# Geometry Operations Needed for CAD/CAM

Everything below is specified in API.md under "CAM WILL NEED" but not yet implemented.
Organized by what drives the need, roughly prioritized.

---

## Priority 1: 2D Curve Operations (Toolpath Foundation)

Nearly all toolpath computation reduces to 2D planar geometry — offset a contour,
intersect the result with stock/boundary, chain the output. These are the workhorses.

### Curve-Curve Intersection (2D)

```rust
fn intersect_lines_2d(a: &LineSegment2D, b: &LineSegment2D) -> Option<Point2>;
fn intersect_line_arc_2d(line: &LineSegment2D, arc: &Arc2D) -> Vec<Point2>;
fn intersect_arc_arc_2d(a: &Arc2D, b: &Arc2D) -> Vec<Point2>;
```

**Who needs it:** Cutter comp, 2D pocketing, contour clipping, stock boundary trimming.
`cavalier_contours` handles polyline-level offset+boolean, but raw geometric intersection
is still needed for profile construction and custom contour ops.

### Curve Offset

```rust
fn offset(curve: &CurveDef, d: f64, axis: &Vec3) -> CurveDef;
```

**Who needs it:** Cutter radius compensation, stock allowance, finish pass offset.
For composite curves (profile boundaries), corners need filleting (convex) or
extension+intersection (concave). This is the single most important CAM geometry op
after what's already implemented.

Note: `cavalier_contours` does polyline offset with native arcs. This function is
for individual analytical/NURBS curve offset — different layer.

### Curve Split / Trim / Extend / Reverse

```rust
fn split(curve: &CurveDef, t: f64) -> (CurveDef, CurveDef);
fn trim(curve: &CurveDef, t0: f64, t1: f64) -> CurveDef;
fn extend(curve: &CurveDef, d: f64, end: CurveEnd) -> CurveDef;
fn reverse(curve: &mut CurveDef);
```

**Who needs it:** Toolpath entry/exit moves, lead-in arcs, contour chaining after
intersection, trimming curves to stock boundaries. Every toolpath operation that
isn't a complete closed contour needs trim/split.

---

## Priority 2: Curve Query & Measurement

### Arc Length + Parameterization

```rust
fn length(curve: &CurveDef) -> f64;
fn length_between(curve: &CurveDef, t0: f64, t1: f64) -> f64;
fn parameter_at_length(curve: &CurveDef, t_start: f64, length: f64) -> f64;
```

**Who needs it:** Constant-feedrate motion along curves (G-code generation needs
equal-arc-length spacing, not equal-parameter spacing). Also needed for scallop
calculation and stepover control on 3D surface passes.

Analytical curves have closed-form arc length (line = distance, circle = r*theta,
ellipse = elliptic integral). NURBS need Gauss-Legendre quadrature.

### Closest Point

```rust
fn closest_point(curve: &CurveDef, point: &Point3) -> (f64, Point3, f64);
```

**Who needs it:** Gouge checking (is the tool closer to the part than the cutter
radius?), snap-to-curve in GUI, toolpath validation.

### Curvature

```rust
fn curvature(curve: &CurveDef, t: f64) -> f64;
fn inflection_points(curve: &CurveDef) -> Vec<f64>;
```

**Who needs it:** Feed rate limiting in corners (high curvature = slow down),
adaptive tessellation, feature recognition (sharp vs smooth transitions).

### Tessellation (Tolerance-Based)

```rust
fn tessellate_by_tolerance(curve: &CurveDef, chord_tol: f64) -> Vec<Point3>;
fn tessellate_by_angle(curve: &CurveDef, angle_tol: f64) -> Vec<Point3>;
```

**Who needs it:** G-code linearization — curved toolpaths need to be broken into
G01 moves that stay within a chord tolerance. Current `sample(n)` is uniform in
parameter space, not adaptive to curvature.

### Other Curve Queries

```rust
fn project_to_plane(curve: &CurveDef, origin: &Point3, normal: &Vec3) -> CurveDef;
fn tangent_parallel_to(curve: &CurveDef, direction: &Vec3) -> Vec<f64>;
fn bounding_box(curve: &CurveDef) -> Aabb;
```

**Who needs it:** `project_to_plane` for flattening 3D curves to workplane.
`tangent_parallel_to` for finding silhouette points, tangent contacts.
`bounding_box` for spatial queries and broad-phase filtering.

---

## Priority 3: Surface Analysis (3D Machining)

### Surface Offset

```rust
fn offset(surface: &SurfaceDef, d: f64) -> SurfaceDef;
```

**Who needs it:** Already partially in the kernel (`offset_solid_geometry`), but
the pure geometry op belongs here. Cutter-contact to cutter-location surface
conversion (offset by tool radius). Shell/hollow. Stock envelope.

Analytical surfaces have exact offsets. NURBS need sample + refit (approximate).

Note: the kernel already has `SurfaceDef::offset()` returning `Option<SurfaceDef>` —
this may already be implemented. Check before duplicating.

### Surface Closest Point

```rust
fn closest_point(surface: &SurfaceDef, point: &Point3) -> (f64, f64, Point3, f64);
```

**Who needs it:** 3D gouge checking, projection of toolpath onto part surface,
rest material detection. Core operation for any 3-axis or 5-axis finishing.

For analytical surfaces, this is a constrained minimization with closed-form
solutions. For NURBS, Newton iteration from grid-search seed (similar to
`inverse_uv` but minimizing distance, not just finding parameter).

### Surface Curvature

```rust
fn curvature(surface: &SurfaceDef, u: f64, v: f64) -> CurvatureInfo;
fn partials(surface: &SurfaceDef, u: f64, v: f64) -> (Vec3, Vec3);
```

**Who needs it:** Scallop height calculation (stepover vs surface finish),
adaptive step size in contouring, feature recognition (flat vs curved regions).

`CurvatureInfo` = `{ k_min, k_max, dir_min, dir_max }` (principal curvatures).

### Scallop Height

```rust
fn scallop_height(surface: &SurfaceDef, u: f64, v: f64, step_over: f64, tool_radius: f64) -> f64;
```

**Who needs it:** Surface finish prediction. Given a ball-end mill of radius R
at stepover S on a surface with curvature k, the scallop height tells you whether
the finish meets the print callout. Drives adaptive stepover in finishing passes.

### Iso-Curves & Projection

```rust
fn iso_curve(surface: &SurfaceDef, param: IsoCurveParam) -> CurveDef;
fn project_curve(surface: &SurfaceDef, curve: &CurveDef) -> Vec<Point2>;
fn bounding_box(surface: &SurfaceDef) -> Aabb;
```

**Who needs it:** `iso_curve` for UV-based surface machining strategies (zig-zag,
spiral). `project_curve` for projecting boundary curves into UV space (trimmed
surface toolpath generation). `bounding_box` for spatial queries.

---

## Priority 4: Ray Casting (Simulation & Verification)

```rust
fn ray_torus(origin: &Point3, dir: &Vec3, torus: &TorusSurface) -> Vec<RayHit>;
fn ray_nurbs_mesh(origin: &Point3, dir: &Vec3, nurbs: &FlippableNurbs, mesh: &TriMesh) -> Option<RayHit>;
```

**Who needs it:** Toolpath simulation (dexel/Z-buffer material removal), collision
detection, visibility testing. `ray_torus` is a quartic — up to 4 intersections.
`ray_nurbs_mesh` uses the tessellation proxy.

Currently `ray_surface()` returns `None` for torus and NURBS. These fill the gaps.

---

## Priority 5: 3D Curve-Curve Intersection

```rust
fn intersect_curves(a: &CurveDef, b: &CurveDef, tol: f64) -> Vec<(f64, f64, Point3)>;
```

**Who needs it:** Trimmed surface boundary computation, STEP import repair,
general 3D modeling operations. Less urgent than 2D intersection since most
toolpath work projects to a plane first.

---

## Priority 6: Mesh Operations

```rust
fn compute_vertex_normals(mesh: &mut TriMesh);
fn mesh_bounding_box(mesh: &TriMesh) -> Aabb;
fn decimate(mesh: &TriMesh, target_count: usize) -> TriMesh;
fn mesh_boolean(a: &TriMesh, b: &TriMesh, op: MeshBooleanOp) -> TriMesh;
```

**Who needs it:** `compute_vertex_normals` for smooth shading in viewer/GUI.
`mesh_bounding_box` for spatial queries. `decimate` for LOD in viewer.
`mesh_boolean` for toolpath simulation (stock material removal via mesh booleans
rather than dexel — alternative approach).

These are lower priority because the kernel already handles most mesh operations
through tessellation, and the viewer can compute normals at render time.

---

## Priority 7: Constant-Engagement Roughing (FACEOM) — 3-Axis and 5-Axis

### Already Implemented (geometry crate)

The core geometric primitive is done: `engagement_angle()` in `engagement.rs`.
Given a tool circle (center + radius) and a 2D stock boundary (`Contour2D`),
it computes the total engagement angle and the individual engaged arcs. Works
with line segments and arcs. Handles full slotting, zero engagement, straight
walls, circular stock, and concave corners.

Helper functions: `engagement_from_radial_depth()`, `radial_depth_from_engagement()`.

### Needed in toolpath crate: FACEOM Path Planner (3-axis)

```rust
/// Generate a constant-engagement roughing toolpath for a 2D pocket.
fn faceom_roughing(
    pocket: &Contour2D,           // pocket boundary
    stock: &Contour2D,            // current stock boundary (may differ from pocket)
    tool_radius: f64,
    target_engagement: f64,       // radians, typically 0.7–1.05 (40°–60°)
    axial_depth: f64,             // Z step
) -> Vec<Contour2D>;             // toolpath contours, innermost first
```

**Algorithm (Lukacs et al. 2019):**
1. Start from innermost offset contour (or seed point).
2. At each toolpath point, call `engagement_angle()` against the current stock boundary.
3. Solve for the offset distance that achieves `target_engagement`.
4. Step along, adjusting offset distance per-step → non-equidistant spiral.
5. Optionally smooth with cubic splines for C2 continuity (machine dynamics).
6. After each pass, update the stock boundary (boolean subtract the swept tool).

This lives in the **toolpath crate**, not geometry. Geometry provides the engagement
primitive; toolpath owns the path planning loop, stock updates, entry/exit moves,
lead-in arcs, and G-code output.

### Needed in toolpath crate: 5-Axis Constant-Engagement

5-axis adaptive clearing is the same 2D engagement problem in a tilted plane.
No one does full FACEOM in 5 axes — this is novel.

**Pipeline:**

1. **Define the cutting plane.** At each toolpath point, the cutting plane is
   perpendicular to the current tool axis. The tool axis comes from the drive
   surface normal (or a lead/lag/tilt strategy).

2. **Project stock boundary.** Slice or project the 3D stock body into the
   tool-perpendicular plane → produces a `Contour2D` in local 2D coordinates.

3. **Compute engagement.** Call `engagement_angle()` in the local 2D plane.
   This is the exact same function used for 3-axis.

4. **Adjust stepover.** FACEOM logic: vary the offset distance so engagement
   stays at the target. Step along the surface, re-projecting at each point.

5. **Orient the tool.** Apply lead/lag/tilt angles for chip thinning, tool life,
   and surface finish. Common strategies:
   - **Surface normal** — tool axis = part surface normal at contact point
   - **Lead/lag** — tilt forward/back along feed direction (typically 3–5°)
   - **Sideways tilt** — lean into the cut for better engagement with ball-nose

6. **Update 3D stock.** After each pass, remove the swept volume of the tilted
   tool from the stock model. This is the hard part — a tilted endmill sweeps
   a different volume than a vertical one. Options:
   - Dexel model (Z-buffer per axis) — fast, approximate
   - Tri-dexel (3 orthogonal Z-buffers) — better for 5-axis
   - Mesh boolean — exact but slow

**Geometry primitives needed (not yet implemented):**

```rust
/// Project a 3D contour (trimmed surface boundary, stock slice) onto a plane,
/// returning a 2D contour in the plane's local coordinate system.
fn project_contour_to_plane(
    points: &[Point3],
    plane_origin: Point3,
    plane_normal: Vec3,
) -> Contour2D;

/// Slice a TriMesh with a plane, returning the cross-section as 2D contour(s)
/// in the plane's local coordinate system.
fn slice_mesh_with_plane(
    mesh: &TriMesh,
    plane_origin: Point3,
    plane_normal: Vec3,
) -> Vec<Contour2D>;
```

The `slice_mesh_with_plane` function is the key enabler: it produces the 2D
stock boundary in any tool-perpendicular plane, which feeds directly into
`engagement_angle()`. The mesh-plane intersection already exists in
`intersection.rs` — this would chain those segments into closed contours.

**Why this matters:** Constant-engagement roughing dramatically improves tool
life and MRR in 3-axis. Extending it to 5-axis (where engagement control is
currently manual or nonexistent in most CAM systems) would be a significant
competitive advantage, especially for impeller/blisk roughing and deep-pocket
aerospace parts where 5-axis access is required.

---

## What's NOT Listed Here

These are already implemented and working:

- All curve/surface struct definitions and enums
- Curve eval, tangent, domain, sample, apply_transform, translate
- Surface eval, normal, domain, inverse_uv, flip_normal, apply_transform, translate, tessellate
- All analytical SSI (plane-plane/sphere/cylinder/cone/torus)
- NURBS SSI pipeline (tessellate + mesh intersect + chain + refine)
- Ray-plane/sphere/cylinder/cone
- Triangle-triangle intersection (Moller)
- Mesh-plane and mesh-mesh intersection + segment chaining
- SSI solver pipeline + default_pipeline()
- TriMesh construction + grid tessellation
- Polygon2D (signed_area, classify_point, intersect_line)
- All utilities (orthonormal_basis, rodrigues, arc generation, AABB, transform helpers)
- Cutter engagement angle computation (engagement.rs)
- EntityRef persistent entity IDs (entity_ref.rs)

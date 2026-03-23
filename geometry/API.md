# forgecam-geometry API Specification

Pure computational geometry library for ForgeCAM. No topology, no B-Rep, no half-edges.
This crate is imported by the kernel, the CAM engine, and any other ForgeCAM component
that needs geometric computation.

**Rule**: nothing in this crate may depend on `rustkernel-topology` or any kernel crate.
Dependency flows one way: `kernel → geometry`, never the reverse.

**Dependencies**: `nalgebra 0.34`, `curvo 0.1.81`, `serde 1`, `tracing 0.1`.

---

## 1. Types (`types.rs`)

Re-export core numeric types so consumers don't need a direct nalgebra dependency.

```rust
pub type Point2 = nalgebra::Point2<f64>;
pub type Point3 = nalgebra::Point3<f64>;
pub type Vec2 = nalgebra::Vector2<f64>;
pub type Vec3 = nalgebra::Vector3<f64>;
pub type Mat4 = nalgebra::Matrix4<f64>;

/// Re-export nalgebra itself for consumers that need it
pub use nalgebra;
```

All types must derive `Debug, Clone`. Where possible, also `Serialize, Deserialize`.

### GPU helpers

```rust
pub fn point3_to_f32(p: &Point3) -> [f32; 3];
pub fn vec3_to_f32(v: &Vec3) -> [f32; 3];
```

---

## 2. Curves (`curves.rs`)

### 2.1 Curve Structs

All derive `Debug, Clone, Serialize, Deserialize`.

```rust
pub struct LineSegment {
    pub start: Point3,
    pub end: Point3,
}

pub struct CircleCurve {
    pub center: Point3,
    pub axis: Vec3,      // unit normal to circle plane
    pub radius: f64,
    pub ref_dir: Vec3,   // unit vector in circle plane, defines angle=0
}

pub struct CircularArcCurve {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub ref_dir: Vec3,
    pub start_angle: f64,  // radians
    pub end_angle: f64,    // radians
}

pub struct EllipseCurve {
    pub center: Point3,
    pub axis: Vec3,         // unit normal to ellipse plane
    pub semi_major: f64,
    pub semi_minor: f64,
    pub major_dir: Vec3,    // unit vector along semi-major axis
}
```

NURBS curves use `curvo::NurbsCurve3D<f64>` directly (re-exported).

### 2.2 CurveDef Enum

```rust
pub enum CurveDef {
    LineSegment(LineSegment),
    Circle(CircleCurve),
    CircularArcCurve(CircularArcCurve),
    Ellipse(EllipseCurve),
    Nurbs(NurbsCurve3D<f64>),
}
```

### 2.3 Curve Trait

Implemented for each concrete curve struct AND for `CurveDef` (dispatches to variant).

```rust
pub trait Curve {
    /// Evaluate position at parameter t.
    /// For analytical curves, t is in [0, 1] (normalized).
    /// For NURBS, t is in knot domain.
    fn eval(&self, t: f64) -> Point3;

    /// Unit tangent vector at parameter t.
    fn tangent(&self, t: f64) -> Vec3;

    /// Parameter domain: (t_min, t_max).
    fn domain(&self) -> (f64, f64);
}
```

### 2.4 Curve Operations (free functions or methods on CurveDef)

**KERNEL REQUIRES** (used today):

```rust
/// Uniform parameter sampling — returns n+1 points from domain start to end.
/// Used by boolean curve_trimming (circle/ellipse → 32-chord polyline).
fn sample(&self, n: usize) -> Vec<Point3>;

/// Apply affine transform (rigid + uniform scale).
/// LineSegment: transform endpoints.
/// Circle/Arc: transform center, normalize axis & ref_dir, scale radius.
/// Ellipse: transform center, normalize axis & major_dir, scale semi-axes.
/// NURBS: delegate to curvo Transformable.
fn apply_transform(&mut self, m: &Mat4);

/// Translate by offset vector.
fn translate(&mut self, delta: &Vec3);
```

**CAM WILL NEED** (from GSNLib feature list — implement when needed):

```rust
/// Offset curve by distance d in the plane perpendicular to `axis`.
/// Returns a new curve. For composite curves, fillet or extend at corners.
fn offset(&self, d: f64, axis: &Vec3) -> CurveDef;

/// Project curve onto a plane, returning a new 3D curve lying on the plane.
fn project_to_plane(&self, plane_origin: &Point3, plane_normal: &Vec3) -> CurveDef;

/// Closest point on curve to a given point. Returns (t, point, distance).
fn closest_point(&self, point: &Point3) -> (f64, Point3, f64);

/// Arc length of entire curve (high precision for NURBS).
fn length(&self) -> f64;

/// Arc length between two parameter values.
fn length_between(&self, t0: f64, t1: f64) -> f64;

/// Parameter at a given arc length from t_start. For constant-speed traversal.
fn parameter_at_length(&self, t_start: f64, length: f64) -> f64;

/// Curvature (1/radius) at parameter t.
fn curvature(&self, t: f64) -> f64;

/// Find inflection points (where curvature changes sign). Returns parameter values.
fn inflection_points(&self) -> Vec<f64>;

/// Split curve at parameter t. Returns (left, right).
fn split(&self, t: f64) -> (CurveDef, CurveDef);

/// Trim curve to parameter range [t0, t1].
fn trim(&self, t0: f64, t1: f64) -> CurveDef;

/// Extend curve beyond its domain by distance d at start (negative) or end (positive).
fn extend(&self, d: f64, end: CurveEnd) -> CurveDef;

/// Reverse curve direction (swap start and end).
fn reverse(&mut self);

/// Tessellate to polyline with chord-height tolerance.
fn tessellate_by_tolerance(&self, chord_tol: f64) -> Vec<Point3>;

/// Tessellate to polyline with angular tolerance (max angle between successive chords).
fn tessellate_by_angle(&self, angle_tol: f64) -> Vec<Point3>;

/// Find parameters where tangent is parallel to a given vector.
fn tangent_parallel_to(&self, direction: &Vec3) -> Vec<f64>;

/// Bounding box.
fn bounding_box(&self) -> Aabb;
```

---

## 3. Surfaces (`surfaces.rs`)

### 3.1 Surface Structs

All derive `Debug, Clone, Serialize, Deserialize`.

```rust
pub struct Plane {
    pub origin: Point3,
    pub normal: Vec3,     // unit normal
}

pub struct CylinderSurface {
    pub origin: Point3,   // point on axis
    pub axis: Vec3,       // unit direction along axis
    pub radius: f64,      // negative = inward-facing normal
}

pub struct SphereSurface {
    pub center: Point3,
    pub radius: f64,      // negative = inward-facing normal
}

pub struct ConeSurface {
    pub apex: Point3,
    pub axis: Vec3,       // unit direction from apex toward base
    pub half_angle: f64,  // radians; negative = inward-facing normal
}

pub struct TorusSurface {
    pub center: Point3,
    pub axis: Vec3,       // unit axis of revolution
    pub major_radius: f64,
    pub minor_radius: f64, // negative = inward-facing normal
}

/// NURBS surface wrapper with orientation flag.
pub struct FlippableNurbs {
    pub surface: NurbsSurface3D<f64>,
    pub flipped: bool,
}
```

### 3.2 SurfaceDef Enum

```rust
pub enum SurfaceDef {
    Plane(Plane),
    Cylinder(CylinderSurface),
    Sphere(SphereSurface),
    Cone(ConeSurface),
    Torus(TorusSurface),
    Nurbs(FlippableNurbs),
}
```

### 3.3 Surface Trait

Implemented for each concrete surface struct AND for `SurfaceDef` (dispatches to variant).

```rust
pub trait Surface {
    /// Evaluate position at parameters (u, v).
    /// Parametrizations:
    ///   Plane: returns origin (u,v unused in current kernel — may change)
    ///   Cylinder: u=angle [0, 2pi], v=axial distance
    ///   Sphere: u=azimuth [0, 2pi], v=elevation [-pi/2, pi/2]
    ///   Cone: u=angle [0, 2pi], v=distance along axis from apex
    ///   Torus: u=major angle [0, 2pi], v=minor angle [0, 2pi]
    ///   NURBS: u,v in knot domain
    fn eval(&self, u: f64, v: f64) -> Point3;

    /// Unit outward normal at (u, v).
    /// Sign convention: negative radius/half_angle flips normal direction.
    /// NURBS: respects FlippableNurbs.flipped flag.
    fn normal(&self, u: f64, v: f64) -> Vec3;

    /// Parameter domain: ((u_min, u_max), (v_min, v_max)).
    /// Analytical: ((0, 2pi), (surface-dependent)).
    /// NURBS: from knot vector.
    fn domain(&self) -> ((f64, f64), (f64, f64));

    /// Inverse parametric mapping: given a 3D point near the surface,
    /// find the closest (u, v) parameters.
    /// Analytical surfaces: closed-form (atan2, asin, etc.).
    /// NURBS: 8x8 grid search + Gauss-Newton refinement (5 iterations).
    fn inverse_uv(&self, point: &Point3) -> (f64, f64);
}
```

### 3.4 Surface Operations (free functions or methods on SurfaceDef)

**KERNEL REQUIRES** (used today):

```rust
/// Flip the surface normal direction.
/// Plane: negate normal. Cylinder/Sphere/Torus: negate radius.
/// Cone: negate half_angle. NURBS: toggle flipped flag.
fn flip_normal(&mut self);

/// Apply affine transform (rigid + uniform scale).
/// Plane: transform origin, direction-normalize normal.
/// Cylinder: transform origin & axis, scale radius by scale_factor(m).
/// Sphere: transform center, scale radius.
/// Cone: transform apex & axis (half_angle unchanged).
/// Torus: transform center & axis, scale both radii.
/// NURBS: delegate to curvo Transformable.
fn apply_transform(&mut self, m: &Mat4);

/// Translate surface by offset vector.
fn translate(&mut self, delta: &Vec3);

/// Generate a triangle mesh of the surface over its domain.
/// Returns None for analytical surfaces that don't need standalone tessellation
/// (their tessellation is driven by face boundaries in the kernel).
/// Returns Some(TriMesh) for NURBS: (divs_u+1)*(divs_v+1) vertices,
/// divs_u*divs_v*2 triangles with positions, normals, and UVs.
fn tessellate(&self, divs_u: usize, divs_v: usize) -> Option<TriMesh>;
```

**CAM WILL NEED** (implement when needed):

```rust
/// Offset surface by distance d along normal direction.
/// Analytical surfaces have closed-form offsets:
///   Plane(n) → Plane(origin + d*n, n)
///   Cylinder(r) → Cylinder(r + d)
///   Sphere(r) → Sphere(r + d)
///   Cone(apex, axis, alpha) → Cone(apex + d/sin(alpha)*axis, axis, alpha)
///   Torus(R, r) → Torus(R, r + d)
///   NURBS → sample S(u,v) + d*N(u,v), refit NURBS (approximate)
fn offset(&self, d: f64) -> SurfaceDef;

/// Curvature at (u, v). Returns (k_min, k_max, dir_min, dir_max).
/// Principal curvatures and their directions.
fn curvature(&self, u: f64, v: f64) -> CurvatureInfo;

/// Partial derivatives at (u, v): dS/du and dS/dv.
fn partials(&self, u: f64, v: f64) -> (Vec3, Vec3);

/// Closest point on surface to a given point.
/// Returns (u, v, point_on_surface, distance).
fn closest_point(&self, point: &Point3) -> (f64, f64, Point3, f64);

/// Project a 3D curve onto the surface.
/// Returns a 2D curve in (u, v) parameter space.
fn project_curve(&self, curve: &CurveDef) -> Vec<Point2>;

/// Iso-parametric curve at constant u or constant v.
fn iso_curve(&self, param: IsoCurveParam) -> CurveDef;

/// Bounding box (tight for analytical, approximate for NURBS).
fn bounding_box(&self) -> Aabb;

/// Scallop height between adjacent toolpaths at distance `step_over` apart.
/// Uses local curvature. Essential for CAM surface finish calculation.
fn scallop_height(&self, u: f64, v: f64, step_over: f64, tool_radius: f64) -> f64;
```

---

## 4. Intersection (`intersection.rs`)

May be split into submodules: `intersection/ssi.rs`, `intersection/ray.rs`,
`intersection/tri_tri.rs`, `intersection/mesh.rs`, `intersection/refine.rs`.

### 4.1 Result Types

```rust
pub enum SurfaceSurfaceResult {
    Curves(Vec<IntersectionCurve>),
    Coincident,
    Empty,
}

pub enum IntersectionCurve {
    Line(IntersectionLine),
    Circle(IntersectionCircle),
    Ellipse(IntersectionEllipse),
    Polyline(IntersectionPolyline),
}

pub struct IntersectionLine {
    pub origin: Point3,
    pub direction: Vec3,  // unit direction
}

pub struct IntersectionCircle {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub ref_dir: Vec3,
}

pub struct IntersectionEllipse {
    pub center: Point3,
    pub axis: Vec3,
    pub semi_major: f64,
    pub semi_minor: f64,
    pub major_dir: Vec3,
}

pub struct IntersectionPolyline {
    pub points: Vec<Point3>,
}

pub struct RayHit {
    pub t: f64,           // parameter along ray
    pub point: Point3,    // hit point
    pub normal: Vec3,     // surface normal at hit
}
```

### 4.2 Surface-Surface Intersection (analytical)

**KERNEL REQUIRES** (used today in rustkernel-solvers):

```rust
/// Plane-Plane: cross product of normals → line direction.
/// Returns Line, Coincident, or Empty.
pub fn intersect_plane_plane(a: &Plane, b: &Plane) -> SurfaceSurfaceResult;

/// Plane-Sphere: signed distance to plane, Pythagorean circle radius.
/// Returns Circle or Empty.
pub fn intersect_plane_sphere(p: &Plane, s: &SphereSurface) -> SurfaceSurfaceResult;

/// Plane-Cylinder: angle analysis → Circle (perpendicular) or Ellipse (oblique).
/// Returns Circle, Ellipse, or Empty.
pub fn intersect_plane_cylinder(p: &Plane, c: &CylinderSurface) -> SurfaceSurfaceResult;

/// Plane-Cone: angle analysis → Circle or Ellipse.
/// Returns Circle, Ellipse, or Empty.
pub fn intersect_plane_cone(p: &Plane, c: &ConeSurface) -> SurfaceSurfaceResult;

/// Plane-Torus: meridional analysis → pair of Circles.
/// Returns Circle pair or Empty.
pub fn intersect_plane_torus(p: &Plane, t: &TorusSurface) -> SurfaceSurfaceResult;
```

### 4.3 NURBS / General SSI (mesh-based)

**KERNEL REQUIRES** (used today):

```rust
/// Plane vs NURBS surface intersection.
/// Algorithm: tessellate NURBS (16x16) → mesh_plane_intersect → chain → refine.
pub fn intersect_plane_nurbs(
    p: &Plane,
    ns: &FlippableNurbs,
    tess_divs: usize,     // default 16
    chain_tol: f64,        // default 0.05
    refine_iters: usize,   // default 5
    refine_tol: f64,        // default 1e-8
) -> SurfaceSurfaceResult;

/// General surface-surface intersection (at least one NURBS, or any pair).
/// Algorithm: tessellate both → mesh_mesh_intersect → chain → refine.
pub fn intersect_surfaces(
    a: &SurfaceDef,
    b: &SurfaceDef,
    tess_divs: usize,
    chain_tol: f64,
    refine_iters: usize,
    refine_tol: f64,
) -> SurfaceSurfaceResult;

/// Refine approximate intersection points to lie on both surfaces.
/// Alternating projection: project onto A → B → midpoint, repeat.
pub fn refine_intersection_points(
    a: &dyn Surface,
    b: &dyn Surface,
    points: &[Point3],
    max_iter: usize,
    tol: f64,
) -> Vec<Point3>;
```

### 4.4 SSI Solver Pipeline

**KERNEL REQUIRES** (used today — the kernel registers solvers and dispatches):

```rust
/// A solver that can handle a specific pair of surface types.
pub trait SsiSolver {
    /// Returns true if this solver can handle (kind_a, kind_b) in either order.
    fn can_solve(&self, kind_a: SurfaceKind, kind_b: SurfaceKind) -> bool;

    /// Solve using concrete surface definitions.
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult;
}

/// Discriminant enum for surface type matching (no geometric data).
pub enum SurfaceKind {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Nurbs,
    Unknown,
}

/// Similarly for curves.
pub enum CurveKind {
    LineSegment,
    Circle,
    CircularArcCurve,
    Ellipse,
    Nurbs,
    Unknown,
}

/// Registry of solvers tried in order. Returns first non-Empty result.
pub struct SsiPipeline {
    solvers: Vec<Box<dyn SsiSolver>>,
}

impl SsiPipeline {
    pub fn new() -> Self;
    pub fn register(&mut self, solver: Box<dyn SsiSolver>);
    pub fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult;
}

/// Pre-built pipeline with all analytical + NURBS solvers registered.
pub fn default_pipeline() -> SsiPipeline;
```

### 4.5 Ray-Surface Intersection

**KERNEL REQUIRES** (used today in face_classifier.rs for point-in-solid ray casting):

```rust
/// Ray-Plane intersection.
/// denom = ray_dir . normal; t = (origin - point) . normal / denom.
pub fn ray_plane(ray_origin: &Point3, ray_dir: &Vec3, plane: &Plane) -> Option<RayHit>;

/// Ray-Sphere intersection.
/// Quadratic: |P + t*d - C|^2 = r^2. Returns nearest positive hit.
pub fn ray_sphere(ray_origin: &Point3, ray_dir: &Vec3, sphere: &SphereSurface) -> Option<RayHit>;

/// Ray-Cylinder intersection (infinite cylinder).
/// Projects perpendicular to axis, solves 2D quadratic.
pub fn ray_cylinder(ray_origin: &Point3, ray_dir: &Vec3, cyl: &CylinderSurface) -> Option<RayHit>;

/// Ray-Cone intersection (infinite cone from apex).
/// Quadratic from cone equation.
pub fn ray_cone(ray_origin: &Point3, ray_dir: &Vec3, cone: &ConeSurface) -> Option<RayHit>;

/// Ray vs any SurfaceDef — dispatches to appropriate function above.
/// NURBS: not yet supported (returns None); add ray-mesh when needed.
pub fn ray_surface(ray_origin: &Point3, ray_dir: &Vec3, surface: &SurfaceDef) -> Option<RayHit>;
```

**CAM WILL NEED** (implement when needed):

```rust
/// Ray-Torus intersection (quartic — 0-4 solutions).
pub fn ray_torus(ray_origin: &Point3, ray_dir: &Vec3, torus: &TorusSurface) -> Vec<RayHit>;

/// Ray-NURBS intersection via mesh proxy.
pub fn ray_nurbs_mesh(ray_origin: &Point3, ray_dir: &Vec3, nurbs: &FlippableNurbs, mesh: &TriMesh) -> Option<RayHit>;
```

### 4.6 Triangle-Triangle Intersection

**KERNEL REQUIRES** (used today in mesh_mesh_intersect):

```rust
pub struct TriTriResult {
    pub intersects: bool,
    pub segment: Option<(Point3, Point3)>,
}

/// Moller's triangle-triangle intersection algorithm.
/// Returns intersection segment endpoints, or None.
/// Handles coplanar, edge-touching, disjoint cases.
/// Snaps near-zero signed distances (1e-12), uses best-magnitude component for projection.
pub fn tri_tri_intersect(
    a0: &Point3, a1: &Point3, a2: &Point3,
    b0: &Point3, b1: &Point3, b2: &Point3,
) -> TriTriResult;
```

### 4.7 Mesh-Level Intersection

**KERNEL REQUIRES** (used today in NURBS SSI pipeline):

```rust
/// Triangle mesh vs plane intersection.
/// Per-triangle: classify vertices by signed distance to plane,
/// find edge crossings via linear interpolation.
/// Returns raw segments.
pub fn mesh_plane_intersect(mesh: &TriMesh, plane: &Plane) -> Vec<(Point3, Point3)>;

/// Triangle mesh vs triangle mesh intersection.
/// O(n*m) brute force — calls tri_tri_intersect per pair.
/// Returns raw segments.
pub fn mesh_mesh_intersect(a: &TriMesh, b: &TriMesh) -> Vec<(Point3, Point3)>;

/// Chain raw intersection segments into continuous polylines.
/// Greedy nearest-endpoint chaining within tolerance.
pub fn chain_segments(segments: Vec<(Point3, Point3)>, tolerance: f64) -> Vec<Vec<Point3>>;
```

### 4.8 Curve-Curve Intersection

**CAM WILL NEED** (implement when needed):

```rust
/// 2D line-line intersection.
pub fn intersect_lines_2d(a: &LineSegment2D, b: &LineSegment2D) -> Option<Point2>;

/// 2D line-arc intersection.
pub fn intersect_line_arc_2d(line: &LineSegment2D, arc: &Arc2D) -> Vec<Point2>;

/// 2D arc-arc intersection.
pub fn intersect_arc_arc_2d(a: &Arc2D, b: &Arc2D) -> Vec<Point2>;

/// General 3D curve-curve closest approach / intersection.
pub fn intersect_curves(a: &CurveDef, b: &CurveDef, tol: f64) -> Vec<(f64, f64, Point3)>;
```

---

## 5. Mesh (`mesh.rs`)

### 5.1 TriMesh Type

**KERNEL REQUIRES** (used today as FaceMesh in topology, but it's pure geometry):

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriMesh {
    pub positions: Vec<Point3>,
    pub normals: Vec<Vec3>,
    pub indices: Vec<[u32; 3]>,  // triangle index triples
    pub uvs: Vec<(f64, f64)>,    // per-vertex UV coordinates
}
```

### 5.2 Mesh Construction

**KERNEL REQUIRES** (used today for NURBS tessellation and SSI):

```rust
/// Tessellate a surface over a uniform parameter grid.
/// Used by NURBS SSI solvers to create mesh proxies for intersection.
/// Evaluates (divs_u+1)*(divs_v+1) grid points via surface.eval() and surface.normal().
/// Generates 2 triangles per quad cell.
pub fn tessellate_surface_grid(
    surface: &dyn Surface,
    divs_u: usize,
    divs_v: usize,
) -> TriMesh;
```

**CAM WILL NEED** (implement when needed):

```rust
/// Compute vertex normals by averaging adjacent face normals.
pub fn compute_vertex_normals(mesh: &mut TriMesh);

/// Compute AABB for mesh.
pub fn mesh_bounding_box(mesh: &TriMesh) -> Aabb;

/// Decimate mesh to target triangle count.
pub fn decimate(mesh: &TriMesh, target_count: usize) -> TriMesh;

/// Mesh boolean (union/difference/intersection) on triangle meshes.
pub fn mesh_boolean(a: &TriMesh, b: &TriMesh, op: MeshBooleanOp) -> TriMesh;
```

---

## 6. 2D Polygon (`polygon2d.rs`)

**KERNEL REQUIRES** (used today in boolean face classifier and curve trimming):

```rust
pub type Point2D = [f64; 2];  // lightweight, no nalgebra overhead

pub struct LineHit {
    pub t_line: f64,      // parameter along the infinite line
    pub edge_index: usize, // which polygon edge was hit (i → i+1)
    pub t_edge: f64,       // parameter along edge [0, 1]
}

pub enum PointClassification {
    Inside,
    Outside,
    OnBoundary,
}

#[derive(Debug, Clone)]
pub struct Polygon2D {
    pub vertices: Vec<Point2D>,
}

impl Polygon2D {
    /// Shoelace formula. Positive = CCW winding, negative = CW.
    pub fn signed_area(&self) -> f64;

    /// Point-in-polygon classification via crossing number + boundary snap.
    /// Checks all edges for boundary proximity (within tol), then ray casts in +x.
    pub fn classify_point(&self, p: Point2D, tol: f64) -> PointClassification;

    /// Intersect an infinite line with the polygon boundary.
    /// Returns hits sorted by t_line parameter.
    /// Solves 2x2 linear system per edge (Cramer's rule), skips parallel edges.
    /// Deduplicates consecutive hits within 1e-9 (shared vertices).
    pub fn intersect_line(&self, line_origin: Point2D, line_dir: Point2D) -> Vec<LineHit>;
}
```

---

## 7. Utilities (`utils.rs`)

### 7.1 Orthonormal Basis

**KERNEL REQUIRES** (used pervasively — geom, sketch, solvers, builders):

```rust
/// Given a unit axis vector, return two perpendicular unit vectors (u, v)
/// forming a right-handed frame (u, v, axis).
/// Picks arbitrary perpendicular: prefers [1,0,0] if |axis.x| < 0.9, else [0,1,0].
/// u = normalize(axis x arbitrary), v = axis x u.
pub fn orthonormal_basis(axis: &Vec3) -> (Vec3, Vec3);
```

### 7.2 Rodrigues Rotation

**KERNEL REQUIRES** (used in revolve_builder, fillet_builder, arc generation):

```rust
/// Rotate a vector around an axis by angle (radians) using Rodrigues' formula:
/// v_rot = v*cos(theta) + (axis x v)*sin(theta) + axis*(axis . v)*(1 - cos(theta))
pub fn rodrigues_rotate(v: &Vec3, axis: &Vec3, angle: f64) -> Vec3;

/// Rotate a point around an axis through `axis_origin` by angle.
pub fn rodrigues_rotate_point(
    point: &Point3,
    axis_origin: &Point3,
    axis_dir: &Vec3,
    angle: f64,
) -> Point3;
```

### 7.3 Arc Generation

**KERNEL REQUIRES** (used in fillet_builder):

```rust
/// Generate n+1 points along a circular arc via Rodrigues rotation.
/// start: point on the arc at angle=0.
/// center: center of the arc.
/// axis: rotation axis (unit).
/// subtended_angle: total angle in radians.
/// n_segments: number of chord segments.
pub fn generate_arc_points(
    center: &Point3,
    start: &Point3,
    axis: &Vec3,
    subtended_angle: f64,
    n_segments: usize,
) -> Vec<Point3>;
```

### 7.4 AABB

**KERNEL REQUIRES** (used in boolean broad_phase):

```rust
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: Point3,
    pub max: Point3,
}

impl Aabb {
    /// Build from a set of points.
    pub fn from_points(pts: &[Point3]) -> Self;

    /// Test overlap with another AABB (all three axes).
    pub fn overlaps(&self, other: &Aabb) -> bool;

    /// Expand box uniformly by `amount` in all directions.
    pub fn inflate(&mut self, amount: f64);

    /// Merge two AABBs.
    pub fn union(&self, other: &Aabb) -> Aabb;

    /// Test if a point is inside.
    pub fn contains_point(&self, p: &Point3) -> bool;
}

/// Compute AABB inflation amount for a curved surface.
/// Plane: 0. Cylinder/Sphere: |radius|. Cone: tan(half_angle)*conservative_factor.
/// Torus: |minor_radius|.
pub fn surface_aabb_inflation(surface: &SurfaceDef) -> f64;
```

### 7.5 Transform Helpers

**KERNEL REQUIRES** (used in geom apply_transform):

```rust
/// Apply the 3x3 rotation/scale sub-matrix of a Mat4 to a direction vector.
/// Normalizes result; returns original if degenerate.
pub fn transform_dir(m: &Mat4, v: &Vec3) -> Vec3;

/// Extract uniform scale factor from a Mat4 (column 0 norm).
pub fn scale_factor(m: &Mat4) -> f64;
```

---

## 8. Module Dependency Map

```
types          (no internal deps — leaf)
  ↑
utils          (depends on types)
  ↑
curves         (depends on types, utils)
surfaces       (depends on types, utils, curves)  [for iso_curve return type]
polygon2d      (depends on types)
  ↑
mesh           (depends on types, surfaces)
  ↑
intersection   (depends on types, curves, surfaces, mesh, polygon2d, utils)
```

---

## 9. Design Conventions

1. **All angles in radians.** No degree/radian conversion in this crate.

2. **Negative radius = flipped normal.** This convention is used by Cylinder, Sphere,
   Torus (negate radius) and Cone (negate half_angle) to indicate inward-facing normals.
   The geometry crate must preserve this convention since the kernel depends on it.

3. **NURBS via curvo.** Re-export `curvo::NurbsCurve3D` and `curvo::NurbsSurface3D`.
   Do not reimplement NURBS evaluation — curvo is the NURBS engine.

4. **Tolerances as parameters, not constants.** Functions that need a tolerance take it
   as an argument. No global tolerance state.

5. **No topology types.** No half-edges, faces, shells, solids, or IDs. The kernel maps
   between topology indices and geometry objects — that mapping lives in the kernel.

6. **Serde on all public types.** The kernel serializes geometry for `.forge` files.

7. **tracing instrumentation.** Add `debug_span!` on expensive operations (NURBS inverse
   mapping, mesh intersection, SSI refinement). `warn!` on convergence failures.

8. **`#[inline]` on trivial helpers.** `point3_to_f32`, `vec3_to_f32`, single-line
   delegations.

---

## 10. Migration Path from Kernel

The kernel currently implements all this geometry inline. Migration order:

1. **Types + Utils** — move Point3/Vec3/Mat4 re-exports, orthonormal_basis, rodrigues,
   AABB, transform helpers. Kernel re-exports from geometry.
2. **Curve/Surface structs** — move struct definitions and CurveDef/SurfaceDef enums.
   Kernel's AnalyticalGeomStore holds geometry::SurfaceDef instead of its own copy.
3. **Curve/Surface trait impls** — move eval/tangent/normal/domain/inverse_uv.
   Kernel's GeomAccess delegates to geometry::Surface trait.
4. **Polygon2D** — move wholesale from rustkernel-math.
5. **TriMesh + mesh ops** — move FaceMesh → geometry::TriMesh.
6. **Intersection** — move all SSI solvers, ray-surface, tri_tri, mesh intersection.
7. **Solver pipeline** — move SsiPipeline from rustkernel-solvers.

After migration, `rustkernel-math` becomes a thin re-export crate (or is deleted),
`rustkernel-geom` becomes just `AnalyticalGeomStore` + `GeomAccess` impl,
and `rustkernel-solvers` is deleted (absorbed into geometry::intersection).

---

## 11. What the Kernel Keeps

These remain in the kernel because they depend on B-Rep topology:

- `GeomAccess` trait (indexes geometry by `u32` IDs from arena)
- `AnalyticalGeomStore` (stores `Vec<SurfaceDef>`, implements `GeomAccess`)
- `inverse_map_from_kind()` (uses topology's SurfaceKind enum)
- Boolean pipeline (face classifier, curve trimming to faces, face splitter)
- All builders (create topology + assign geometry)
- Euler operators (pure topology)
- Sketch constraint solver (outputs profiles for extrude)
- Diagnostics/validation (topology checks)
- Persistence (.forge, STL, OBJ export)

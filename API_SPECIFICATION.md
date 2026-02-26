# RustKernel API Specification

> **Status:** Design document — not yet implemented.
> This is the north-star public API. Implementation follows this spec, not the other way around.

---

## 1. Philosophy & Principles

**Handle-based kernel context.** `Kernel` owns all topological and geometric data.
Callers receive small, opaque, `Copy` handles. No data is borrowed out; every
query goes through a `Kernel` method.

**Immutable operations.** Booleans, transforms, and feature ops return *new*
handles. The originals remain valid. This enables undo/redo, parametric history,
and safe concurrent reads.

**No global state.** Multiple `Kernel` instances coexist. Each carries its own
arenas, geometry store, and intersection pipeline. Essential for tests and
multi-document apps.

**`[f64; 3]` at the boundary.** The public API never exposes `nalgebra` types
or internal `Idx<T>` indices. All points, vectors, and normals are plain `f64`
arrays. This makes FFI (C, Python/PyO3, WASM/wasm-bindgen) trivial.

**Rust-first, binding-friendly.** The canonical API is Rust. Python and WASM
bindings are thin wrappers generated with `#[pyclass]`/`#[pymethods]` and
`#[wasm_bindgen]`. A C FFI layer exposes handles as plain `u64`.

**Viewer is a consumer.** Visualization is a separate crate that calls
`tessellate` and feeds the result to a renderer. The kernel has no opinion on
display — it produces triangle soup.

---

## 2. Core Types

### Kernel

```rust
pub struct Kernel { /* private: TopoStore, AnalyticalGeomStore, IntersectionPipeline */ }

impl Kernel {
    pub fn new() -> Self;
}
```

`Kernel` is the sole entry point. All methods below are `impl Kernel`.

### Handles

| Handle | Backed by | Description |
|--------|-----------|-------------|
| `SolidHandle` | `u64` | A closed, manifold volume |
| `ShellHandle` | `u64` | A connected set of faces bounding a solid |
| `FaceHandle` | `u64` | A bounded region of a surface |
| `WireHandle` | `u64` | An ordered chain of edges (loop) |
| `EdgeHandle` | `u64` | A bounded curve between two vertices |
| `VertexHandle` | `u64` | A point in 3-space |
| `SketchHandle` | `u64` | A 2D constraint sketch (Phase 5) |
| `Sketch2DHandle` | `u64` | An element within a sketch |

Every handle is `Copy + Eq + Hash + Send + Sync`.

**Encoding:** `[32-bit generation | 32-bit index]`. The low 32 bits are the
arena index (maps to current `Idx<T>.raw: u32`). The high 32 bits are a
generation counter that detects use-after-delete.

### KernelError

```rust
pub enum KernelError {
    InvalidHandle,      // Handle is stale or was never issued
    GeometryError(String),
    BooleanFailed(String),
    TopologyError(String),
    TessellationError(String),
    ConstraintError(String),
    IoError(std::io::Error),
    NotImplemented,
}
```

All fallible methods return `Result<T, KernelError>`.

---

## 3. The Three API Layers

| Layer | Purpose | Example |
|-------|---------|---------|
| **Builder Ops** | Atomic B-Rep construction | `make_vertex`, `make_edge`, `make_face`, `sew_shell`, `make_solid` |
| **Primitives** | Convenience shapes (composed from builder ops) | `make_box`, `make_cylinder`, `make_sphere` |
| **Feature Ops** | High-level modeling operations | `fuse`, `cut`, `fillet`, `extrude` |

Builder ops are the foundation. `make_box` is syntactic sugar that internally
calls `make_vertex` → `make_edge` → `make_wire` → `make_face` → `sew_shell` →
`make_solid`. Power users and custom shape builders use builder ops directly.

---

## 4. Phase 1 — Foundation

> What exists in the codebase today, wrapped in a proper `Kernel` API.

### Builder Operations

```rust
fn make_vertex(&mut self, point: [f64; 3]) -> VertexHandle

fn make_edge(
    &mut self, v1: VertexHandle, v2: VertexHandle,
) -> Result<EdgeHandle, KernelError>

fn make_wire(
    &mut self, edges: &[EdgeHandle],
) -> Result<WireHandle, KernelError>

fn make_face(
    &mut self, wire: WireHandle, surface: SurfaceDef,
) -> Result<FaceHandle, KernelError>

fn sew_shell(
    &mut self, faces: &[FaceHandle],
) -> Result<ShellHandle, KernelError>

fn make_solid(
    &mut self, shell: ShellHandle,
) -> Result<SolidHandle, KernelError>
```

#### SurfaceDef

`SurfaceDef` is the public enum users pass to `make_face`. It mirrors
`SurfaceKind` (the internal query result) but is the *input* type at the API
boundary.

```rust
pub enum SurfaceDef {
    Plane { origin: [f64; 3], normal: [f64; 3] },
    // Phase 3 adds: Cylinder, Sphere, Cone, Torus
}
```

### Primitives

```rust
fn make_box(
    &mut self, dx: f64, dy: f64, dz: f64,
) -> Result<SolidHandle, KernelError>

fn make_box_at(
    &mut self, origin: [f64; 3], dx: f64, dy: f64, dz: f64,
) -> Result<SolidHandle, KernelError>
```

`make_box` decomposes internally to:
1. 8 × `make_vertex` (corner points)
2. 12 × `make_edge` (box edges)
3. 6 × `make_wire` (face loops, 4 edges each)
4. 6 × `make_face` (with `SurfaceDef::Plane`)
5. 1 × `sew_shell` (twin-matches half-edges)
6. 1 × `make_solid`

### Transforms

All transforms return new handles; originals stay valid.

```rust
fn translate(
    &mut self, solid: SolidHandle, delta: [f64; 3],
) -> Result<SolidHandle, KernelError>

fn rotate(
    &mut self, solid: SolidHandle, axis: [f64; 3], angle_rad: f64,
) -> Result<SolidHandle, KernelError>

fn transform(
    &mut self, solid: SolidHandle, matrix: [[f64; 4]; 4],
) -> Result<SolidHandle, KernelError>

fn mirror(
    &mut self, solid: SolidHandle,
    plane_origin: [f64; 3], plane_normal: [f64; 3],
) -> Result<SolidHandle, KernelError>
```

### Topology Queries

```rust
fn solid_faces(&self, solid: SolidHandle) -> Result<Vec<FaceHandle>, KernelError>
fn solid_edges(&self, solid: SolidHandle) -> Result<Vec<EdgeHandle>, KernelError>
fn face_edges(&self, face: FaceHandle) -> Result<Vec<EdgeHandle>, KernelError>
fn edge_vertices(&self, edge: EdgeHandle) -> Result<(VertexHandle, VertexHandle), KernelError>
fn vertex_position(&self, vertex: VertexHandle) -> Result<[f64; 3], KernelError>
fn face_normal(&self, face: FaceHandle, uv: Option<[f64; 2]>) -> Result<[f64; 3], KernelError>
```

### Tessellation

```rust
fn tessellate(
    &mut self, solid: SolidHandle,
) -> Result<TessellationResult, KernelError>

fn tessellate_with_options(
    &mut self, solid: SolidHandle, opts: TessellationOptions,
) -> Result<TessellationResult, KernelError>
```

```rust
pub struct TessellationResult {
    pub positions: Vec<f32>,   // [x,y,z, x,y,z, …]  flat
    pub normals:   Vec<f32>,   // [nx,ny,nz, …]       flat
    pub indices:   Vec<u32>,   // triangle indices
    pub face_ranges: Vec<FaceRange>,  // per-face picking
}

pub struct FaceRange {
    pub face: FaceHandle,
    pub index_offset: u32,
    pub index_count: u32,
}

pub struct TessellationOptions {
    pub angular_tolerance: f64,   // radians, default π/18 (10°)
    pub chordal_tolerance: f64,   // absolute, default 0.01
}
```

Positions and normals are `f32` for GPU consumption. Internal math stays `f64`.

### Copy

```rust
fn copy_solid(
    &mut self, solid: SolidHandle,
) -> Result<SolidHandle, KernelError>
```

Deep-copies all topology and geometry under the solid into fresh arena slots.

---

## 5. Phase 2 — Booleans

```rust
fn fuse(
    &mut self, a: SolidHandle, b: SolidHandle,
) -> Result<SolidHandle, KernelError>

fn cut(
    &mut self, a: SolidHandle, b: SolidHandle,
) -> Result<SolidHandle, KernelError>

fn common(
    &mut self, a: SolidHandle, b: SolidHandle,
) -> Result<SolidHandle, KernelError>

fn fuse_many(
    &mut self, solids: &[SolidHandle],
) -> Result<SolidHandle, KernelError>

fn section(
    &mut self, solid: SolidHandle,
    plane_origin: [f64; 3], plane_normal: [f64; 3],
) -> Result<Vec<WireHandle>, KernelError>
```

Naming follows FreeCAD convention (`fuse`/`cut`/`common` instead of
`union`/`subtract`/`intersect`). Initially plane-only; extends as solvers land.

Internally dispatches through `IntersectionPipeline`, which holds a registry of
`SurfaceSurfaceSolver` implementations (currently: `PlanePlaneSolver`).

---

## 6. Phase 3 — Extended Primitives & Analytical Surfaces

```rust
fn make_cylinder(
    &mut self, radius: f64, height: f64,
) -> Result<SolidHandle, KernelError>

fn make_sphere(
    &mut self, radius: f64,
) -> Result<SolidHandle, KernelError>

fn make_cone(
    &mut self, r1: f64, r2: f64, height: f64,
) -> Result<SolidHandle, KernelError>

fn make_torus(
    &mut self, major_r: f64, minor_r: f64,
) -> Result<SolidHandle, KernelError>
```

Each new primitive extends `SurfaceDef` and `SurfaceKind`:

```rust
pub enum SurfaceDef {
    Plane    { origin: [f64; 3], normal: [f64; 3] },
    Cylinder { origin: [f64; 3], axis: [f64; 3], radius: f64 },
    Sphere   { center: [f64; 3], radius: f64 },
    Cone     { apex: [f64; 3], axis: [f64; 3], half_angle: f64 },
    Torus    { center: [f64; 3], axis: [f64; 3], major_r: f64, minor_r: f64 },
}
```

Each analytical surface pair requires a new `SurfaceSurfaceSolver`
(Plane-Cylinder, Cylinder-Cylinder, etc.). Builder ops use `SurfaceDef`
variants to construct faces on these surfaces.

---

## 7. Phase 4 — Feature Operations

```rust
fn fillet(
    &mut self, solid: SolidHandle, edges: &[EdgeHandle], radius: f64,
) -> Result<SolidHandle, KernelError>

fn chamfer(
    &mut self, solid: SolidHandle, edges: &[EdgeHandle], distance: f64,
) -> Result<SolidHandle, KernelError>

fn extrude(
    &mut self, profile: FaceHandle, direction: [f64; 3], distance: f64,
) -> Result<SolidHandle, KernelError>

fn revolve(
    &mut self, profile: FaceHandle,
    axis_origin: [f64; 3], axis_dir: [f64; 3], angle_rad: f64,
) -> Result<SolidHandle, KernelError>

fn sweep(
    &mut self, profile: FaceHandle, path: WireHandle,
) -> Result<SolidHandle, KernelError>

fn loft(
    &mut self, profiles: &[WireHandle], solid: bool,
) -> Result<SolidHandle, KernelError>

fn shell(
    &mut self, solid: SolidHandle, faces_to_remove: &[FaceHandle], thickness: f64,
) -> Result<SolidHandle, KernelError>

fn draft(
    &mut self, solid: SolidHandle, faces: &[FaceHandle],
    neutral_normal: [f64; 3], angle_rad: f64,
) -> Result<SolidHandle, KernelError>
```

### Selection Helpers

These feed into feature ops — select entities by geometric criteria.

```rust
fn find_sharp_edges(
    &self, solid: SolidHandle, angle_rad: f64,
) -> Result<Vec<EdgeHandle>, KernelError>

fn find_face_by_normal(
    &self, solid: SolidHandle, direction: [f64; 3], tol_rad: f64,
) -> Result<Vec<FaceHandle>, KernelError>

fn find_nearest_face(
    &self, solid: SolidHandle, point: [f64; 3],
) -> Result<FaceHandle, KernelError>
```

---

## 8. Phase 5 — 2D Sketching

```rust
fn begin_sketch(
    &mut self, face: FaceHandle,
) -> Result<SketchHandle, KernelError>

fn begin_sketch_on_plane(
    &mut self, origin: [f64; 3], normal: [f64; 3],
) -> Result<SketchHandle, KernelError>

fn sketch_add_line(
    &mut self, sk: SketchHandle, start: [f64; 2], end: [f64; 2],
) -> Result<Sketch2DHandle, KernelError>

fn sketch_add_arc(
    &mut self, sk: SketchHandle,
    center: [f64; 2], r: f64, start_angle: f64, end_angle: f64,
) -> Result<Sketch2DHandle, KernelError>

fn sketch_add_circle(
    &mut self, sk: SketchHandle, center: [f64; 2], r: f64,
) -> Result<Sketch2DHandle, KernelError>

fn sketch_add_constraint(
    &mut self, sk: SketchHandle, constraint: ConstraintKind,
) -> Result<(), KernelError>

fn sketch_solve(
    &mut self, sk: SketchHandle,
) -> Result<SketchSolveResult, KernelError>

fn sketch_to_wire(
    &mut self, sk: SketchHandle,
) -> Result<WireHandle, KernelError>

fn sketch_to_face(
    &mut self, sk: SketchHandle,
) -> Result<FaceHandle, KernelError>
```

### ConstraintKind

```rust
pub enum ConstraintKind {
    Coincident(Sketch2DHandle, Sketch2DHandle),
    Horizontal(Sketch2DHandle),
    Vertical(Sketch2DHandle),
    Distance(Sketch2DHandle, Sketch2DHandle, f64),
    Angle(Sketch2DHandle, Sketch2DHandle, f64),
    Tangent(Sketch2DHandle, Sketch2DHandle),
    Parallel(Sketch2DHandle, Sketch2DHandle),
    Perpendicular(Sketch2DHandle, Sketch2DHandle),
    FixedPoint(Sketch2DHandle, [f64; 2]),
    Radius(Sketch2DHandle, f64),
    Equal(Sketch2DHandle, Sketch2DHandle),
}
```

### SketchSolveResult

```rust
pub enum SketchSolveResult {
    FullyConstrained,
    UnderConstrained { dof: usize },
    OverConstrained,
    SolveFailed(String),
}
```

---

## 9. Phase 6 — NURBS & STEP I/O

```rust
fn make_nurbs_curve(
    &mut self,
    control_points: &[[f64; 3]], weights: &[f64],
    knots: &[f64], degree: usize,
) -> Result<WireHandle, KernelError>

fn interpolate_curve(
    &mut self, points: &[[f64; 3]],
) -> Result<WireHandle, KernelError>

fn import_step(&mut self, path: &str) -> Result<Vec<SolidHandle>, KernelError>
fn export_step(&self, solids: &[SolidHandle], path: &str) -> Result<(), KernelError>

fn import_step_from_bytes(&mut self, data: &[u8]) -> Result<Vec<SolidHandle>, KernelError>
fn export_step_to_bytes(&self, solids: &[SolidHandle]) -> Result<Vec<u8>, KernelError>

fn export_stl(
    &self, solid: SolidHandle, path: &str,
) -> Result<(), KernelError>
```

`_from_bytes` / `_to_bytes` variants are critical for WASM (no filesystem
access in the browser).

---

## 10. Phase 7 — Scripting Bindings

### Python (PyO3)

`Kernel` → `#[pyclass]`, handles → `#[pyclass]`, methods → `#[pymethods]`.
Mesh data returned as NumPy arrays via buffer protocol.

### WASM (wasm-bindgen)

`Kernel` → JS class, handles → JS objects, mesh data → `Float32Array`/`Uint32Array`.

### C FFI

`extern "C"` functions, handles are plain `u64`, errors via return code +
`rk_last_error() -> *const c_char`.

### Example: Python workflow

```python
import rustkernel as rk

k = rk.Kernel()

# Create geometry
box_  = k.make_box(10.0, 10.0, 10.0)
cyl   = k.make_cylinder(3.0, 15.0)
cyl   = k.translate(cyl, [5.0, 5.0, 0.0])

# Boolean cut
result = k.cut(box_, cyl)

# Fillet the top edges
edges = k.find_sharp_edges(result, 1.2)
result = k.fillet(result, edges, 1.0)

# Tessellate for display
mesh = k.tessellate(result)
print(f"triangles: {len(mesh.indices) // 3}")

# Export
k.export_step([result], "part.step")
```

---

## 11. Measurement & Analysis

```rust
fn edge_length(&self, edge: EdgeHandle) -> Result<f64, KernelError>
fn face_area(&self, face: FaceHandle) -> Result<f64, KernelError>
fn solid_volume(&self, solid: SolidHandle) -> Result<f64, KernelError>
fn solid_surface_area(&self, solid: SolidHandle) -> Result<f64, KernelError>

fn bounding_box(
    &self, solid: SolidHandle,
) -> Result<BoundingBox, KernelError>

fn point_classify(
    &self, solid: SolidHandle, point: [f64; 3],
) -> Result<PointClassification, KernelError>

fn mass_properties(
    &self, solid: SolidHandle, density: f64,
) -> Result<MassProperties, KernelError>

fn distance_solid_solid(
    &self, a: SolidHandle, b: SolidHandle,
) -> Result<f64, KernelError>
```

### Supporting types

```rust
pub struct BoundingBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

pub enum PointClassification {
    Inside,
    Outside,
    OnFace(FaceHandle),
    OnEdge(EdgeHandle),
    OnVertex(VertexHandle),
}

pub struct MassProperties {
    pub mass: f64,
    pub center_of_gravity: [f64; 3],
    pub inertia_tensor: [[f64; 3]; 3],
}
```

---

## Appendix A: Handle Type Catalog

| RustKernel | OpenCASCADE | FreeCAD (Python) | Notes |
|------------|-------------|-----------------|-------|
| `SolidHandle` | `TopoDS_Solid` | `Part.Solid` | Closed manifold volume |
| `ShellHandle` | `TopoDS_Shell` | `Part.Shell` | Connected face set |
| `FaceHandle` | `TopoDS_Face` | `Part.Face` | Bounded surface region |
| `WireHandle` | `TopoDS_Wire` | `Part.Wire` | Ordered edge chain |
| `EdgeHandle` | `TopoDS_Edge` | `Part.Edge` | Bounded curve (maps to internal `Edge` + half-edge pair) |
| `VertexHandle` | `TopoDS_Vertex` | `Part.Vertex` | Point in 3-space |
| `SketchHandle` | — | `Sketcher.Sketch` | 2D constraint sketch |

---

## Appendix B: KernelError Catalog

| Variant | When raised | Recovery |
|---------|-------------|----------|
| `InvalidHandle` | Handle generation mismatch or out-of-range index | Obtain fresh handle via query |
| `GeometryError` | Degenerate input (zero-length edge, zero-area face) | Fix input parameters |
| `BooleanFailed` | Intersection pipeline cannot solve or classify | Simplify geometry, check overlaps |
| `TopologyError` | Non-manifold result, open shell where solid expected | Check input topology |
| `TessellationError` | Tessellator produces degenerate output | Adjust tolerances |
| `ConstraintError` | Over-constrained sketch, solver divergence | Remove conflicting constraints |
| `IoError` | File read/write failure (STEP, STL export) | Check path and permissions |
| `NotImplemented` | Method exists in API but solver/feature not yet coded | Wait for implementation |

---

## Appendix C: Internal Mapping

> For kernel developers. How public types relate to internal types.

| Public type | Internal type | Location |
|-------------|--------------|----------|
| `SolidHandle` (u64) | `SolidIdx = Idx<Solid>` (u32) | `rustkernel-topology/src/arena.rs` |
| `FaceHandle` (u64) | `FaceIdx = Idx<Face>` (u32) | same |
| `EdgeHandle` (u64) | `EdgeIdx = Idx<Edge>` (u32) | same |
| `VertexHandle` (u64) | `VertexIdx = Idx<Vertex>` (u32) | same |
| `SurfaceDef` | `Plane` struct / `SurfaceKind` | `rustkernel-primitives/src/geom.rs`, `rustkernel-topology/src/geom_store.rs` |
| `TessellationResult` | `FaceMesh` per face | `rustkernel-topology/src/mesh_cache.rs` |
| `Kernel` | `TopoStore` + `AnalyticalGeomStore` + `IntersectionPipeline` | `rustkernel-topology/src/store.rs`, `rustkernel-primitives/src/geom.rs`, `rustkernel-topology/src/intersection.rs` |

Handle encoding: `(generation: u32) << 32 | (Idx<T>.raw(): u32)`. The
`Kernel` maintains a parallel generation array per arena. On deletion the
generation increments; stale handles fail with `InvalidHandle`.

`[f64; 3]` ↔ `Point3`/`Vec3` conversion happens at the `Kernel` method
boundary. Internal code continues to use `nalgebra` types freely.

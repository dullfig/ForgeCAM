# Kernel TODO

Last updated: 2025-03-10
Tests: 190 passing, 1 known failure (test_volume_sphere)

---

## Legend

- `[x]` Done
- `[ ]` Not started
- `[~]` Partial / has caveats

---

## Primitives & Builders

- [x] Box (make_box, make_box_at)
- [x] Cylinder (16-face lateral, top/bottom caps)
- [x] Sphere (lat/lon tessellation)
- [x] Cone / frustum
- [x] Torus
- [x] Extrude (planar profile → solid)
- [x] Revolve (profile around axis)
- [x] NURBS extrude (NURBS profile → solid)
- [x] NURBS revolve (NURBS profile → solid)
- [x] Loft surface (NURBS loft between profiles)
- [ ] **Solid loft** — close two profile loops with lofted skin + caps. NURBS loft surface exists; need topology builder.
- [ ] **Sweep along path** — extrude profile along 3D curve. Frenet frame or fixed orientation. Produces ruled/NURBS side surfaces.
- [ ] **Pattern (linear)** — copy solid N times with translation. Trivial: transform + fuse.
- [ ] **Pattern (circular)** — copy solid N times around axis. Trivial: rotate + fuse.

## Booleans

- [x] Fuse (union)
- [x] Cut (subtract)
- [x] Common (intersect)
- [x] fuse_many
- [x] Section (cross-section plane → loops)
- [x] Plane-plane intersections (box-box)
- [x] Plane-cylinder intersections (box-cylinder)
- [x] Plane-sphere intersections (box-sphere)
- [x] Interior loop splitting (slit-based, for closed curves inside a face)
- [x] Segment chaining (chain chord segments into polylines)
- [x] SSI result caching + surface-pair dedup
- [x] Genus > 0 Euler validation (V-E+F = 2-2g)
- [ ] **Cylinder-cylinder** — needs solver (plane_cylinder handles plane-cyl only)
- [ ] **Sphere-cylinder, cone-cylinder, etc.** — needs general SSI solver or NURBS fallback
- [ ] **Torus intersections** — plane-torus solver exists; torus-anything else missing

## Local Operations

- [x] Chamfer (planar faces, straight edges)
- [x] Fillet (planar faces, straight edges, equal radii)
- [x] Euler chamfer (topology-preserving)
- [x] Euler fillet (topology-preserving)
- [x] Dissimilar-radius fillets (per-edge radii, 4-tier dispatch)
- [ ] **Fillet on curved edges** — current fillet requires straight edges between planar faces. Need rolling-ball center path on cylinder-cylinder, sphere-cylinder, arbitrary surface pairs. Hard.
- [ ] **Surface offset** — parallel surface at distance d. Plane→Plane, Cylinder→Cyl(r±d), Sphere→Sph(r±d), Cone→Cone, NURBS→NURBS (control point displacement). Foundation for shell.
- [ ] **Shell / hollow** — offset all faces inward, subtract from original. Depends on surface offset + robust booleans.
- [ ] **Draft / taper** — tilt selected faces by angle relative to pull direction. Modifies face surface + recomputes edges. Critical for injection molding.

## Transforms

- [x] Translate
- [x] Rotate
- [x] Scale
- [x] Mirror
- [x] General 4x4 affine transform

## Analysis & Queries

- [x] Volume (divergence theorem on tessellated shell)
- [x] Surface area
- [x] Center of mass
- [x] Inertia tensor
- [x] Bounding box
- [~] **Sphere volume 25% off** — test_volume_sphere fails. Likely coarse tessellation for divergence theorem. Needs finer mesh or analytical path for analytical surfaces.
- [ ] **Curvature analysis** — principal curvatures at surface point, Gaussian/mean curvature fields
- [ ] **Interference / clearance check** — two solids: intersecting? min distance?
- [ ] **Mesh quality metrics** — triangle aspect ratio, min angle, self-intersection
- [ ] **Watertight validation** — beyond Euler: manifold check, consistent normals, no T-junctions

## Spatial Queries & Interaction

- [~] Ray-cast point-in-solid (exists in face_classifier, not exposed as API)
- [ ] **Closest point on face/edge/vertex** — spatial query
- [ ] **BVH acceleration** — persistent bounding volume hierarchy for shells
- [ ] **Face/edge filtering** — "all faces with normal within 10deg of Z", predicate-based selection
- [ ] **Topology walking API** — given face→adjacent faces, edge→two faces, vertex→all edges. Exists internally, not exposed.
- [ ] **Pick / hit-test** — screen coords + camera → face under cursor

## Export & Persistence

- [x] STL export (ASCII + binary)
- [x] OBJ export
- [x] .forge native save/load (bincode, round-trip tested)
- [ ] **STEP AP214 export** — industry standard. Enormous spec. Likely subset first.
- [ ] **STEP import** — parse STEP → build geometry → B-Rep assembly → repair. Interactive repair needs GUI.
- [ ] **IGES import/export** — legacy format, still used in shops
- [ ] **3MF export** — modern mesh format, better than STL for color/materials

## Sketcher

- [x] 2D constraint solver (Newton-Raphson + SVD)
- [x] Distance, angle, coincident, horizontal, vertical, tangent constraints
- [x] Profile extraction (sketch → wire → extrude)
- [ ] **Sketch on face** — project sketch plane onto existing face. Needs face→workplane.
- [ ] **Fillet/chamfer in sketch** — 2D fillet/chamfer on sketch profiles
- [ ] **Offset curve in sketch** — parallel curve at distance

## Robustness & Tech Debt

- [ ] **Kill catch_unwind in topology_builder** — boolean failures panic instead of returning Result. #1 reliability issue. Return Result from build_result_solid, propagate through ops.rs.
- [ ] **Global tolerance struct** — `KernelTolerance { linear: f64, angular: f64 }`. Replace 6+ hardcoded epsilons (1e-6, 1e-8, 1e-10). Thread through solvers, booleans, tessellation, snap.
- [ ] **Arena garbage collection** — Euler kill ops leave orphaned elements. Mark-sweep or compact after operations.
- [ ] **Entity naming / persistent IDs** — `HashMap<EntityId, String>`. Prerequisite for parametric features, UI, scripting. Hardest unsolved architectural problem.
- [ ] **Vertex dedup uses exact f64 bit pattern** — topology_builder's point_hash() has no tolerance-based merge. Works but fragile.

## Performance

- [ ] **BVH caching** — persistent BVH for shells (current AABB broad phase rebuilds every boolean)
- [ ] **Parallel tessellation** — rayon over faces (each face independent)
- [ ] **Parallel boolean broad phase** — face-pair generation is embarrassingly parallel
- [ ] **WASM compilation** — verify all crates compile to wasm32. curvo + nalgebra should work.

## Viewer (dev tool, not production)

- [x] 16 scenes: box, cylinder, sphere, cone, torus, extrude, revolve, chamfer, fillet, booleans, section, NURBS, dissimilar fillet
- [ ] **Boolean result scenes** — show box-cylinder cut, box-sphere cut in viewer
- [ ] **Interactive camera** — orbit/pan/zoom (currently fixed views?)

---

## Suggested Attack Order

```
QUICK WINS (hours each):
  Pattern linear/circular    — transforms + fuse already work
  Solid loft                 — NURBS loft exists, need builder wrapper
  Fix sphere volume          — finer tessellation or analytical path
  Boolean viewer scenes      — show off new curved-surface booleans

MEDIUM (days each):
  Kill catch_unwind          — Result propagation through boolean pipeline
  Surface offset (analytical)— plane/cyl/sphere/cone are trivial formulas
  Shell / hollow             — offset + subtract, huge modeling value
  Global tolerance struct    — pay down debt, user-controllable precision

HARD (weeks each):
  Fillet on curved edges     — rolling-ball geometry on arbitrary surfaces
  Draft / taper              — face surface modification + edge recompute
  Sweep along path           — Frenet frame, NURBS side surfaces
  STEP export subset         — massive spec, start with AP214 geometry

LONG TERM:
  STEP import + repair
  Entity naming / persistent IDs
  Spatial queries + BVH
  Scripting layer (PyO3 / WASM)
```

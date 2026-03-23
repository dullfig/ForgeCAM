# RustKernel Hit List

Status: 251 tests, 9 crates, solid prototype. Time to make it real.

---

## Decision Block: File Format & Persistence

Everything downstream depends on this. Can't demo, can't test workflows, can't
integrate with anything until we can save and load.

**Decision needed:** What format(s) to support?

| Option | Pros | Cons |
|--------|------|------|
| **Custom JSON/bincode** | Fast to implement, full fidelity to our B-Rep | Island — nobody else reads it |
| **STEP AP214** | Industry standard, every CAD tool reads it | Enormous spec, 6-12 months for full compliance |
| **STEP via opencascade-sys** | Leverage OCC's STEP reader/writer from Rust FFI | Defeats "pure Rust" goal, massive C++ dep |
| **glTF 2.0** | Web-native, viewer-friendly, Rust crates exist | Mesh-only, loses exact geometry |
| **STL/OBJ export only** | Trivial (days), immediate value | One-way — can't round-trip |

Likely answer: **custom binary for native save/load** + **STL/OBJ export** now,
**STEP subset** later. But this shapes everything — discuss before coding.

---

## Phase 7: Immediate Wins (days each)

These are low-hanging fruit that make the kernel usable *today*.

- [ ] **7a. STL export** — Walk shell faces, emit triangles from `mesh_cache`.
      Literally just a header + triangle normals + vertices. ASCII and binary.
      Depends on: nothing (tessellation already works).

- [ ] **7b. OBJ export** — Same source data, different format. Vertex sharing
      across faces for smaller files.

- [ ] **7c. General affine transforms** — `rotate`, `scale`, `mirror` on solids.
      Extend existing `translate`. Surface/curve types already have `translate`;
      need to generalize to `Mat4`. curvo's `Transformable` already handles NURBS.
      Unblocks: pattern, mirror-body, assembly.

- [ ] **7d. Mass properties** — Volume via divergence theorem on tessellated shell
      (sum signed tetrahedra volumes). Surface area = sum triangle areas.
      Center of mass, moments of inertia follow from same loop.
      Depends on: tessellation (done).

- [ ] **7e. Native save/load** — Serialize `TopoStore` + `AnalyticalGeomStore` to
      binary (bincode or postcard). Round-trip test: build box → save → load →
      `validate_solid` passes. Design the format version header now.

---

## Phase 8: Robustness & Technical Debt

These aren't features but they're blocking trust in the kernel.

- [ ] **8a. Kill `catch_unwind` in booleans** — `topology_builder::build_result_solid`
      panics on twin-matching failure. Needs to return `Result`, propagate errors
      through `ops.rs` → kernel API. This is the #1 reliability issue.

- [ ] **8b. Global tolerance struct** — `KernelTolerance { linear: f64, angular: f64 }`.
      Thread it through solvers, boolean pipeline, tessellation, snap.
      Replace the 6+ hardcoded epsilons. User can set "coarse" vs "fine".

- [ ] **8c. Arena garbage collection** — Euler kill ops leave orphaned elements.
      After a sequence of operations, dead HEs/edges/faces accumulate.
      Compact or mark-sweep. Not urgent until models get large.

- [ ] **8d. Entity naming / tagging** — `HashMap<EntityId, String>` on the kernel.
      "top face", "fillet edge #3". Needed for any UI or scripting layer.
      Prerequisite for parametric features.

---

## Phase 9: Core Missing Operations

These are the "any self-respecting kernel" items, in dependency order.

- [ ] **9a. Surface offset** — Given a surface + distance, produce the parallel
      surface. Plane→Plane, Cylinder→Cylinder(r±d), Sphere→Sphere(r±d),
      NURBS→NURBS (control point displacement). Foundation for shell/hollow.

- [ ] **9b. Shell / hollow** — Offset all faces inward, subtract from original.
      One of the most-used CAD operations (make a box into a box with walls).
      Depends on: 9a (surface offset) + robust booleans (8a).

- [ ] **9c. Draft / taper** — Tilt selected faces by angle relative to pull direction.
      Critical for injection molding. Modifies face surface + recomputes edges.

- [ ] **9d. Fillet on curved edges** — Current fillet requires straight edges between
      planar faces. Need: cylinder-cylinder fillet, sphere-cylinder, arbitrary.
      The Euler topology surgery works regardless — it's the geometry computation
      (rolling-ball center path) that needs generalization.

- [ ] **9e. Sweep along path** — Extrude a profile along a 3D curve. Frenet frame
      or fixed orientation. Produces ruled/NURBS side surfaces.

- [ ] **9f. Loft between profiles** — Already have NURBS loft for surfaces; need
      the solid builder: close two profile loops with a lofted skin + caps.

- [ ] **9g. Pattern (linear / circular)** — Copy a feature N times with transform.
      Depends on: 7c (transforms) + booleans (to fuse copies).

---

## Phase 10: Query & Interaction Layer

Needed before any UI or scripting can work.

- [ ] **10a. Spatial queries** — Closest edge/face/vertex to a point. Ray-cast into
      solid (already partially exists in face_classifier). BVH acceleration.

- [ ] **10b. Face/edge filtering** — "All faces with normal within 10° of Z."
      "All edges longer than 5mm." Predicate-based selection.

- [ ] **10c. Topology walking API** — Given a face, get adjacent faces. Given an
      edge, get the two faces. Given a vertex, get all edges. These exist
      internally but aren't exposed as clean kernel API.

- [ ] **10d. Pick / hit-test** — Given screen coordinates + camera, find the face
      under the cursor. Depends on: ray-cast (10a) + viewer integration.

---

## Phase 11: Analysis & Validation

- [ ] **11a. Mass properties** (covered in 7d)

- [ ] **11b. Curvature analysis** — Principal curvatures at a surface point.
      Gaussian and mean curvature fields. For surface quality inspection.

- [ ] **11c. Interference / clearance check** — Given two solids, are they
      intersecting? What's the minimum distance? Uses existing AABB + mesh.

- [ ] **11d. Mesh quality metrics** — Triangle aspect ratio, minimum angle,
      self-intersection detection on tessellated output.

- [ ] **11e. Watertight validation** — Beyond Euler check: verify mesh is
      manifold, all normals consistent, no T-junctions.

---

## Phase 12: Performance

Not urgent yet but will matter at scale.

- [ ] **12a. BVH caching** — Persistent bounding volume hierarchy for shells.
      Current AABB broad phase rebuilds every boolean.

- [ ] **12b. Parallel tessellation** — Rayon over faces. Each face tessellates
      independently. Easy win.

- [ ] **12c. Parallel boolean broad phase** — Face-pair candidate generation
      is embarrassingly parallel.

- [ ] **12d. WASM compilation** — Already a stated goal. Verify all crates
      compile to wasm32-unknown-unknown. curvo + nalgebra should be fine.
      three-d viewer won't (needs WebGL backend).

---

## Open Questions

1. **Scripting layer?** — Python bindings (PyO3)? Lua? WASM API? Or pure Rust CLI?
   Affects how we expose the kernel API.

2. **Viewer direction** — Keep three-d? Move to egui + wgpu? Web viewer via WASM?
   Affects export format priorities.

3. **Target user** — CAD developers integrating the kernel? End users modeling parts?
   Shapes every API decision.

4. **STEP timeline** — Is STEP import/export a hard requirement, or can we ship
   with STL/OBJ + native format and add STEP later?

---

## Suggested Attack Order

```
NOW   7a STL export          (1-2 days, immediate demo value)
      7b OBJ export          (half day, same data)
      7c Transforms          (2-3 days, unblocks everything)
      7d Mass properties     (1-2 days, easy win)
 |
 v
NEXT  8a Kill catch_unwind   (3-5 days, hard but critical)
      8b Tolerance struct    (2 days, pay down debt)
      7e Native save/load    (2-3 days, needs format decision)
      8d Entity naming       (1 day)
 |
 v
THEN  9a Surface offset      (1 week, foundational)
      9b Shell/hollow        (1 week, depends on 9a + 8a)
      9d Curved-edge fillet  (1-2 weeks, geometry is hard)
      9e Sweep along path    (1 week)
 |
 v
LATER 10a-d Query layer      (ongoing)
      9c Draft/taper
      9f Solid loft
      9g Pattern
      STEP subset
```

The file format decision (custom vs STEP vs both) should happen before 7e.
The viewer direction question can wait until after Phase 8.

# rustkernel-boolean — Sub-Agent

You own the boolean pipeline: fuse (union), cut (subtract), common (intersect), and section.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `ops.rs` | ~701 | Public API: `fuse`, `cut`, `common`, `fuse_many` + orchestration |
| `curve_trimming.rs` | ~615 | Trim intersection curves to face boundaries, circle/ellipse/NURBS |
| `face_splitter.rs` | ~426 | Split faces along trimmed intersection curves |
| `face_classifier.rs` | ~412 | Classify faces as inside/outside/on using ray-casting |
| `topology_builder.rs` | ~374 | Reconstruct solid from selected faces (deep copy + twin match) |
| `section.rs` | ~338 | Cross-section: intersect solid with plane, return closed loops |
| `broad_phase.rs` | ~157 | AABB overlap test for face-pair candidates |
| `face_selector.rs` | ~132 | Select faces based on classification + boolean operation type |
| `lib.rs` | ~8 | Module declarations |

## Pipeline

```
broad_phase → SSI solve + curve_trimming → face_splitter → face_classifier → face_selector → topology_builder
```

1. **broad_phase**: AABB overlap test filters face pairs
2. **SSI solve**: intersection solvers (from `rustkernel-solvers`) compute raw intersection curves
3. **curve_trimming**: trim curves to face boundaries (handles line, circle, ellipse, NURBS)
4. **face_splitter**: split faces along trimmed curves, creating new sub-faces
5. **face_classifier**: ray-cast point-in-solid test to classify each face as IN/OUT/ON
6. **face_selector**: select faces based on boolean op type (union keeps OUT+ON, subtract keeps OUT_A+IN_B_flipped, etc.)
7. **topology_builder**: deep-copy selected faces into new solid, merge coincident vertices, match twins, validate

## Key implementation details

- Face classification uses mesh centroids for ray-cast origin
- Coplanar faces detected + overlap-tested separately
- `topology_builder` returns `Result<BuildResult, BuildError>` — no panics, no catch_unwind
- Vertex dedup uses `point_hash()` (exact f64 bit pattern) — no tolerance-based merge
- `snap_to_boundary()` uses 1e-6 tolerance in face splitting

## Conventions

- `tracing::debug_span!` on solver calls, `warn!` on classification failures
- All trimming respects face boundary orientation (CCW outer loop)

## Dependencies

- `rustkernel_topology` — TopoStore, topo types, mesh_cache, GeomAccess, IntersectionPipeline
- `rustkernel_geom` — AnalyticalGeomStore, SurfaceDef, CurveDef
- `rustkernel_math` — Point3, Vec3, Polygon2D, tri_tri
- `rustkernel_solvers` — intersection solvers (via pipeline)

## Future work (from HITLIST.md)

- ~~**8a. Kill catch_unwind**~~ — DONE. `topology_builder` returns `Result`, all tests use proper error handling.
- **8b. Global tolerance struct** — Thread `KernelTolerance` through the pipeline instead of hardcoded epsilons.
- **12b/12c. Parallel broad phase + tessellation** — Rayon over face pairs and faces.

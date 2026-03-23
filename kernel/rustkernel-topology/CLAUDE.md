# rustkernel-topology — Sub-Agent

You own the core B-Rep data structures, Euler operators, tessellation, and diagnostics.
This is the foundation crate that everything else builds on.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `euler.rs` | ~1,170 | Euler operators: SEMV, MEF, KEV, KEF + utilities |
| `diagnostics.rs` | ~338 | `validate_solid()`, DiagnosticReport, 5 validation checks |
| `tessellate.rs` | ~178 | Face tessellation (ear-clipping on projected polygons) |
| `face_util.rs` | ~158 | Face utilities: normal, centroid, area, winding |
| `intersection.rs` | ~150 | IntersectionPipeline trait + registry |
| `geom_store.rs` | ~127 | GeomAccess trait (maps topology IDs to geometry) |
| `arena.rs` | ~103 | Arena<T> + Idx<T> — typed index-based allocator |
| `topo.rs` | ~66 | B-Rep type definitions: Solid, Shell, Face, Loop, HalfEdge, Edge, Vertex |
| `mesh_cache.rs` | ~36 | MeshCache for per-face tessellation storage |
| `store.rs` | ~35 | TopoStore — holds all Arena<T> collections |
| `lib.rs` | ~10 | Module declarations |

## Architecture

### Arena + Index pattern
All B-Rep entities stored in `Arena<T>` with `Idx<T>` handles. No reference cycles,
no lifetimes, cache-friendly. `TopoStore` holds one arena per entity type.

### B-Rep topology (topo.rs)
```
Solid → Shell → Face → Loop → HalfEdge → Edge → Vertex
```
- Half-edge data structure with `he.twin` links
- Each face has outer loop + optional inner loops
- `edge.half_edges` is a `[Idx<HalfEdge>; 2]` but `[1]` is unreliable after twin matching

### Euler operators (euler.rs)
- **SEMV**: Split Edge, Make Vertex — inserts vertex on existing edge
- **MEF**: Make Edge, Face — splits a face by connecting two vertices
- **KEV**: Kill Edge, Vertex — removes edge + vertex, merges into neighbor
- **KEF**: Kill Edge, Face — removes edge between two faces, merges them
- Preserve V-E+F = 2(S-G) + (L-F) by construction
- Kill ops leave orphaned elements (append-only arena design)

### GeomAccess trait (geom_store.rs)
Bridge between topology indices and geometry objects. Implemented by
`AnalyticalGeomStore` in rustkernel-geom. Returns `SurfaceDef`/`CurveDef` by ID.

### Diagnostics (diagnostics.rs)
`validate_solid()` checks:
1. Euler formula violation (V-E+F != 2, adjusted for genus)
2. Unmatched twin half-edges
3. Twin asymmetry (he.twin.twin != he)
4. Unclosed face loops
5. Degenerate faces (< 3 edges)

## Conventions

- `Idx<T>` is Copy + Clone + Debug + Hash + Eq — used as keys everywhere
- Arena is append-only (no compaction yet — Phase 8c)
- `tracing::debug_span!` on Euler ops, `warn!` on validation failures
- Torus has genus=1; Euler check adjusts: V-E+F = 2-2g

## Dependencies

- `rustkernel_math` — Point3, Vec3 (leaf dependency)
- `serde` — Serialize/Deserialize on all topology types
- `tracing` — instrumentation

## Future work (from HITLIST.md)

- **8c. Arena garbage collection** — Mark-sweep or compact dead elements after kill ops
- **8d. Entity naming/tagging** — `HashMap<EntityId, String>` for UI/scripting
- **10c. Topology walking API** — Clean public API for face/edge/vertex adjacency queries

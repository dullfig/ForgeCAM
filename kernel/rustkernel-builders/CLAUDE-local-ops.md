# rustkernel-builders — Local Ops Sub-Agent

You own chamfer, fillet, and Euler-based local topology surgery. This is the most complex
and actively evolving part of the builders crate (~4,000 lines).

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `chamfer_builder.rs` | ~540 | Rebuild-based chamfer (single edges, no shared vertices) |
| `fillet_builder.rs` | ~602 | Rebuild-based fillet (single edges, no shared vertices) |
| `euler_chamfer.rs` | ~1,261 | Euler-based chamfer: multi-edge, shared corners, modular 9-op sequence |
| `euler_fillet.rs` | ~1,665 | Euler-based fillet: multi-edge, 3-tier corner blends |

You also read from (but don't own) `edge_analysis.rs` — that belongs to the sweeps agent.

## NOT your files

Primitive builders and sweep builders belong to other agents.
Do not modify `lib.rs` unless adding a new local-ops module.

## Architecture

### Rebuild vs Euler approach
- **Rebuild** (chamfer_builder, fillet_builder): reconstructs the entire solid from scratch.
  Works for isolated edges. Breaks on shared vertices (Euler violation).
- **Euler** (euler_chamfer, euler_fillet): local topology surgery using SEMV/MEF/KEV/KEF.
  Only touches affected edges/faces. Handles multi-edge operations and shared corners.
  This is the production approach — the rebuild versions exist for simple cases.

### Euler chamfer/fillet sequence (9 ops)
1. Compute contact points on each edge
2. SEMV: split edges at contact points
3. MEF: split adjacent faces along rail curves
4. Insert chamfer/blend face with appropriate surface
5. KEV: remove original sharp edge segment
6. Corner resolution (shared vertices)
7. Cross-edge insertion at corners
8. Corner face creation
9. Mesh cache update

### 3-tier corner blend system (euler_fillet)
- **Tier 1**: Per-edge blend surfaces (cylinder for fillet, plane for chamfer)
  with per-edge distance parameters
- **Tier 2**: Sphere patch for orthogonal corners with equal radii
  (analytical SphereSurface — exact, fast)
- **Tier 3**: NURBS loft for general angle/mixed-radius corners
  (skinned surface through boundary curves)

### Corner classification (from SMLib NxM framework)
- **3x1**: 1 of 3 edges blended — trim end face at unfilleted edges
- **3x2**: 2 of 3 — join if same radius+tangent, otherwise bevel
- **3x3 univex**: all same convexity — sphere (orthogonal) or NURBS patch
- **NxN mixed**: mixed convexity — additional spring edges + (N+k)-sided patch

## Key implementation details

- Euler operators live in `rustkernel_topology::euler` (SEMV, MEF, KEV, KEF)
- `edge.half_edges[1]` is UNRELIABLE after twin matching — always use `he.twin`
- Kill ops leave orphaned elements (append-only arena) — this is by design
- Contact point placement: `t * edge_length` from each vertex, where `t = distance / edge_length`
- Uses `edge_analysis::adjacent_faces()` and `edge_analysis::edge_convexity()` heavily

## Conventions

- `tracing::info_span!` on entry, `warn!` on topology failures
- All angles in radians
- Validate Euler invariant after each operation during development

## Dependencies

- `rustkernel_topology` — TopoStore, euler ops, arena, topo types, mesh_cache, GeomAccess
- `rustkernel_geom` — AnalyticalGeomStore, SurfaceDef (Plane, Cylinder, Sphere), CurveDef
- `rustkernel_math` — Point3, Vec3
- `curvo` — NurbsSurface3D for Tier 3 corner blends

## Future work (from HITLIST.md)

- **9d. Fillet on curved edges** — Current fillet requires straight edges between planar
  faces. Needs rolling-ball center path on curved edges (geometry is the hard part;
  the Euler surgery works regardless). This is your domain.
- **9a. Surface offset** — Foundation for shell/hollow. May interact with fillet.
- **Phase 6 Solid Offset** — See `offset_solid.md` in memory. Related but distinct.

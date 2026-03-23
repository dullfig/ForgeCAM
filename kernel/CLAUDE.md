# forgecam-kernel

B-Rep solid modeling kernel for the ForgeCAM CAM system.

## What this crate is

The topological core — half-edge B-Rep data structures, boolean operations, solid builders,
Euler operators, chamfer/fillet, persistence, and export. Analogous to SMLib's TSNLib+SMLib
layers: everything at and above trimmed-surface topology.

## What this crate is NOT

Not a geometry library. Pure geometric computation (curve/surface evaluation, intersection
algorithms, meshing) will migrate to `forgecam-geometry`. This crate maps between topology
indices and geometry objects via the `GeomAccess` trait.

## Key files

- `HITLIST.md` — Prioritized backlog (Phases 7-12)
- `API_SPECIFICATION.md` — North-star public API design
- 9 sub-crates (see workspace structure below)

## Workspace Structure

```
rustkernel-math        — nalgebra types, Polygon2D, tri_tri intersection
rustkernel-topology    — arena, TopoStore, B-Rep types, tessellation, GeomAccess trait
rustkernel-geom        — SurfaceDef/CurveDef enums, AnalyticalGeomStore
rustkernel-solvers     — 7 intersection solvers (5 analytical + 2 NURBS)
rustkernel-builders    — primitive builders, extrude/revolve, chamfer/fillet, euler_*
rustkernel-boolean     — boolean pipeline (8 modules)
rustkernel-sketch      — 2D parametric sketcher, constraint solver (NR+SVD)
rustkernel-kernel      — thin Kernel orchestrator
rustkernel-viewer      — three-d viewer (dev tool, not production GUI)
```

Dep graph: math → topology → geom → solvers/builders/boolean/sketch → kernel → viewer

## Sub-Agent Delegation

This CLAUDE.md is the **kernel coordinator**. When work touches a specific sub-crate,
delegate to the appropriate sub-agent. Each has its own CLAUDE.md with file ownership,
patterns, conventions, and future work.

| Domain | CLAUDE.md | Owns |
|--------|-----------|------|
| Primitives | `rustkernel-builders/CLAUDE-primitives.md` | box, cylinder, sphere, cone, torus builders |
| Sweeps | `rustkernel-builders/CLAUDE-sweeps.md` | extrude, revolve, nurbs_extrude, nurbs_revolve, edge_analysis |
| Local ops | `rustkernel-builders/CLAUDE-local-ops.md` | chamfer, fillet, euler_chamfer, euler_fillet |
| Boolean | `rustkernel-boolean/CLAUDE.md` | boolean pipeline (8 modules) |
| Topology | `rustkernel-topology/CLAUDE.md` | arena, store, topo types, Euler ops, tessellation, diagnostics |
| Solvers | `rustkernel-solvers/CLAUDE.md` | intersection solvers (5 analytical + 2 NURBS) |
| Sketch | `rustkernel-sketch/CLAUDE.md` | 2D constraint solver, workplane, profile extraction |
| Orchestrator | `rustkernel-kernel/CLAUDE.md` | Kernel API, persistence, transforms, mass properties, export |

### Delegation rules

1. **Read the sub-agent's CLAUDE.md** before working in its files
2. **Single-crate work**: delegate directly to the sub-agent
3. **Cross-crate work** (new shared type, API change): design the interface here first,
   then update both sub-agent CLAUDE.md files
4. **builders is split 3 ways**: primitives, sweeps, and local-ops each have separate
   CLAUDE.md files within `rustkernel-builders/`. Route to the right one based on which
   files are affected
5. **math, geom, viewer** are small enough (~130-830 lines) to not need sub-agents.
   Handle directly from this coordinator level

## Status

- 278 tests passing
- Phases 1-7 complete (foundation through persistence/export/transforms/mass properties)
- See HITLIST.md for next priorities (Phase 8: robustness, Phase 9: core missing operations)

## Dependencies

- `nalgebra = "0.34"` — must match forgecam-geometry
- `curvo = "0.1.81"` — must match forgecam-geometry
- `serde = "1"` — serialization for .forge format
- `tracing = "0.1"` — in 6 of 9 crates

## Conventions

- Arena + index pattern for B-Rep topology (Idx<T> handles, no reference cycles)
- Negative radius = flipped normal (Cylinder, Sphere, Torus; Cone uses negative half_angle)
- All angles in radians
- `tracing` instrumentation: info_span! on builders, debug_span! on solvers, warn! on errors
- `build_face_from_vert_idxs()` and `match_twins_from_map()` shared builder helpers
- Euler operators preserve V-E+F=2 by construction; kill ops leave orphans (append-only arena)
- `edge.half_edges[1]` is unreliable after twin matching — always use `he.twin`

## Migration plan (future)

When forgecam-geometry is implemented, pure geometry will migrate out:
- `rustkernel-math` → re-exports from geometry (or deleted)
- `rustkernel-geom` → just AnalyticalGeomStore + GeomAccess impl (structs move to geometry)
- `rustkernel-solvers` → absorbed into geometry::intersection
- Kernel keeps: topology, boolean pipeline, builders, Euler ops, sketch, persistence

# rustkernel-builders — Primitives Sub-Agent

You own the primitive solid builders. These create closed B-Rep solids from scratch.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `box_builder.rs` | ~243 | Axis-aligned box from center + dimensions |
| `cylinder_builder.rs` | ~354 | Cylinder from center + radius + height + segments |
| `sphere_builder.rs` | ~309 | UV-sphere from center + radius + slices/stacks |
| `cone_builder.rs` | ~238 | Cone/frustum from center + r1 + r2 + height |
| `torus_builder.rs` | ~226 | Torus from center + major/minor radii |

## NOT your files

Everything else in `rustkernel-builders/src/` belongs to the sweeps or local-ops agents.
Do not modify `lib.rs`, `edge_analysis.rs`, or any euler/chamfer/fillet/extrude/revolve files.

## Patterns

All primitive builders follow the same structure:

1. **Signature**: `make_X_into(topo: &mut TopoStore, geom: &mut AnalyticalGeomStore, ...) -> SolidIdx`
2. **Geometry first**: push `SurfaceDef` and `CurveDef` entries into `geom`
3. **Topology**: create vertices, edges, half-edges, faces, loops, shell, solid in `topo`
4. **Helpers**: `build_face_from_vert_idxs()` builds a face from a polygon of vertex indices;
   `match_twins_from_map()` matches half-edge twins via vertex-pair map
5. **Mesh cache**: call `topo.mesh_cache_mut().tessellate_solid(...)` at the end
6. **Tests**: every builder has inline `#[cfg(test)]` tests validating Euler (V-E+F=2),
   vertex count, face count, and `validate_solid()` passes

## Conventions

- Negative radius = flipped surface normal (Cylinder, Sphere, Torus). Cone uses negative half_angle.
- All angles in radians
- `tracing::info_span!` on each builder function
- Polygon-approximated edges with exact analytical surfaces

## Dependencies

These files import from:
- `rustkernel_topology` — TopoStore, arena types, topo types, mesh_cache
- `rustkernel_geom` — AnalyticalGeomStore, SurfaceDef, CurveDef
- `rustkernel_math` — Point3, Vec3, Mat4
- `nalgebra` — rotation/transform math

## Future work (from HITLIST.md)

No immediate backlog items target primitives. They're stable. If new primitive types
are needed (wedge, prism, etc.), follow the existing pattern exactly.

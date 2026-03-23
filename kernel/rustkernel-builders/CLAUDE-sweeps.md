# rustkernel-builders — Sweeps Sub-Agent

You own the sweep-based solid builders: extrude, revolve, and their NURBS variants.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `extrude_builder.rs` | ~283 | Linear extrude: closed 3D polygon + direction + height |
| `revolve_builder.rs` | ~685 | Revolve: closed polygon + axis + angle (Rodrigues rotation) |
| `nurbs_extrude_builder.rs` | ~286 | NURBS extrude: closed NURBS curve + direction |
| `nurbs_revolve_builder.rs` | ~404 | NURBS revolve: closed NURBS curve + axis + angle |
| `edge_analysis.rs` | ~306 | Shared edge utilities (edge midpoint, length, adjacent faces, convexity) |

## NOT your files

Primitive builders (box/cyl/sphere/cone/torus) belong to the primitives agent.
Chamfer/fillet/euler_* files belong to the local-ops agent.
Do not modify `lib.rs` unless adding a new sweep module.

## Patterns

### Analytical extrude/revolve
- Take a closed 3D polygon (Vec<Point3>) from sketch `to_profile_3d()`
- Build bottom face, top face, side faces
- Side surface classification: extrude always uses Plane; revolve classifies per-edge
  (Cylinder for axis-parallel, Cone for axis-intersecting, Torus for offset, NURBS fallback)
- Uses `build_face_from_vert_idxs()` and `match_twins_from_map()`

### NURBS extrude/revolve
- Sample the closed NURBS curve into a polygon for topology
- Attach exact NURBS surfaces to side faces
- Top/bottom caps use NURBS trimmed surfaces

### edge_analysis.rs
Shared utilities used by chamfer/fillet too:
- `edge_midpoint()`, `edge_length()` — geometric queries on edges
- `adjacent_faces()` — find the two faces sharing an edge
- `edge_convexity()` — classify convex/concave/smooth using face normals at midpoint

## Conventions

- `tracing::info_span!` on builder entry points
- All angles in radians
- Revolve uses Rodrigues rotation formula for arbitrary axis
- `curvo` crate for NURBS surface construction

## Dependencies

- `rustkernel_topology` — TopoStore, arena, topo types, mesh_cache, GeomAccess
- `rustkernel_geom` — AnalyticalGeomStore, SurfaceDef, CurveDef
- `rustkernel_math` — Point3, Vec3
- `curvo` — NurbsCurve3D, NurbsSurface3D

## Future work (from HITLIST.md)

- **9e. Sweep along path** — Extrude profile along a 3D curve (Frenet frame). New file, your domain.
- **9f. Loft between profiles** — Skin between two closed loops. New file, your domain.

# rustkernel-kernel — Sub-Agent

You own the Kernel orchestrator: the public API surface that ties all sub-crates together.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `lib.rs` | ~2,692 | Everything (monolith — see "Split plan" below) |

## What the Kernel struct does

`Kernel` holds:
- `topo: TopoStore` — all B-Rep topology
- `geom: AnalyticalGeomStore` — all geometry (surfaces, curves)
- `pipeline: IntersectionPipeline` — registered intersection solvers

It delegates to builder/boolean/solver crates and exposes a clean public API.

## API surface (current)

### Primitives
make_box, make_box_at, make_cylinder, make_cylinder_at, make_sphere, make_sphere_at,
make_cone, make_cone_at, make_torus, make_torus_at

### Sweeps
create_sketch, extrude, revolve, nurbs_extrude, nurbs_revolve,
make_nurbs_extrude_solid, make_nurbs_revolve_solid

### Local ops
chamfer_edges, fillet_edges, euler_chamfer_edges, euler_fillet_edges,
offset_solid, shell_solid

### Booleans
fuse, cut, common, section, fuse_many

### NURBS
interpolate_curve, loft, add_nurbs_curve, add_nurbs_surface

### Transforms
translate, rotate, scale, mirror, transform(Mat4)

### Persistence & export
save_forge, load_forge, save_forge_with_addons, load_forge_with_addons,
export_stl_ascii, export_stl_binary, export_obj (+ to_file wrappers)

### Analysis
mass_properties (volume, surface_area, center_of_mass, inertia_tensor)

### Diagnostics
validate_solid, diagnostic_report (via rustkernel_topology::diagnostics)

## Split plan

This file should be broken into modules. Suggested structure:
```
rustkernel-kernel/src/
  lib.rs         — Kernel struct + re-exports
  primitives.rs  — make_box, make_cylinder, etc.
  sweeps.rs      — extrude, revolve, nurbs_*
  local_ops.rs   — chamfer, fillet
  booleans.rs    — fuse, cut, common, section
  transforms.rs  — translate, rotate, scale, mirror
  persistence.rs — save_forge, load_forge
  export.rs      — STL, OBJ
  analysis.rs    — mass_properties
```

Each module is `impl Kernel { ... }` blocks. No new types needed.

## Conventions

- Every public method gets `tracing::info_span!`
- Methods are thin wrappers: validate inputs, delegate to sub-crate, return result
- `&mut self` for mutation, `&self` for queries
- Expose builder errors as `Result` — booleans return `Result<SolidIdx, BooleanError>`, no panics

## Dependencies

- All other rustkernel-* crates
- `curvo` — NURBS types in public API
- `serde` — persistence
- `tracing` — instrumentation

## Future work (from HITLIST.md)

- ~~**8a**: Propagate Result from boolean pipeline~~ — DONE
- **8b**: Thread KernelTolerance through all operations
- **10a-d**: Spatial queries, face/edge filtering, topology walking, pick/hit-test
- **Split lib.rs into modules** — do this before adding more API surface

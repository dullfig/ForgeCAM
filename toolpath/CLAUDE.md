# forgecam-toolpath

2D/3D toolpath computation for the ForgeCAM CAM system.

## Key Dependency: cavalier_contours

Use `cavalier_contours` (https://github.com/jbuckmccready/cavalier_contours) for all
2D contour operations. It provides:

- **2D polyline offset with native arc segments** — arcs stay as arcs through offset,
  not approximated as line segments. Critical for G-code quality (G02/G03 output).
- **Boolean operations** (union, intersection, difference) on closed polylines —
  needed for pocket/contour region computation.
- **Multi-polyline offset** for shapes with holes/islands.
- **Winding number** point-in-polygon, containment testing.
- **Area, path length, closest point** on polylines with arc segments.

Core type is `Polyline` with `PlineVertex` (x, y, bulge). Bulge encodes arc segments:
0 = line, nonzero = arc (bulge = tan(angle/4)).

```toml
[dependencies]
cavalier_contours = "0.3"  # check crates.io for latest
forgecam-geometry = { path = "../geometry" }
```

## Dependencies

- `forgecam-geometry` — 3D surface queries, projection, curve evaluation
- `cavalier_contours` — 2D contour offset and booleans (native arc segments)
- Do NOT depend on the kernel directly for toolpath math. Read the solid model
  through a shared interface, then work in geometry-space.

## Scope

This crate handles:
- Cutter compensation (2D offset via cavalier_contours)
- Pocket clearing (contour-parallel, zigzag)
- Contour/profile toolpaths
- 3D surface toolpaths (project onto surfaces via forgecam-geometry)
- Lead-in / lead-out / linking moves
- Gouge detection
- Step-over calculation (scallop height from surface curvature)

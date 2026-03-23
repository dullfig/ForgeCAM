# ForgeCAM

A full-featured CAM system built in Rust, from the geometry engine up.

ForgeCAM is a ground-up reimagining of what a modern CAM system should be — built for
aerospace job shops that need reliable STEP import, model-based definition, and toolpath
generation they can trust. No legacy baggage, no geometry/solid split, no bolt-on features.

## Architecture

```
geometry/              Pure computational geometry — curves, surfaces, intersection, meshing
kernel/                B-Rep solid modeling kernel — topology, booleans, builders, Euler ops
mrsev/                 Machining feature recognition + closure surface generation
toolpath/              2D/3D toolpath computation — cutter comp, pocketing, contouring
annotations/           Dimensions, GD&T, PMI — semantic + presentation (AP242-aligned)
document/              Per-part document model — multi-document session, single active editor
parameter-registry/    Parametric feature tree, persistent entity IDs, associativity
param-cascade/         Hierarchical parameter dictionary with inheritance and bulk edit
translators/           STEP/IGES import/export, geometry repair, translation validation
gui/                   Native Win32 shell + wgpu viewports — Windows-only
```

### Dependency graph

```
geometry                         ← foundation, zero ForgeCAM deps
    ↑
kernel                           ← imports geometry
    ↑
mrsev                            ← imports kernel + geometry
toolpath                         ← imports geometry only (not kernel)
annotations                      ← imports geometry only (not kernel)
translators                      ← imports kernel + geometry + annotations
    ↑
parameter-registry               ← imports kernel + annotations
    ↑
document                         ← imports kernel + annotations + toolpath + mrsev + param-registry
    ↑
gui                              ← imports everything
```

Strict acyclic — enforced by convention and Cargo. Toolpath and annotations intentionally
do not depend on the kernel. They work in geometry-space, referencing model entities through
persistent IDs rather than kernel types.

## Current Status

| Crate | Status |
|-------|--------|
| **geometry** | Complete — curves, surfaces, intersection, meshing, 2D contour booleans |
| **kernel** | Functional — 290+ tests, booleans, fillets, chamfers, sweeps, Euler ops, sketcher |
| **mrsev** | Phase 1 — type skeleton and stubs |
| **param-cascade** | Complete — hierarchical parameter dictionary, zero deps |
| **parameter-registry** | Phase 1 — NamingJournal + EntityRef types |
| **document** | Skeleton — Document/Session structure |
| **gui** | Framework — Slint UI markup, wgpu viewport stubs |
| **toolpath** | Spec only |
| **annotations** | Spec only |
| **translators** | Spec only |

### Kernel capabilities

The kernel is a workspace of 9 sub-crates implementing a half-edge B-Rep solid modeler:

- **Primitives**: box, cylinder, sphere, cone, torus
- **Sweeps**: linear extrude, revolve (analytical + NURBS)
- **Booleans**: union, intersection, subtraction — works on curved surfaces (cylinder, sphere, torus, NURBS)
- **Local ops**: constant and variable-radius fillet, chamfer (Euler-operator based)
- **Sketcher**: 2D parametric constraint solver (Newton-Raphson + SVD)
- **Topology**: arena-indexed half-edge structure, Euler operators (SEMV/MEF/KEV/KEF)
- **Persistence**: shape evolution tracking for persistent entity IDs

### Geometry engine

Built on [nalgebra](https://nalgebra.org/) and [curvo](https://github.com/mattatz/curvo) (NURBS).
Analytical curves and surfaces with full NURBS support. Surface-surface intersection,
curve-surface intersection, mesh generation. 2D contour booleans via
[cavalier_contours](https://github.com/jbuckmccready/cavalier_contours) — arcs stay as arcs,
critical for G-code quality.

## Design Decisions

**Geometry as separate crate.** Everything below trimmed-surface topology is pure computational
geometry. CAM and annotations import geometry directly without pulling in the kernel. Modeled
after SMLib's NLib+GSNLib layering.

**Arena + index topology.** Cache-friendly, no reference cycles, no garbage collection.
Euler operators guarantee manifold validity by construction.

**Negative radius = flipped normal.** Convention across geometry and kernel for inside-out
surfaces (cylinder bores, sphere cavities). Cone uses negative half-angle.

**Annotations own presentation.** The annotation crate computes 2D lines, arcs, and text
layout (DimensionPresentation). The GUI just renders what it's given. Matches STEP AP242's
split between semantic and presentation PMI.

**Native Windows GUI.** Win32 for application chrome, wgpu for 3D viewports. Target market
is aerospace job shops — CMMC/ITAR environments where everything runs on Windows workstations.
No cross-platform, no Electron.

**Agentic AI integration.** ForgeCAM exposes its kernel as an
[AgentOS](https://github.com/dullfig/agentos) WASM tool component. An LLM can create geometry,
run booleans, recognize features, and inspect parts — all running locally (ITAR-safe with Ollama
on local GPU). The WASM component is a thin dispatch layer; all state lives in the native host.

## Building

Requires Rust stable (2021 edition or later).

```bash
# Build the geometry crate
cd geometry && cargo build

# Build and test the kernel
cd kernel && cargo test

# Run the kernel viewer (dev tool)
cd kernel && cargo run -p rustkernel-viewer
```

## License

Proprietary. All rights reserved.

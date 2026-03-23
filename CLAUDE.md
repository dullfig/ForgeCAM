# ForgeCAM — Architect

You are the top-level architect for ForgeCAM, a full-featured CAM system built in Rust.
Each subdirectory is an independent crate with its own Claude agent working from its own
CLAUDE.md. Your job is system-level coordination: architecture decisions, dependency flow,
cross-crate API contracts, and delegating work to the right sub-agent.

## Project Structure

```
ForgeCAM/
├── geometry/           Pure computational geometry (curves, surfaces, intersection, meshing)
├── kernel/             B-Rep solid modeling kernel (topology, booleans, builders, Euler ops)
├── document/           Single-document model — owns all per-part state
├── toolpath/           2D/3D toolpath computation (cutter comp, pocketing, contouring)
├── annotations/        Dimensions, GD&T, PMI, visual annotation + presentation layout
├── mrsev/              Machining feature recognition + closure surface generation
├── gui/                GUI (renders everything, dumb renderer for annotations)
├── parameter-registry/ Parametric feature tree, persistent entity IDs, associativity
├── translators/        STEP/IGES import/export, geometry repair, translation validation
└── CLAUDE.md           ← you are here
```

## Dependency Graph (strict — no cycles)

```
geometry                         (leaf — no ForgeCAM deps)
    ↑
kernel                           (imports geometry)
    ↑
mrsev                            (imports kernel + geometry)
toolpath                         (imports geometry, NOT kernel)
annotations                      (imports geometry, NOT kernel)
translators                      (imports kernel + geometry + annotations)
    ↑
parameter-registry               (imports kernel + annotations)
    ↑
document                         (imports kernel + annotations + toolpath + mrsev + parameter-registry)
    ↑
gui                              (imports document + everything above)
```

External integration:
```
AgentOS (separate repo: agentos/)
  └── tools/forgecam-tool/       WASM component (thin dispatch)
        imports → kernel-host     backed by rustkernel-kernel (native)
        imports → document-host   backed by forgecam-document (native)
        imports → mrsev-host      backed by forgecam-mrsev (native)
  └── wit/forgecam.wit           WIT interface contract
```

Key rules:
- `geometry` depends on NOTHING in ForgeCAM. It is the foundation.
- `toolpath` does NOT depend on kernel. It reads solid model data through a shared
  interface, then works in geometry-space.
- `annotations` does NOT depend on kernel. References model entities via persistent IDs
  (EntityRef), not kernel types.
- `document` owns per-part state (`Document`) and the session (`Session`). Multiple
  documents can be open; one is active for editing, others run background computation.
- `gui` is native Win32 shell + wgpu viewports. It borrows the Session/Document. It does
  NOT compute annotation layout, toolpath geometry, or feature recognition. Those crates
  emit presentation-ready data. Windows-only — no cross-platform.

## Sub-Agent CLAUDE.md Status

| Crate | CLAUDE.md | API spec | Cargo.toml | Code |
|-------|-----------|----------|------------|------|
| geometry | Yes | API.md (complete) | Yes | Stubs only |
| kernel | Has HITLIST.md | Has API_SPECIFICATION.md | Yes | 278 tests, 9 sub-crates, functional |
| toolpath | Yes | No | No | Empty |
| annotations | Yes | In CLAUDE.md | No | Empty |
| mrsev | Yes | In CLAUDE.md | Yes | Phase 1 (stubs) |
| document | Yes | In CLAUDE.md | Yes | Skeleton |
| gui | Yes | In CLAUDE.md | No | Empty |
| parameter-registry | Yes | No | No | Empty |
| translators | No | No | No | Empty |

## Architecture Decisions (settled)

### Geometry as separate crate
Modeled after SMLib's NLib+GSNLib layers. Everything below trimmed-surface topology is
pure computational geometry. Everything at/above is kernel. CAM and annotations import
geometry directly without pulling in the kernel.

### Kernel architecture
- Arena + index pattern for B-Rep topology (cache-friendly, no reference cycles)
- Dual representation: exact B-Rep + persistent background mesh
- Euler operators (SEMV/MEF/KEV/KEF) for local topology surgery
- GeomAccess trait bridges topology IDs to geometry objects
- AnalyticalGeomStore holds Vec<SurfaceDef>/Vec<CurveDef> from geometry crate

### Annotations own presentation
Annotation crate computes 2D lines/arcs/text (DimensionPresentation). GUI just renders.
Matches STEP AP242 PMI model (semantic + presentation representations).

### cavalier_contours for 2D contour ops
Used by toolpath crate. 2D polyline offset with native arc segments — arcs stay as arcs,
critical for G-code G02/G03 quality. It's a cargo dependency, not a ForgeCAM component.

### Document model: multi-document session, single active editor
Multiple documents can be open (tabs), but only one is active for editing. Others may run
background computation (MRSEV, toolpaths). `Document` struct owns per-part state (kernel,
annotations, features, toolpaths, metadata). `Session` struct holds all open documents.
No cross-document references — documents are fully independent. GUI borrows the Session.

### Native Windows GUI: Win32 shell + wgpu viewports
Windows-only. No cross-platform (target market is aerospace job shops, CMMC/ITAR
environments). Win32 for application chrome (menus, toolbars, docking panels, tree views).
wgpu for 3D viewport rendering (DX12/Vulkan backend). Toolpath simulation uses dexel/Z-buffer
material removal, not real booleans. Not Bevy, not egui, not Tauri/Electron.

### Negative radius = flipped normal
Convention across geometry and kernel: Cylinder/Sphere/Torus use negative radius for
inward-facing normals. Cone uses negative half_angle. Do not change — kernel depends on it.

### Agentic AI: ForgeCAM as AgentOS Tool
ForgeCAM's kernel operations are exposed as an AgentOS WASM tool component. AgentOS owns
the agent runtime (LLM client, state machine, permissions, TUI/chat). ForgeCAM provides
the domain tools (solid modeling, feature recognition, part inspection).

Architecture:
```
LLM → AgentOS agent loop → forgecam.wasm → host imports → native Kernel
```

The WASM component (`agentos/tools/forgecam-tool/`) is a thin dispatch layer. It parses
JSON tool input, calls host-imported functions, and formats results for the LLM. All state
(Kernel, Document, Session) lives in the native host — zero-copy, no serialization overhead
for geometry data.

WIT interface (`agentos/wit/forgecam.wit`) defines three host import interfaces:
- `kernel-host`: solid queries, primitive builders, booleans, transforms, local ops
- `document-host`: document metadata, material, stock
- `mrsev-host`: feature recognition, machining sequence suggestions

ITAR-safe: Ollama runs locally on RTX 5090, part geometry never leaves the workstation.
Optional Anthropic cloud fallback for non-ITAR shops.

## Open Architecture Problems

### 1. Persistent Entity IDs
Annotations, parameter-registry, and mrsev all need to reference "this edge" or "this face"
in a way that survives kernel operations (chamfer, fillet, boolean). No solution chosen yet.
Options: topological naming (TNaming from OCCT), persistent hash, user-assigned tags.
This is the hardest unsolved problem in the system. Likely lives in parameter-registry
with a shared ID type exported for all crates.

### 2. Document Model — SETTLED
Multi-document session with single active editor. `Document` owns per-part state.
`Session` holds all open documents with tab-style switching. One active for editing,
others can run background computation (MRSEV, toolpaths). No cross-document references.
GUI borrows the Session. Lives in `document/` crate.

### 3. STEP Import Pipeline
Translators crate needs to: parse STEP → build geometry → attempt B-Rep assembly →
diagnose failures → offer repair. The repair workflow is interactive (needs GUI).
Pipeline: translators parses → kernel assembles → geometry repairs → gui shows what broke.

### 4. mrsev ↔ toolpath Interface
mrsev identifies features and generates closure surfaces. Toolpath needs to consume those
features (pocket walls, floor, entry faces, stock shape) to generate cutting paths.
The interface between them is TBD — probably a shared FeatureVolume type.

## Shared Version Constraints

These dependency versions MUST match across all crates:
- `nalgebra = "0.34"` (curvo requires this specific major version)
- `curvo = "0.1.81"` (NURBS engine)
- `serde = "1"` (all public types must be Serialize + Deserialize)
- `tracing = "0.1"` (instrumentation — add to every crate except pure-math leaves)

## Conventions (all crates)

1. All angles in radians
2. Tolerances as function parameters, not globals
3. `tracing` instrumentation on new modules (debug_span! on expensive ops, warn! on failures)
4. Serde derives on all public types
5. No `unsafe` in application code (okay in FFI boundaries if needed)
6. Tests for every public function

## User Context

Dan is a machinist / programmer with aerospace job shop background (Mastercam power user).
Pain points: STEP/IGES import failures, geometry repair, MBD/PMI translation fidelity.
Practical perspective: features that save shop floor time matter more than theoretical
elegance. If something works for 95% of parts, ship it and handle the edge cases later.

## How to Delegate

When work needs to happen in a sub-crate:
1. Check that crate's CLAUDE.md for context and conventions
2. Ensure the work doesn't violate the dependency graph above
3. If the work crosses crate boundaries (new shared type, API change), design the
   interface here first, then update both crates' CLAUDE.md files
4. If an open architecture problem (persistent IDs, document model) is blocking progress,
   flag it rather than making an ad-hoc decision that locks in the wrong approach

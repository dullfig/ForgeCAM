# forgecam-translators

STEP import/export, geometry repair, and translation validation.

## What this crate does

### Phase 1 — STEP Import (the priority)

Parse a STEP Part 21 file → build geometry → assemble B-Rep in the kernel →
diagnose failures → offer repair. This is the critical path — can't machine a
part you can't import.

### Phase 2 — Geometry Repair

Gap healing, normal consistency, sliver face removal, degenerate edge cleanup.
Interactive — the GUI shows what's wrong and the user fixes incrementally.

### Phase 3 — STEP Export

Serialize kernel B-Rep + annotations PMI to STEP AP242. Must include both
semantic and presentation PMI representations.

### Phase 4 — IGES Import (legacy)

IGES files are notoriously broken. Lower priority — most shops have moved to
STEP, but some legacy data (pre-2000) is still IGES-only.

### Phase 5 — Translation Validation

Deviation report between source file and imported model (for MBD compliance).

## Architecture

```
STEP file (.stp/.step)
    │
    ▼
┌──────────────────────┐
│  Parser              │  ISO 10303-21 text format → entity graph
│  (step_parser module) │  Handles AP203, AP214, AP242
└──────────┬───────────┘
           │  HashMap<u64, StepEntity>
           ▼
┌──────────────────────┐
│  Geometry Builder    │  STEP entities → forgecam-geometry types
│  (geom_builder)      │  B-spline curves/surfaces, planes, cylinders, etc.
└──────────┬───────────┘
           │  Vec<SurfaceDef>, Vec<CurveDef>, Vec<Point3>
           ▼
┌──────────────────────┐
│  Topology Assembler  │  STEP topology entities → kernel B-Rep
│  (topo_assembler)    │  Closed_shell → Shell, Advanced_face → Face, etc.
└──────────┬───────────┘
           │  SolidIdx (or ImportResult with diagnostics)
           ▼
┌──────────────────────┐
│  Repair Engine       │  Diagnose and fix: gaps, flips, slivers, degen edges
│  (repair)            │  Interactive — returns issues list, user picks fixes
└──────────────────────┘
```

## STEP Entity Coverage (Phase 1 minimum)

### Geometry entities (→ forgecam-geometry types)

| STEP Entity | Our Type | Priority |
|-------------|----------|----------|
| CARTESIAN_POINT | Point3 | Must |
| DIRECTION | Vec3 | Must |
| AXIS2_PLACEMENT_3D | (origin, axis, ref_dir) | Must |
| PLANE | SurfaceDef::Plane | Must |
| CYLINDRICAL_SURFACE | SurfaceDef::Cylinder | Must |
| CONICAL_SURFACE | SurfaceDef::Cone | Must |
| SPHERICAL_SURFACE | SurfaceDef::Sphere | Must |
| TOROIDAL_SURFACE | SurfaceDef::Torus | Must |
| B_SPLINE_SURFACE_WITH_KNOTS | SurfaceDef::NurbsSurface | Must |
| LINE | CurveDef::LineSegment | Must |
| CIRCLE | CurveDef::Circle | Must |
| ELLIPSE | CurveDef::Ellipse | Must |
| B_SPLINE_CURVE_WITH_KNOTS | CurveDef::NurbsCurve | Must |
| SURFACE_OF_REVOLUTION | SurfaceDef (convert to analytical or NURBS) | Phase 2 |
| SURFACE_OF_LINEAR_EXTRUSION | SurfaceDef (convert to analytical or NURBS) | Phase 2 |
| OFFSET_SURFACE | Evaluate + approximate as NURBS | Phase 2 |

### Topology entities (→ kernel B-Rep)

| STEP Entity | Kernel Operation | Priority |
|-------------|-----------------|----------|
| CLOSED_SHELL | Shell (one per solid) | Must |
| ADVANCED_FACE | Face (surface_id + loops) | Must |
| FACE_BOUND / FACE_OUTER_BOUND | Loop (outer vs inner) | Must |
| EDGE_LOOP | Loop of half-edges | Must |
| ORIENTED_EDGE | HalfEdge (with direction flag) | Must |
| EDGE_CURVE | Edge (curve + vertex endpoints) | Must |
| VERTEX_POINT | Vertex (point_id) | Must |

### Product structure (metadata)

| STEP Entity | What we extract |
|-------------|-----------------|
| PRODUCT_DEFINITION | Part name, ID |
| SHAPE_REPRESENTATION | Units, coordinate system |
| REPRESENTATION_CONTEXT | Length unit (mm vs inch) |
| APPLICATION_PROTOCOL_DEFINITION | AP203/AP214/AP242 detection |

## Parser Strategy

Use the **`ruststep`** crate (or a minimal hand-written parser if `ruststep`
doesn't handle our files). STEP Part 21 format is line-oriented:

```
#1 = CARTESIAN_POINT('', (0.0, 0.0, 0.0));
#2 = DIRECTION('', (0.0, 0.0, 1.0));
#3 = AXIS2_PLACEMENT_3D('', #1, #2, #4);
#4 = DIRECTION('', (1.0, 0.0, 0.0));
#5 = PLANE('', #3);
```

Each line is `#id = ENTITY_TYPE(parameters);` where parameters can be
literals, references to other entities (`#id`), or lists. We need:
1. Tokenizer: split into entity records
2. Entity parser: parse each record into typed struct
3. Graph builder: resolve `#id` references into a navigable entity graph

## Key Design Decisions

### Units
STEP files specify units (mm, inch, radian, degree). We must detect and
convert on import. Kernel uses inches internally. Convert all geometry
during import — don't carry unit metadata through the kernel.

### Coordinate Systems
STEP may use a different coordinate system. Detect from
AXIS2_PLACEMENT_3D and transform all geometry during import.

### Error Handling
Import must NEVER panic. Return `ImportResult`:
```rust
pub struct ImportResult {
    /// Successfully imported solids.
    pub solids: Vec<SolidIdx>,
    /// Warnings (non-fatal issues that were auto-repaired).
    pub warnings: Vec<ImportWarning>,
    /// Errors (entities that couldn't be imported).
    pub errors: Vec<ImportError>,
    /// Entities that were skipped (unsupported type).
    pub skipped: Vec<(u64, String)>,
}
```

### Tolerance
STEP geometry often has gaps at edge boundaries (vertices don't exactly
match curve endpoints, curves don't exactly lie on surfaces). We need a
sewing tolerance (default ~1e-4 for mm files, ~4e-6 for inch files).

## Test Corpus

940 STEP files available for testing:
- Q:\001 -- PRINTS\ (537 files, 2013-2026, 5KB to 107MB, aerospace)
- M:\ (403 files)
- Mix of AP203, AP214, AP242
- Real-world files from CATIA, NX, SOLIDWORKS, Creo origins

## Dependencies

```toml
[dependencies]
forgecam-geometry = { path = "../geometry" }
rustkernel-kernel = { path = "../kernel/rustkernel-kernel" }
rustkernel-topology = { path = "../kernel/rustkernel-topology" }
rustkernel-geom = { path = "../kernel/rustkernel-geom" }
rustkernel-math = { path = "../kernel/rustkernel-math" }
serde = { version = "1", features = ["derive"] }
tracing = "0.1"
```

## Module Structure

```
src/
  lib.rs              Public API: import_step(), ImportResult
  step_parser/
    mod.rs            Tokenizer + entity record parser
    entities.rs       STEP entity type definitions
    graph.rs          Entity reference resolution
  geom_builder.rs     STEP geometry → forgecam-geometry types
  topo_assembler.rs   STEP topology → kernel B-Rep
  repair/
    mod.rs            Repair engine API
    gap_heal.rs       Edge gap healing
    normal_fix.rs     Normal consistency
    sliver.rs         Sliver face detection + removal
  export/
    mod.rs            STEP AP242 export (Phase 3)
  iges/
    mod.rs            IGES import (Phase 4)
  validation.rs       Translation deviation report (Phase 5)
```

## Implementation Priority

1. Parser — read entity graph from .stp file
2. Geometry builder — convert STEP geometry entities to our types
3. Topology assembler — build kernel solid from STEP topology
4. Basic import test — load a simple STEP box, render in GUI
5. Test against real files — iterate on the 940-file corpus
6. Repair — gap healing, normal fixing
7. Export, IGES, validation — later phases

## Status

Architecture designed. No code yet.

# forgecam-mrsev

Material Removal Shape Element Volume (MRSEV) recognition for ForgeCAM.

Based on the MRSEV framework from Kramer, Gupta, and Nau (1991-1994).
Reference papers in `docs/`.

## What this crate does

1. **Feature recognition** — Traverses the delta volume (stock minus part), classifies
   faces by surface type, constructs maximal parametric MRSEV volumes (holes, pockets,
   grooves, edge cuts, ramps, rotation pockets). Air is allowed in MRSEVs.

2. **Model generation** — Finds irredundant set covers of the delta volume. Proposes
   default machining sequences grouped by orientation, large-to-small. User reorders
   in GUI; system recalculates metrics and warns about accessibility issues.

3. **Closure surface generation** — Creates cap surfaces to "plug" feature openings,
   reconstructing the stock volume above each feature for toolpath consumption.

4. **Toolpath bridge** — `FeatureVolume` type contains everything toolpath needs to
   generate cutting paths without importing the kernel.

## Module structure

```
src/
  lib.rs              Crate root, re-exports, MrsevError, FeatureVolume/FeatureType
  volume.rs           MrsevVolume enum and all parametric types
  feature.rs          RecognizedFeature, FeatureSet, FeatureId
  recognize.rs        Recognition algorithm (face traversal -> features)
  model.rs            MrsevModel, set cover, model generation, sequencing
  evaluate.rs         Metric computation (time, setups, tool changes)
  closure.rs          Closure surface generation
  access.rs           AccessDirection, accessibility analysis
```

## Dependencies

```toml
rustkernel-math       # Point3, Vec3
rustkernel-topology   # TopoStore, SolidIdx, FaceIdx, SurfaceKind, GeomAccess
rustkernel-geom       # AnalyticalGeomStore, SurfaceDef
rustkernel-builders   # edge_analysis (adjacency, convexity, dihedral_angle)
rustkernel-kernel     # Kernel struct (for boolean ops)
nalgebra = "0.34"
serde = "1"
tracing = "0.1"
thiserror = "2"
```

## Key types

### MrsevVolume (volume.rs)
Enum with variants: `Hole`, `Pocket`, `Groove`, `EdgeCut`, `Ramp`, `RotationPocket`,
`BRepVolume`. Each is a closed parametric solid with planar top surface.

### Profile2D (volume.rs)
Ordered sequence of `ProfileSegment::Line` / `ProfileSegment::Arc` with bulge convention
matching cavalier_contours (`bulge = tan(angle/4)`).

### RecognizedFeature (feature.rs)
MRSEV volume + source faces + access direction + effective volume + through flag.

### FeatureSet (feature.rs)
All primary MRSEVs found in one recognition pass. Point-in-time snapshot.

### MrsevModel (model.rs)
Ordered machining sequence + ModelMetrics. User-reorderable.

### FeatureVolume (lib.rs)
Bridge to toolpath: volume + access direction + closure + depth + corner radii.

### ClosureSurface (closure.rs)
Cap surface (SurfaceDef from kernel) + trim boundary (Profile2D).

## Public API

### Recognition
```rust
fn recognize_features(topo, geom, stock, part, config) -> Result<FeatureSet, MrsevError>
```

### Model generation
```rust
fn generate_models(feature_set, topo, geom, config) -> Result<Vec<MrsevModel>, MrsevError>
fn recalculate(model, feature_set, topo, geom) -> Result<(), MrsevError>
fn propose_sequence(feature_ids, feature_set) -> Vec<SequenceEntry>
```

### Evaluation
```rust
fn estimate_feature_time(feature, mrr) -> f64
fn compute_metrics(sequence, feature_set, mrr, setup_time, tool_change_time) -> ModelMetrics
```

### Closure surfaces
```rust
fn generate_closures(model, feature_set, topo, geom) -> Result<Vec<ClosureSurface>, MrsevError>
```

### Accessibility
```rust
fn check_accessibility(feature, workpiece, direction, tool_radius, topo, geom) -> Result<bool, MrsevError>
fn find_access_directions(feature, workpiece, tool_radius, topo, geom) -> Result<Vec<AccessDirection>, MrsevError>
```

### Toolpath bridge
```rust
fn to_feature_volumes(model, feature_set, closures) -> Vec<FeatureVolume>
```

## Kernel APIs consumed

| API | Source crate |
|-----|-------------|
| `TopoStore` (face/edge/vertex traversal) | `rustkernel-topology` |
| `AnalyticalGeomStore` (surface/curve queries) | `rustkernel-geom` |
| `GeomAccess` trait (eval, normal, kind) | `rustkernel-topology` |
| `SurfaceKind`, `CurveKind` enums | `rustkernel-topology` |
| `SurfaceDef`, `CurveDef` | `rustkernel-geom` |
| `SolidIdx`, `FaceIdx`, `EdgeIdx` | `rustkernel-topology` |
| `edge_adjacency()`, `edge_convexity()`, `dihedral_angle()` | `rustkernel-builders::edge_analysis` |
| `Kernel::cut()` (boolean subtraction) | `rustkernel-kernel` |

## Design decisions

- **Shape != operation**: MRSEVs describe geometry only; machining strategy is toolpath's job
- **Air allowed**: Features can overlap air, avoiding complex booleans at intersections
- **Snapshot model**: Recognition is point-in-time; re-run after model changes. No persistent ID dependency
- **User-driven sequencing**: System proposes; user reorders; system recalculates and warns
- **BRepVolume escape hatch**: Arbitrary B-Rep for sculptured features NURBS can't parametrize

## Implementation phases

| Phase | What | Status |
|-------|------|--------|
| 1 | All types + crate skeleton + stubbed functions | **Done** |
| 2 | Hole recognition (cylindrical/conical face traversal) | Not started |
| 3 | Pocket recognition (planar floor, wall walking, islands) | Not started |
| 4 | Groove + edge cut recognition + pattern detection | Not started |
| 5 | Model generation (set cover, sequencing, metrics) | Not started |
| 6 | Closure surfaces + toolpath bridge + ramp/rotation pocket | Not started |

## Conventions

- All angles in radians
- Tolerances as function parameters, not globals
- Negative radius = flipped normal (inherited from kernel)
- `tracing` instrumentation: `debug_span!` on expensive ops, `warn!` on failures
- All public types: `Debug + Clone + Serialize + Deserialize`

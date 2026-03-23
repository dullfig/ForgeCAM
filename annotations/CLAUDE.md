# forgecam-annotations

Dimensions, GD&T, PMI, and visual annotation for the ForgeCAM CAM system.

## What this crate is

This crate owns two things:
1. **Semantic annotation data** — what is being measured, the value, tolerances, datum references
2. **Presentation generation** — computing the 2D lines, arcs, arrows, and text placement
   that represent an annotation on screen or on a drawing sheet

The GUI does not compute annotation layout. This crate generates a `DimensionPresentation`
(line segments, arcs, text placements, pick regions) and the GUI just renders it.

This matches how production CAD systems work and how STEP AP242 stores PMI — both a
semantic representation (what the tolerance means) and a presentation representation
(where the lines go).

## Dependencies

```toml
[dependencies]
forgecam-geometry = { path = "../geometry" }   # 2D line/arc types, projection math
serde = { version = "1", features = ["derive"] }
```

Does NOT depend on the kernel. Annotations reference kernel entities (edges, faces, vertices)
via stable persistent IDs, not by importing kernel types directly. The ID scheme is TBD —
probably lives in a shared types crate or the parameter-registry.

## Module structure

```
src/
  dimension.rs       — Linear, angular, radial, diametral, ordinate dimensions
  gdt.rs             — Feature control frames, datum references, composite tolerances
  surface_finish.rs  — Surface finish symbols (Ra, Rz, etc.)
  notes.rs           — Text notes, balloons, flags
  presentation.rs    — Layout engine: semantics + view → lines/arcs/text
  pmi.rs             — AP242 PMI data model for STEP import/export
```

## Semantic types

```rust
/// Reference to a model entity (edge, face, vertex) via persistent ID.
/// The ID scheme must survive kernel operations (chamfer, fillet, boolean).
pub struct EntityRef { /* TBD — persistent naming / stable ID */ }

pub enum DimensionType {
    Linear { from: EntityRef, to: EntityRef },
    Angular { edge_a: EntityRef, edge_b: EntityRef },
    Radial { arc_or_circle: EntityRef },
    Diametral { arc_or_circle: EntityRef },
    Ordinate { from_datum: EntityRef, to: EntityRef, axis: Axis },
}

pub struct Dimension {
    pub dim_type: DimensionType,
    pub nominal: f64,
    pub tolerance: Option<Tolerance>,
    pub text_override: Option<String>,
    pub placement: AnnotationPlane,  // which plane/view the dim lives in
}

pub enum Tolerance {
    Symmetric(f64),                    // ± value
    Bilateral { plus: f64, minus: f64 },
    Limits { upper: f64, lower: f64 },
    FitClass(String),                  // e.g. "H7/g6"
}
```

## Presentation output

```rust
pub struct DimensionPresentation {
    pub lines: Vec<LineSegment2D>,     // witness lines, dimension line, leaders
    pub arcs: Vec<Arc2D>,              // arc arrows, radius leaders
    pub text: Vec<TextPlacement>,      // string + position + height + angle
    pub pick_regions: Vec<Aabb2D>,     // for mouse hit testing in GUI
}

pub struct TextPlacement {
    pub content: String,
    pub position: Point2,   // anchor point
    pub height: f64,        // text height in annotation-space units
    pub angle: f64,         // rotation in radians
    pub anchor: TextAnchor, // center, left, right
}
```

The presentation engine handles:
- Witness line extension rules and minimum gap from geometry
- Arrow-inside vs arrow-outside based on available space
- Tolerance stacking and text formatting
- Leader routing for offset dimensions
- Multiple views (same dimension can present differently in 3D view vs 2D drawing sheet)

## GD&T / Feature Control Frames

```rust
pub struct FeatureControlFrame {
    pub characteristic: GdtCharacteristic,  // position, profile, flatness, etc.
    pub tolerance_value: f64,
    pub modifier: Option<MaterialCondition>,  // MMC, LMC, RFS
    pub datums: Vec<DatumRef>,
    pub attached_to: EntityRef,
}

pub enum GdtCharacteristic {
    Flatness, Straightness, Circularity, Cylindricity,        // form
    Parallelism, Perpendicularity, Angularity,                 // orientation
    Position, Concentricity, Symmetry,                         // location
    CircularRunout, TotalRunout,                                // runout
    ProfileOfLine, ProfileOfSurface,                            // profile
}
```

## Key design decisions

1. **Annotations don't own geometry, they reference it.** An annotation holds an EntityRef,
   not a copy of the edge. If the model changes, the annotation follows (or invalidates)
   based on the persistent ID resolution.

2. **Presentation is computed, not stored.** The semantic data is the source of truth.
   Presentation is regenerated when the view changes, the model changes, or display
   settings change. Exception: AP242 import may store presentation directly when
   semantic data isn't available.

3. **Associativity is someone else's problem.** When a chamfer deletes an edge that a
   dimension references, the annotation crate can detect the broken reference (EntityRef
   no longer resolves). But deciding what to do about it — re-attach, warn, delete —
   is a session/document-model concern, likely in parameter-registry or a layer above.

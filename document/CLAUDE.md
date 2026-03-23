# forgecam-document

Document model and session management for ForgeCAM.

## What this crate does

1. **`Document`** — per-part state container. Owns the kernel, annotations, feature tree,
   toolpaths, MRSEV results, and metadata. The unit of save/load.

2. **`Session`** — holds all open documents. One is **active** (being edited, rendered by
   GUI). Others may be running background computation (feature recognition, toolpath
   generation) or idle. Tab-style switching between documents.

## What this crate is NOT

- Not a GUI. The GUI borrows the session/document and renders.
- Not a feature tree. The parameter-registry crate owns the feature tree; the document
  holds an instance of it.
- Not a kernel. The kernel is a field on the document.

## Architecture: Multi-Document Session, Single Active Editor

Design decision (settled): Multiple documents can be open simultaneously, but only one
is active for editing at a time. Background computation (MRSEV, toolpaths) runs on
inactive documents. The user can switch tabs to check progress.

This handles the real workflow: start MRSEV on part A, switch to part B to set it up,
switch back to A when recognition finishes. Aerospace job shops run 3-5 jobs at a time.

Key rules:
- **One active document** — the GUI edits/renders this one
- **Others are background** — running computation or showing results
- **No cross-document references** — documents are fully independent
- **Documents are independently saveable** — each is its own .forge file

## Document struct

```rust
pub struct Document {
    pub metadata: PartMetadata,       // name, units, material, stock
    pub kernel: Kernel,               // B-Rep solid modeling state
    // pub annotations: AnnotationSet, // when annotations crate is ready
    // pub features: FeatureTree,      // when parameter-registry is ready
    // pub toolpaths: ToolpathSet,     // when toolpath crate is ready
    // pub mrsev: FeatureMap,          // when mrsev has results to store
}
```

## Session struct

```rust
pub struct Session {
    slots: Vec<DocumentSlot>,   // document + compute status
    active: Option<usize>,      // which slot is being edited
}
```

`ComputeStatus` tracks per-document background state: `Idle`, `Recognizing`,
`GeneratingToolpaths`, `Ready`, `Failed(String)`.

## Dependencies

```toml
rustkernel-math       # Point3, Vec3
rustkernel-topology   # TopoStore types
rustkernel-geom       # AnalyticalGeomStore
rustkernel-kernel     # Kernel struct
nalgebra = "0.34"
serde = "1"
tracing = "0.1"
thiserror = "2"
```

Future dependencies (add when crates are ready):
- `forgecam-annotations`
- `forgecam-toolpath`
- `forgecam-mrsev`
- `forgecam-parameter-registry`

## Module structure

```
src/
  lib.rs          Document struct, constructors
  metadata.rs     PartMetadata, Units, Material, StockShape
  session.rs      Session, DocumentId, DocumentSlot, ComputeStatus
```

Future modules:
- `persistence.rs` — .forge file format, save/load, versioned migration
- `undo.rs` — undo/redo stack (command pattern over document mutations)
- `dirty.rs` — change tracking ("unsaved changes" flag, notification to GUI)

## Key design decisions

1. **Multiple documents open, one active.** Tab-style UX. Background computation on
   inactive documents. No cross-document references.
2. **Document owns everything by value.** Kernel, annotations, features, toolpaths
   are fields, not Arc/Rc pointers. The document is the single owner of its state.
3. **Serialize = save.** The entire Document derives Serialize/Deserialize (once Kernel
   supports it). File format is the serialized Document (.forge file).
4. **GUI borrows, doesn't own.** The GUI holds a reference to the Session, renders the
   active document, shows tab bar with status for all documents.
5. **ComputeStatus is observable.** GUI polls or subscribes to status changes to update
   the tab bar (spinner, checkmark, error icon).

## Conventions

- All public types: `Debug + Clone + Serialize + Deserialize` (except Document and Session
  which skip Clone — you don't clone the world)
- Units stored in metadata, not implicit. Functions that care about units read
  `document.metadata.units`.
- Stock shape is optional — some workflows start from a STEP import with no stock defined.
- DocumentId is a simple u32 wrapper, stable for the lifetime of the session.

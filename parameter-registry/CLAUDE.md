# forgecam-parameter-registry

Parametric feature tree, persistent entity IDs, and associativity.

## What this crate does

1. **Persistent entity IDs** — `EntityRef` (opaque u64, lives in geometry crate) +
   `NamingJournal` that maps EntityRef ↔ live FaceIdx/EdgeIdx/VertexIdx. The journal
   processes kernel `ShapeEvolution` records to maintain stable identity across operations
   and feature-tree replays.

2. **Feature tree** — (Not yet implemented) Ordered list of modeling operations that can
   be replayed to regenerate the model.

3. **Associativity** — (Not yet implemented) Change propagation, broken reference detection.

## Architecture

```
geometry          ← EntityRef type (opaque u64, usable by all crates)
    ↑
kernel            ← ShapeEvolution (populated by all ops)
    ↑
parameter-registry ← NamingJournal (maps EntityRef ↔ FaceIdx, processes evolution)
    ↑
document          ← (future) Owns the journal, calls process_evolution after every op
```

### Key types

- **`EntityRef(u64)`** — in `geometry/src/entity_ref.rs`. Top 2 bits = EntityKind
  (Face/Edge/Vertex), lower 62 bits = monotonic counter. Copy + Hash + Eq + Serialize.

- **`FeatureId(u32)`** — persistent identifier for a feature in the tree.

- **`NamingKey`** — deterministic name derived from creation context:
  - `Primitive { feature, kind, ordinal }` — new entity from a builder
  - `FromEdge { feature, source_edge_key }` — fillet/chamfer face
  - `FromIntersection { feature, face_a_key, face_b_key }` — boolean seam
  - `SplitFrom { feature, parent_key, split_ordinal }` — one piece of a split face

- **`NamingJournal`** — runtime mapping engine:
  - `process_evolution(feature_id, &ShapeEvolution)` — core method, call after every op
  - `resolve_face/edge/vertex(EntityRef) -> ResolveResult` — Live / Split / Deleted / Unknown
  - `entity_for_face/edge/vertex(Idx) -> Option<EntityRef>` — reverse lookup
  - `clear()` — reset live state for replay (preserves key_to_ref for stable IDs)

### Identity rules

- **Modified/CopiedFrom** entities inherit their parent's EntityRef (no new key).
  This is why a dimension on a face survives a chamfer on an adjacent edge.
- **SplitFrom** entities get new EntityRefs (distinct from parent and from each other).
  Split ordering uses sorted FaceIdx raw values for determinism.
- **Deleted** entities resolve to `ResolveResult::Deleted { by_feature }`.

## Dependencies

```toml
[dependencies]
forgecam-geometry = { path = "../geometry" }
rustkernel-topology = { path = "../kernel/rustkernel-topology" }
serde = { version = "1", features = ["derive"] }
tracing = "0.1"
```

## Key files

- `src/lib.rs` — module structure, re-exports
- `src/naming.rs` — NamingKey, FeatureId, NamingJournal, ResolveResult, all tests

## Status

Phase 1 complete: NamingJournal + EntityRef implemented with 8 passing tests.
Next steps: feature tree, document integration, MRSEV/annotations migration.

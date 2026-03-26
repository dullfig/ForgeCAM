# param-cascade

Hierarchical machining parameter dictionary with folder inheritance and lock/protect.

## The Problem

In Mastercam (and most CAM), every operation stores its own private copy of every
parameter. Want to change the stepover on 15 pocket operations? Edit each one. Want to
know which operations leave .010" on the bottom? Open them one by one and check.

NX does better: operations live in folders, folders define default parameters, children
inherit. But even NX's model is limited.

## The Solution

A **parameter data dictionary** with **cascading inheritance** and **lock/protect**.

### Core concepts

1. **ParamDef** — a canonical parameter definition (name, type, unit, valid range).
   The dictionary of "what parameters exist." Every toolpath that has a depth-of-cut
   uses the SAME `depth_of_cut` key, not its own bespoke field.

2. **OpTree** — a tree of folders and operations. Folders group operations (ROUGHING,
   FINISH_CUTS, DRILL_OPS). Operations are leaves.

3. **Cascade** — parameters resolve by walking up the tree:
   operation → parent folder → grandparent folder → global defaults.
   First defined value wins.

4. **Three param states at any node:**
   - **Inherit** (not set locally) — value comes from ancestor cascade
   - **Override** (set locally) — explicit value, bulk edit CAN change it
   - **Locked** (set locally + locked) — explicit value, bulk edit SKIPS it

5. **Bulk edit** — select N operations, set a param. All selected ops get the value
   EXCEPT locked ones (skipped, reported back to caller).

6. **Query** — "which ops have `bottom_stock > 0.005`?" resolves every op's value
   through the cascade, then filters.

### The "clone rough to finish" workflow

```
ROUGHING/                     (side_stock=0.010, bottom_stock=0.010, stepover=0.500)
  ├─ OP1 adaptive pocket      (inherits all)
  ├─ OP2 adaptive pocket      (inherits all)
  └─ OP3 contour              (side_stock=0.015 🔒)  ← customer callout

  [User copies ROUGHING → FINISH_CUTS, changes two folder params]

FINISH_CUTS/                  (side_stock=0.000, bottom_stock=0.000, stepover=0.030)
  ├─ OP4 adaptive pocket      (inherits → now 0.000/0.000/0.030)
  ├─ OP5 adaptive pocket      (inherits → now 0.000/0.000/0.030)
  └─ OP6 contour              (side_stock=0.015 🔒, rest inherits)
```

OP6 kept its locked side stock. Everything else cascaded. Zero manual editing.

### AI companion integration

The dictionary makes natural-language queries trivial:

- "Where am I leaving stock on the bottom?" → query `bottom_stock > 0`
- "Change all roughing stepovers to 60% of tool diameter" → bulk edit on ROUGHING folder
- "Which operations use more than 0.250 depth of cut?" → query `depth_of_cut > 0.250`
- "Set all finish passes to climb milling" → bulk edit `cut_direction = Climb` on FINISH_CUTS

The AI doesn't need to understand toolpath internals — just the dictionary keys and
the cascade query API.

## Dependency Position

```
param-cascade               (imports NOTHING from ForgeCAM — pure data structure)
    ↑
document                    (imports param-cascade, owns the OpTree per-part)
    ↑
gui                         (renders the tree, dispatches edits)
toolpath                    (reads resolved params for path computation)
```

This crate has **zero ForgeCAM dependencies**. It's a generic hierarchical parameter
cascade engine. The machining-specific parameter definitions (what keys exist, their
types and defaults) are defined by the consumer (document or toolpath crate), not here.

## External Dependencies

```toml
[dependencies]
serde = { version = "1", features = ["derive"] }
tracing = "0.1"
```

That's it. No geometry, no kernel, no nalgebra.

## Core Types

```rust
/// Canonical parameter key — interned string for fast comparison
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ParamKey(pub String);

/// What kind of value a parameter holds
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ValueType {
    F64,
    Bool,
    Int,
    Enum(Vec<String>),   // named choices: "Climb", "Conventional"
    String,
}

/// A concrete parameter value
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ParamValue {
    F64(f64),
    Bool(bool),
    Int(i64),
    Enum(String),
    String(String),
}

/// Units for display and conversion
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ParamUnit {
    Length,          // inches or mm (document decides)
    Percent,
    Speed,          // SFM or m/min
    FeedRate,       // IPM or mm/min
    FeedPerTooth,   // IPT or mm/tooth
    SpindleRpm,
    Angle,          // degrees (display) — stored as radians internally
    Count,
    None,
}

/// Definition of a parameter in the dictionary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamDef {
    pub key: ParamKey,
    pub display_name: String,
    pub description: String,
    pub unit: ParamUnit,
    pub value_type: ValueType,
    pub default: ParamValue,
    pub range: Option<(f64, f64)>,  // min/max for numeric types
}

/// State of a parameter at a tree node
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParamState {
    Inherit,    // no local value — cascade from parent
    Override,   // local value — bulk edit CAN change it
    Locked,     // local value — bulk edit SKIPS it
}

/// A parameter value at a node
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamEntry {
    pub value: ParamValue,
    pub state: ParamState,  // Override or Locked (Inherit nodes have no entry)
}

/// Index into the node arena
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NodeId(pub u32);

/// What kind of tree node
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NodeKind {
    Folder,
    Operation(String),  // operation type: "adaptive_pocket", "contour", "drill", etc.
}

/// A node in the operation tree (folder or operation)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpNode {
    pub id: NodeId,
    pub name: String,
    pub kind: NodeKind,
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,          // non-empty only for Folder
    pub params: HashMap<ParamKey, ParamEntry>,  // locally set params only
}

/// Result of resolving a parameter through the cascade
#[derive(Debug, Clone)]
pub struct ResolvedParam {
    pub value: ParamValue,
    pub source: NodeId,         // which node provided the value
    pub source_name: String,    // "ROUGHING" or "OP3" — for display
    pub is_locked: bool,        // locked at the resolved source
    pub is_default: bool,       // fell through to global default
}

/// Result of a bulk edit operation
#[derive(Debug, Clone)]
pub struct BulkEditResult {
    pub updated: Vec<NodeId>,   // nodes that were changed
    pub skipped: Vec<NodeId>,   // nodes that were locked (not changed)
}

/// Declares which dictionary params a toolpath type uses
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpTypeRegistration {
    pub op_type: String,                  // "adaptive_pocket", "contour", "drill"
    pub display_name: String,             // "Adaptive Pocket"
    pub required_params: Vec<ParamKey>,   // must resolve to a value before posting
    pub optional_params: Vec<ParamKey>,   // may be set, but not required
}

/// The parameter dictionary + operation tree
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpTree {
    dictionary: HashMap<ParamKey, ParamDef>,
    op_types: HashMap<String, OpTypeRegistration>,  // registered op types
    nodes: Vec<OpNode>,
    roots: Vec<NodeId>,  // top-level folders/ops
    next_id: u32,
}
```

## Public API

```rust
impl OpTree {
    // ── Construction ──────────────────────────────────────────────
    pub fn new() -> Self;

    // ── Dictionary management ─────────────────────────────────────
    pub fn register_param(&mut self, def: ParamDef);
    pub fn param_def(&self, key: &ParamKey) -> Option<&ParamDef>;
    pub fn all_param_defs(&self) -> impl Iterator<Item = &ParamDef>;

    // ── Op type registration ──────────────────────────────────────
    /// Register a toolpath type and the dictionary keys it uses.
    /// Toolpath code calls this at startup — "I'm an adaptive pocket,
    /// I need stepover, depth_of_cut, side_stock, bottom_stock, ..."
    pub fn register_op_type(&mut self, reg: OpTypeRegistration);
    /// What params does this op type require/support?
    pub fn op_type_params(&self, op_type: &str) -> Option<&OpTypeRegistration>;
    /// All registered op types
    pub fn all_op_types(&self) -> impl Iterator<Item = &OpTypeRegistration>;
    /// Which op types use a given param? (reverse lookup)
    /// "Who uses stepover?" → ["adaptive_pocket", "pocket", "facing"]
    pub fn op_types_using_param(&self, key: &ParamKey) -> Vec<&str>;
    /// Validate an operation node: are all required params resolvable?
    pub fn validate_op(&self, node: NodeId) -> Vec<ValidationError>;

    // ── Tree structure ────────────────────────────────────────────
    pub fn add_folder(&mut self, name: &str, parent: Option<NodeId>) -> NodeId;
    pub fn add_operation(&mut self, name: &str, op_type: &str, parent: Option<NodeId>) -> NodeId;
    pub fn remove_node(&mut self, id: NodeId);  // removes children too
    pub fn move_node(&mut self, id: NodeId, new_parent: Option<NodeId>);
    pub fn clone_subtree(&mut self, id: NodeId, new_parent: Option<NodeId>) -> NodeId;
    pub fn node(&self, id: NodeId) -> Option<&OpNode>;
    pub fn children(&self, id: NodeId) -> &[NodeId];
    pub fn roots(&self) -> &[NodeId];

    // ── Parameter access ──────────────────────────────────────────
    /// Set a param on a node (Override state by default)
    pub fn set_param(&mut self, node: NodeId, key: ParamKey, value: ParamValue);
    /// Set a param and lock it
    pub fn set_param_locked(&mut self, node: NodeId, key: ParamKey, value: ParamValue);
    /// Lock an existing param (must already have a local value)
    pub fn lock_param(&mut self, node: NodeId, key: &ParamKey);
    /// Unlock a param (keeps the value, changes state to Override)
    pub fn unlock_param(&mut self, node: NodeId, key: &ParamKey);
    /// Clear a local param (reverts to Inherit)
    pub fn clear_param(&mut self, node: NodeId, key: &ParamKey);
    /// Get the local entry (None if inheriting)
    pub fn local_param(&self, node: NodeId, key: &ParamKey) -> Option<&ParamEntry>;

    // ── Cascade resolution ────────────────────────────────────────
    /// Resolve a param for a node: walk up tree, fall back to default
    pub fn resolve(&self, node: NodeId, key: &ParamKey) -> Option<ResolvedParam>;
    /// Resolve ALL params for a node (full effective parameter set)
    pub fn resolve_all(&self, node: NodeId) -> HashMap<ParamKey, ResolvedParam>;

    // ── Bulk operations ───────────────────────────────────────────
    /// Set a param on multiple nodes. Locked nodes are skipped.
    pub fn bulk_set(&mut self, nodes: &[NodeId], key: ParamKey, value: ParamValue) -> BulkEditResult;
    /// Set a param on a folder + all descendant ops. Locked nodes skipped.
    pub fn cascade_set(&mut self, folder: NodeId, key: ParamKey, value: ParamValue) -> BulkEditResult;

    // ── Query ─────────────────────────────────────────────────────
    /// Find all leaf operations where a resolved param matches a predicate
    pub fn query<F>(&self, key: &ParamKey, pred: F) -> Vec<NodeId>
        where F: Fn(&ParamValue) -> bool;
    /// Find all leaf operations where a resolved param equals a value
    pub fn query_eq(&self, key: &ParamKey, value: &ParamValue) -> Vec<NodeId>;
    /// Find all leaf operations where a resolved numeric param is in a range
    pub fn query_range(&self, key: &ParamKey, min: f64, max: f64) -> Vec<NodeId>;
}

/// Validation error for an operation node
#[derive(Debug, Clone)]
pub enum ValidationError {
    /// Required param has no value anywhere in the cascade
    MissingRequired { key: ParamKey, op_type: String },
    /// Param value doesn't match the dictionary's ValueType
    TypeMismatch { key: ParamKey, expected: ValueType, got: ParamValue },
    /// Param value is outside the dictionary's valid range
    OutOfRange { key: ParamKey, value: f64, min: f64, max: f64 },
    /// Op uses a param key not in the dictionary (custom/unknown)
    UnknownParam { key: ParamKey },
}
```

## Standard Machining Parameters

These are defined by the consumer (document crate), not by param-cascade itself.
But here's the standard set for reference:

### Cut parameters
| Key | Display Name | Unit | Type | Typical Default |
|-----|-------------|------|------|----------------|
| `stepover` | Stepover | Length | F64 | 0.500 |
| `stepover_pct` | Stepover % | Percent | F64 | 50.0 |
| `depth_of_cut` | Depth of Cut | Length | F64 | 0.250 |
| `final_cut` | Final Cut | Length | F64 | 0.010 |
| `side_stock` | Side Stock to Leave | Length | F64 | 0.010 |
| `bottom_stock` | Bottom Stock to Leave | Length | F64 | 0.010 |
| `cut_direction` | Cut Direction | None | Enum | "Climb" |

### Feeds and speeds
| Key | Display Name | Unit | Type |
|-----|-------------|------|------|
| `spindle_rpm` | Spindle Speed | SpindleRpm | F64 |
| `feed_rate` | Feed Rate | FeedRate | F64 |
| `plunge_rate` | Plunge Rate | FeedRate | F64 |
| `retract_rate` | Retract Rate | FeedRate | F64 |
| `surface_speed` | Surface Speed | Speed | F64 |
| `feed_per_tooth` | Feed Per Tooth | FeedPerTooth | F64 |

### Heights / clearance
| Key | Display Name | Unit | Type |
|-----|-------------|------|------|
| `clearance_height` | Clearance Height | Length | F64 |
| `retract_height` | Retract Height | Length | F64 |
| `feed_height` | Feed Height | Length | F64 |
| `top_of_stock` | Top of Stock | Length | F64 |
| `depth` | Total Depth | Length | F64 |

### Coolant / safety
| Key | Display Name | Unit | Type |
|-----|-------------|------|------|
| `coolant` | Coolant | None | Enum |
| `compensation` | Cutter Compensation | None | Enum |
| `compensation_direction` | Comp Direction | None | Enum |

## Material Integration

Material is NOT a parameter — it's a **context** that determines parameter defaults.
The cascade has an additional root above the folder tree:

```
Material (Brass 360)
  → Tool defaults (1/2 Carbide EM in Brass: SFM=600, FPT=0.004)
    → Folder (ROUGHING: DOC=0.5D, WOC=65%)
      → Operation (OP1 pocket: inherits all)
```

### How material change works

"Change the material from 6061-T6 to Brass 360":

1. **Resolve new defaults** — look up each tool's cutting data for Brass.
   This comes from the tool library's `CuttingData` table (tool + material → SFM, FPT, etc.).
2. **Inject into cascade root** — the material defaults become the bottom-most
   fallback in the cascade, below the folder tree.
3. **Cascade resolves** — every operation's feeds/speeds are re-resolved.
   Folder overrides survive. Operation overrides survive. Locked values survive.
4. **Warn on gaps** — if a tool has no cutting data for Brass, flag it:
   "Tool T3 (1/4 HSS drill) has no data for Brass 360 — using aluminum defaults."
5. **Warn on mismatches** — if a tool isn't suitable for the material:
   "Tool T7 (HSS rougher) is not recommended for Brass — consider carbide."
6. **Re-post** — toolpath regenerates with new feeds/speeds, post outputs new G-code.

### What changes on material switch

| Parameter | Typically changes? | Source |
|-----------|-------------------|--------|
| `surface_speed` | Yes | tool-library CuttingData |
| `feed_per_tooth` | Yes | tool-library CuttingData |
| `spindle_rpm` | Yes (derived from SFM + diameter) | computed |
| `feed_rate` | Yes (derived from RPM × FPT × flutes) | computed |
| `depth_of_cut` | Maybe | tool-library CuttingData |
| `stepover` / `stepover_pct` | Maybe | tool-library CuttingData |
| `coolant` | Maybe (Ti needs through-spindle) | tool-library CuttingData |
| `cut_direction` | Rarely | unchanged |
| `side_stock` / `bottom_stock` | No | unchanged |
| `clearance_height` etc. | No | unchanged |

### The "almost" in "almost just change the material and re-post"

Three things that prevent fully automatic material change:

1. **Missing data** — tool has no cutting data entry for the new material.
   The cascade uses the previous material's values and warns. The machinist
   must add cutting data or accept the defaults.

2. **Tool unsuitability** — HSS tools in hardened steel, aluminum-specific
   geometry in titanium. The tool library can flag this if tools have
   `material_compatibility` tags, but it's advisory, not blocking.

3. **Strategy change** — roughing aluminum is different from roughing titanium
   (chip thinning, axial depth limits, dwell sensitivity). The CutScript
   sequence might need to change, not just the parameters. This is the one
   thing the cascade CAN'T automate — it requires the programmer's judgment.

### Material definition

Materials live in a simple lookup table (separate from param-cascade):

```rust
pub struct Material {
    pub id: MaterialId,
    pub name: String,                    // "6061-T6 Aluminum"
    pub category: String,                // "Aluminum", "Steel", "Titanium"
    pub hardness: Option<f64>,           // HRC or HB
    pub machinability_rating: Option<f64>, // relative to B1112 steel = 100%
    pub notes: String,
}
```

The tool library stores `CuttingData` entries keyed by `(ToolId, MaterialId)`.
The cascade queries the tool library to get material-specific defaults.

## Module Structure

```
src/
  lib.rs          Re-exports, OpTree constructor
  types.rs        ParamKey, ParamValue, ParamDef, ParamUnit, ValueType, ParamState, ParamEntry
  tree.rs         OpTree, OpNode, NodeId, NodeKind — structure and mutation
  registry.rs     OpTypeRegistration, register/lookup/reverse-lookup, validation
  cascade.rs      Resolution logic (walk-up), resolve/resolve_all
  bulk.rs         Bulk edit and cascade_set operations
  query.rs        Query engine (predicate, equality, range)
```

## Design Decisions

1. **Zero ForgeCAM dependencies.** This is a pure data structure crate. Machining-specific
   definitions live in the consumer. Could theoretically be used for any hierarchical
   parameter system.

2. **Arena storage** (Vec<OpNode> with NodeId indices). Consistent with kernel's pattern.
   No reference cycles, cache-friendly, serializable.

3. **Inherit = absence.** A node that inherits a param simply doesn't have it in its
   local HashMap. Resolution walks up. This means cloning a subtree and moving it to a
   new folder automatically re-inherits from the new parent — no bookkeeping.

4. **Lock protects against bulk/cascade edits only.** Direct `set_param` on a locked
   param is allowed (explicit intent). `bulk_set` and `cascade_set` skip locked params.
   This matches the UX: you can always manually override a lock if you go to that
   specific operation, but mass operations respect it.

5. **ParamKey is a String, not an enum.** The dictionary is extensible — users, plugins,
   or the AI companion can register new parameter types at runtime. The standard keys
   are conventional, not compiler-enforced.

6. **No unit conversion in this crate.** The crate stores values in document units
   (inches or mm). Conversion between inch/metric is the document crate's job.

## Conventions

- All public types: `Debug, Clone, Serialize, Deserialize`
- `tracing` instrumentation on bulk operations (debug_span! on cascade_set, warn! on lock skips)
- Tests for every public method
- No unsafe

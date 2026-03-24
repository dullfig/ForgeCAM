# CutScript — First Draft

A domain-specific language for describing machining sequences. Captures the
strategy a senior programmer uses to go from recognized features to finished
part — stock allowances, sequencing, approach strategies, tool selection,
and shop floor experience that would otherwise live in someone's head.

## Design Goals

1. **Readable by machinists** — looks like shop talk, not code
2. **Writable by LLMs** — structured enough that an AI can generate valid programs
3. **Composable** — sequences can be saved, shared, and applied to similar features
4. **Parameterized** — stock, tools, speeds are values, not hardcoded
5. **Feature-aware** — references mrsev features and EntityRef persistent IDs

## File Extension

`.cuts`

## Example: Simple Pocket

```
# Rough and finish a pocket
# Applied to: pocket with floor + walls

sequence rough_finish_pocket(pocket: Pocket) {

    # Phase 1: Hog out material
    rough pocket.volume leave 0.010 {
        tool flat 0.500
        stepover 65%
        ramp spiral angle 3
        feed 40.0
        speed 6000
    }

    # Phase 2: Finish the floor, staying away from walls
    finish pocket.floor avoid pocket.walls by 0.050 {
        tool flat 0.500
        stepover 30%
        climb
    }

    # Phase 3: Finish the walls
    finish pocket.walls {
        tool flat 0.500
        climb
        spring_pass
    }

    # Phase 4: One last cleanup on the floor corners
    finish pocket.floor {
        tool flat 0.500
        stepover 30%
        climb
    }
}
```

## Example: Drilling

```
sequence drill_holes(holes: HoleGroup) {

    # Center drill first
    spot holes {
        tool spot_drill 0.500
        depth 0.100
        speed 4000
        feed 10.0
    }

    # Peck drill
    drill holes {
        tool drill 0.250
        peck 0.200 retract 0.050
        speed 3000
        feed 8.0
    }

    # Ream to size
    ream holes {
        tool reamer 0.2500
        speed 800
        feed 4.0
        dwell 0.5 at_bottom
    }
}
```

## Example: Contour with Depth Passes

```
sequence profile_outside(contour: Contour, stock: Stock) {

    rough contour leave 0.015 {
        tool flat 0.375
        depth_per_pass 0.250
        stepover 50%
        conventional    # rough conventional for chip load
        lead_in arc radius 0.100
        lead_out arc radius 0.100
    }

    finish contour {
        tool flat 0.375
        depth_per_pass full
        climb           # finish climb for surface quality
        lead_in arc radius 0.100
        lead_out arc radius 0.100
        spring_pass
    }
}
```

## Language Elements

### Operations

The core verbs — what the tool is doing:

| Operation | Meaning |
|-----------|---------|
| `rough`   | Remove bulk material, leave stock |
| `finish`  | Cut to final dimension |
| `semi`    | Intermediate pass (e.g., leave 0.002 for grind) |
| `spot`    | Center drill / spot drill |
| `drill`   | Hole drilling (peck, through, blind) |
| `ream`    | Ream to size |
| `bore`    | Boring bar operation |
| `tap`     | Thread tapping |
| `engrave` | Shallow marking / engraving |
| `chamfer_op` | Edge break / chamfer milling |

### Geometry References

How you point at what to cut. These reference mrsev features and resolve
through EntityRef to actual geometry:

```
pocket.floor          # bottom face(s) of a pocket
pocket.walls          # side faces of a pocket
pocket.volume         # the material to remove (closure volume)
contour               # a profile boundary
holes                 # a hole group
face(F42)             # explicit EntityRef
edge(E17)             # explicit EntityRef
```

### Stock Control

The critical concept that separates roughing from finishing:

```
leave 0.010           # leave this much stock on all surfaces
leave 0.010 on walls  # leave stock on specific surfaces only
avoid walls by 0.050  # don't cut within this distance of a surface
to_size               # cut to final dimension (default for finish)
```

### Tool Specification

```
tool flat 0.500       # flat endmill, 0.5" diameter
tool ball 0.250       # ball endmill
tool bull 0.500 0.030 # bull nose, diameter + corner radius
tool drill 0.250      # twist drill
tool spot_drill 0.500 # spot/center drill
tool reamer 0.2500    # reamer
tool tap 0.250-20     # tap, diameter-TPI
```

### Strategy

How the tool moves:

```
climb                 # climb milling (conventional = opposite)
conventional

stepover 65%          # WOC as percentage of tool diameter
stepover 0.100        # WOC as absolute value
depth_per_pass 0.250  # DOC per Z level
depth_per_pass full   # full depth in one pass

ramp spiral angle 3   # spiral ramp entry, 3° ramp angle
ramp zigzag angle 5   # zigzag ramp
ramp plunge           # straight plunge (drills, center drills)
ramp helix            # helical entry

lead_in arc radius 0.100   # arc lead-in
lead_in tangent length 0.050
lead_out arc radius 0.100

spring_pass           # repeat last pass with no new depth
                      # (lets tool deflection recover)

peck 0.200 retract 0.050   # peck drilling
dwell 0.5 at_bottom        # pause at bottom (ream, bore)
```

### Feeds and Speeds

```
feed 40.0             # IPM (inches per minute)
speed 6000            # RPM
feed_per_tooth 0.003  # alternative to feed — calculated from RPM + flutes
```

These can also come from param-cascade (material + tool → lookup table),
so they're optional in the script:

```
rough pocket.volume leave 0.010 {
    tool flat 0.500
    # feeds/speeds inherited from param-cascade for this material+tool
}
```

### Control Flow

```
# Conditional — depth-dependent strategy
if pocket.depth > tool.diameter * 2 {
    ramp helix
} else {
    ramp plunge
}

# Repeat — spring passes with count
repeat 2 {
    finish pocket.walls { climb }
}

# For each — iterate over feature groups
for hole in holes {
    drill hole { ... }
}
```

### Sequences and Composition

```
# Named sequence = reusable recipe
sequence my_pocket_recipe(p: Pocket) { ... }

# Apply to a specific feature
apply rough_finish_pocket to pocket1
apply rough_finish_pocket to pocket2

# Include another sequence
sequence full_part {
    apply rough_finish_pocket to pocket1
    apply drill_holes to hole_group_1
    apply profile_outside to outer_contour, stock
}
```

## Type System

Minimal — just enough to catch mistakes:

| Type | What it is |
|------|-----------|
| `Pocket` | mrsev pocket feature (floor + walls + volume) |
| `Contour` | Profile boundary (open or closed) |
| `HoleGroup` | Set of holes (same diameter) |
| `Face` | Single face reference (EntityRef) |
| `Edge` | Single edge reference |
| `Stock` | Material block definition |
| `Tool` | Tool specification |
| `number` | Dimensional value (inches by default) |
| `percent` | 0-100% value |
| `angle` | Degrees (converted to radians internally) |

## Compilation Target

CutScript compiles to a `MachiningPlan` — a list of `MachiningOp` structs
that the toolpath crate consumes:

```rust
struct MachiningPlan {
    ops: Vec<MachiningOp>,
}

struct MachiningOp {
    kind: OpKind,           // Rough, Finish, Drill, etc.
    geometry: GeometryRef,  // What to cut (feature + faces)
    tool: ToolSpec,         // What to cut with
    strategy: Strategy,     // How to cut (stepover, DOC, ramp, etc.)
    stock: StockControl,    // Leave / avoid / to_size
    feeds: FeedSpeed,       // IPM, RPM (or inherited from param-cascade)
}
```

The toolpath crate turns each `MachiningOp` into actual cutter paths
(polylines/arcs with feeds). The post processor turns those into G-code.

## Open Questions

1. **Units** — inches vs mm? Default to inches with `units mm` override?
   Or require explicit units always?

2. **Tool library** — should tools be defined inline or referenced from
   an external tool library? Probably both: `tool flat 0.500` inline,
   `tool "T12"` from library.

3. **Workholding / fixtures** — does the DSL need to express "flip part"
   or "move to op 2 fixture"? Or is that a higher-level concern?

4. **Multi-axis** — 3-axis first. 4th/5th axis adds coordinate systems,
   tool vectors. Handle later.

5. **Error handling** — what happens when a sequence references a feature
   that doesn't exist after a model change? The EntityRef system handles
   this (ResolveResult::Deleted), but the DSL needs a way to express
   fallbacks or warnings.

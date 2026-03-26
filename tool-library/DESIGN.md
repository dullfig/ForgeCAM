# Tool Library

First-class tool management for ForgeCAM. Every tool is a real object with
geometry, cutting data, and shop-specific metadata — not a form with 30 fields
where half are blank and the other half are wrong.

## The Mastercam Problem

Mastercam's tool library was designed for standard catalog tools. Custom tools
(form tools, step drills, specials, brazed-tip cutters) are second-class
citizens bolted on with "custom profile" fields that don't integrate with
toolpath computation. The library is per-machine-group, tools don't travel
between files cleanly, and the database format has changed multiple times.

The result: machinists maintain their own shadow tool list (spreadsheet, sticky
notes on the monitor, a dog-eared catalog) because the CAM system's library
doesn't reflect what's actually in the crib.

## Design Principles

1. **Tool geometry is real geometry.** A tool isn't a set of parameters — it's
   a 2D profile that gets revolved. Any tool shape is expressible, from a
   standard 4-flute endmill to a custom dovetail cutter to a form tool with
   an arbitrary profile. The profile IS the tool definition.

2. **One library, multiple scopes.** Shop-level (every machine can see it),
   machine-level (tools loaded in this machine's carousel), job-level (tools
   used in this program). Same tool object, different visibility.

3. **Cutting data lives with the tool.** Not in a separate feeds/speeds table.
   A tool knows what it can do in each material: "1/2 endmill in 6061: 8000 RPM,
   40 IPM, 0.250 DOC, 65% WOC." This integrates with param-cascade.

4. **Custom tools are just tools.** No special "custom" flag. Define the profile,
   define the cutting data, done. A standard 1/2 endmill and a custom form tool
   are represented identically — they just have different profiles.

5. **Tools have real-world identity.** Crib number, purchase date, measured
   runout, remaining life, regrind count. This is shop management data, not
   CAM data, but it lives alongside the geometry because that's where people
   look for it.

## Tool Definition

```rust
pub struct Tool {
    // Identity
    pub id: ToolId,
    pub crib_number: Option<String>,     // "T-0347" (shop's tracking)
    pub description: String,             // "1/2 4FL SE CARBIDE ENDMILL"
    pub manufacturer: Option<String>,    // "Kennametal"
    pub part_number: Option<String>,     // "F3AA0500AWL"

    // Geometry — the actual shape
    pub profile: ToolProfile,            // 2D profile (revolved = 3D tool)
    pub flutes: u32,
    pub helix_angle: Option<f64>,        // radians
    pub overall_length: f64,
    pub shank_diameter: f64,

    // Holder
    pub holder: Option<HolderId>,        // reference to holder library
    pub gauge_length: f64,               // stickout from gauge line
    pub gauge_diameter: f64,             // at the gauge line

    // Cutting data per material
    pub cutting_data: Vec<CuttingData>,

    // Shop data
    pub status: ToolStatus,              // Available, InUse, NeedsRegrind, Retired
    pub regrind_count: u32,
    pub notes: String,
}
```

## Tool Profile

The key insight: every tool is a 2D profile revolved around the tool axis.

```rust
pub enum ToolProfile {
    /// Standard parametric tools — defined by a few dimensions.
    FlatEndmill {
        diameter: f64,
        corner_radius: f64,    // 0.0 = sharp, >0 = bull nose
        cut_length: f64,
    },
    BallEndmill {
        diameter: f64,
        cut_length: f64,
    },
    ChamferMill {
        diameter: f64,
        angle: f64,            // half angle from axis (radians)
        tip_diameter: f64,     // 0.0 = pointed
    },
    Drill {
        diameter: f64,
        point_angle: f64,      // full included angle (radians)
    },
    Tap {
        diameter: f64,
        pitch: f64,            // thread pitch
    },

    /// Custom profile — arbitrary 2D contour revolved around axis.
    /// Handles form tools, step drills, specials, dovetails, anything.
    Custom {
        /// 2D polyline/arc profile from centerline to OD, bottom to top.
        /// Uses the same Contour2D type as the geometry crate.
        contour: Vec<ProfileSegment>,
        cut_length: f64,
    },
}

/// A segment of a custom tool profile.
pub enum ProfileSegment {
    Line { to: [f64; 2] },                          // [r, z]
    Arc { to: [f64; 2], center: [f64; 2], cw: bool },
}
```

Standard tools are just shortcuts — `FlatEndmill { diameter: 0.5, corner_radius: 0.0 }`
generates the same profile as a Custom with two line segments. The toolpath engine
sees only the profile; it doesn't care which variant created it.

## Cutting Data

```rust
pub struct CuttingData {
    pub material: MaterialId,      // reference to material library
    pub speed: f64,                // RPM (or SFM, converted at runtime)
    pub feed_per_tooth: f64,       // IPT (feed = RPM * flutes * FPT)
    pub axial_depth: f64,          // max DOC
    pub radial_depth_pct: f64,     // max WOC as % of diameter
    pub ramp_angle: Option<f64>,   // max ramp angle (radians)
    pub plunge_capable: bool,      // can this tool plunge-cut?
    pub coolant: CoolantType,      // flood, mist, through, air, dry
    pub notes: String,
}

pub enum CoolantType {
    Flood,
    Mist,
    ThroughSpindle,
    AirBlast,
    Dry,
}
```

This integrates with param-cascade: the cutting data can inherit from
material defaults, be overridden per tool, and further overridden per
operation. The cascade resolves: material default → tool-specific →
operation-specific.

## Library Scopes

```
Shop Library (central, all machines)
  └── Machine Library (tools loaded in carousel, positions assigned)
        └── Job Library (tools used in this program, with operation assignments)
```

A tool exists once in the shop library. It can appear in multiple machine
libraries (assigned to a carousel position). A job references tools from
the machine library.

## Storage

JSON or MessagePack on disk. One file per library scope:
- `tools.json` — shop library (hundreds of tools)
- `machine-haas-vf4.json` — machine-specific carousel
- Embedded in `.forge` file — job-level tool assignments

## Tool Geometry for Simulation

The tool profile is also used for:
- **Toolpath simulation** — the dexel material removal engine needs the tool
  shape to compute what material is removed at each position
- **Collision detection** — holder + tool assembly checked against part/fixture
- **Toolpath display** — tool marker in the viewport shows actual tool shape

The profile → revolved solid conversion uses the same geometry crate
primitives as any other solid. A tool IS a solid (just a special one).

## Integration Points

```
tool-library
  → toolpath (tool geometry for offset computation)
  → post (tool number, description, gauge length for G-code)
  → param-cascade (cutting data feeds into operation parameters)
  → gui (tool editor panel, carousel visualization)
  → simulation (tool profile for dexel material removal)
```

## Status

Design phase. No code yet.

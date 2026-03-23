//! MRSEV volume types — parametric descriptions of material removal volumes.
//!
//! Each variant is a closed parametric solid with a planar top surface.
//! Air is allowed inside the volume (it need not be entirely within material).

use rustkernel_math::{Point3, Vec3};
use serde::{Deserialize, Serialize};

// ── Primary volume enum ──

/// A Material Removal Shape Element Volume — a closed parametric volume
/// describing material to be removed by a machining operation.
///
/// Each variant is a parametric solid with a planar top surface and a
/// native coordinate system. Air is allowed inside the volume (it need
/// not be entirely within material).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MrsevVolume {
    Hole(Hole),
    Pocket(Pocket),
    Groove(Groove),
    EdgeCut(EdgeCut),
    Ramp(Ramp),
    RotationPocket(RotationPocket),
    /// Escape hatch: arbitrary B-Rep volume when parametric types don't fit.
    /// Toolpath must handle this via general 3D strategies.
    BRepVolume(BRepVolumeRef),
}

// ── Hole ──

/// Cylindrical or conical-bottomed hole.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hole {
    pub location: Point3,
    /// Axis direction (into material).
    pub orientation: Vec3,
    pub radius: f64,
    pub depth: f64,
    pub bottom: HoleBottom,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HoleBottom {
    Flat,
    /// Full included angle, in radians.
    Conical { tip_angle: f64 },
    Spherical { radius: f64 },
}

// ── Pocket ──

/// Linear sweep pocket — a closed profile swept perpendicular to its plane.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pocket {
    pub location: Point3,
    /// Normal to top surface (into material).
    pub orientation: Vec3,
    pub depth: f64,
    /// Outer boundary.
    pub profile: Profile2D,
    pub corner_radius: Option<f64>,
    pub floor_radius: Option<f64>,
    pub islands: Vec<Island>,
    pub is_through: bool,
}

/// A 2D profile: ordered sequence of line segments and circular arcs.
/// Compatible with cavalier_contours polylines.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Profile2D {
    pub segments: Vec<ProfileSegment>,
    pub closed: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ProfileSegment {
    Line { end: [f64; 2] },
    /// bulge = tan(angle/4), cavalier_contours convention.
    Arc { end: [f64; 2], bulge: f64 },
}

/// Island left standing inside a pocket.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Island {
    pub profile: Profile2D,
    /// Height from pocket floor.
    pub height: f64,
    pub kind: IslandKind,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum IslandKind {
    /// Flat-topped island.
    Mesa,
    /// Opposite of a groove.
    AntiGroove,
    /// Opposite of a rotation pocket.
    Rotation { axis: Vec3 },
}

// ── Groove ──

/// Swept cross-section along a profile path.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Groove {
    pub location: Point3,
    pub orientation: Vec3,
    pub cross_section: GrooveCrossSection,
    pub path: Profile2D,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GrooveCrossSection {
    Square { width: f64, depth: f64 },
    Vee { width: f64, depth: f64, angle: f64 },
    Round { radius: f64 },
    Custom(Profile2D),
}

// ── Edge cut ──

/// Chamfer or fillet along an edge path.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeCut {
    pub path: Profile2D,
    pub kind: EdgeCutKind,
    pub orientation: Vec3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EdgeCutKind {
    Chamfer { leg1: f64, leg2: f64 },
    Fillet { radius: f64 },
}

// ── Ramp ──

/// Helical or angled ramp.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ramp {
    pub location: Point3,
    pub orientation: Vec3,
    pub kind: RampKind,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RampKind {
    Straight {
        direction: Vec3,
        depth: f64,
        tool_radius: f64,
    },
    Helical {
        axis: Vec3,
        helix_radius: f64,
        pitch: f64,
        turns: f64,
        tool_radius: f64,
    },
}

// ── Rotation pocket ──

/// Revolved pocket (turning-style feature on a milling part).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RotationPocket {
    pub location: Point3,
    pub axis: Vec3,
    /// Profile revolved around axis.
    pub profile: Profile2D,
    pub kind: RotationKind,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RotationKind {
    /// Axis parallel to spindle.
    Vertical,
    /// Axis perpendicular, clipped by plane.
    Horizontal { cut_plane_offset: f64 },
}

// ── BRep escape hatch ──

/// Reference to an arbitrary B-Rep solid in the kernel when
/// parametric MRSEV types are insufficient.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRepVolumeRef {
    /// Index into kernel's solid arena.
    pub solid_idx: u32,
    pub description: String,
}

// ── Copy patterns ──

/// A complete MRSEV: one or more volumes + optional pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mrsev {
    pub name: String,
    /// Usually 1; >1 = grouped volumes.
    pub volumes: Vec<MrsevVolume>,
    pub pattern: Option<CopyPattern>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CopyPattern {
    Rectangular {
        rows: u32,
        cols: u32,
        row_spacing: f64,
        col_spacing: f64,
    },
    Circular {
        count: u32,
        center: Point3,
        axis: Vec3,
        /// Total angular span in radians.
        angular_span: f64,
    },
    Mirror {
        plane_origin: Point3,
        plane_normal: Vec3,
    },
}

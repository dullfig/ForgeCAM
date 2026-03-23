//! Cutter engagement angle computation for constant-engagement roughing.
//!
//! The core primitive for adaptive clearing / dynamic motion toolpaths (Mastercam
//! Dynamic, Fusion 360 Adaptive, etc.). Given a circular tool and a 2D stock
//! boundary, computes how much of the tool circumference is buried in material.
//!
//! # Algorithm
//!
//! 1. Intersect the tool circle with every segment of the stock boundary contour.
//! 2. Collect all intersection angles on the tool circle.
//! 3. Sort them, then classify each arc between consecutive angles as inside or
//!    outside the stock (via winding number of the arc midpoint).
//! 4. Sum the inside arcs → engagement angle.
//!
//! # References
//!
//! - Ibaraki, Yamaji, Matsubara — "Geometric Algorithms for Computing Cutter
//!   Engagement Functions in 2.5D Milling Operations"
//! - Lukacs et al. — "The Fast Constant Engagement Offsetting Method" (FACEOM), 2019

use std::f64::consts::TAU;

use crate::contour2d::{Contour2D, Segment2D};
use crate::types::*;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Result of an engagement angle computation.
#[derive(Debug, Clone)]
pub struct EngagementResult {
    /// Total engagement angle in radians, in [0, 2π].
    /// 0 means the tool is entirely outside the stock.
    /// 2π means the tool is entirely buried (full slotting).
    pub angle: f64,

    /// Entry and exit points where the tool circle crosses the stock boundary,
    /// ordered CCW around the tool circle. Each pair (entry, exit) defines one
    /// engaged arc.
    pub arcs: Vec<EngagedArc>,
}

/// A contiguous arc of the tool circle that is inside the stock.
#[derive(Debug, Clone)]
pub struct EngagedArc {
    /// Start angle on the tool circle (radians, CCW from +X).
    pub start_angle: f64,
    /// End angle on the tool circle (radians, CCW from +X).
    pub end_angle: f64,
    /// Angular span in radians (always positive).
    pub span: f64,
}

// ---------------------------------------------------------------------------
// Angle helpers
// ---------------------------------------------------------------------------

/// Normalize angle to [0, 2π).
fn normalize(a: f64) -> f64 {
    let mut r = a % TAU;
    if r < 0.0 {
        r += TAU;
    }
    r
}

/// Point on the tool circle at angle θ.
fn circle_point(center: Point2, radius: f64, angle: f64) -> Point2 {
    Point2::new(
        center.x + radius * angle.cos(),
        center.y + radius * angle.sin(),
    )
}

// ---------------------------------------------------------------------------
// Circle–segment intersection (returns angles on the tool circle)
// ---------------------------------------------------------------------------

/// Intersect a full circle (tool) with a line segment. Returns angles on the circle.
fn circle_line_angles(
    center: Point2,
    radius: f64,
    p0: Point2,
    p1: Point2,
) -> Vec<f64> {
    let d = p1 - p0;
    let f = p0 - center;

    let a = d.dot(&d);
    let b = 2.0 * f.dot(&d);
    let c = f.dot(&f) - radius * radius;

    if a.abs() < 1e-30 {
        return vec![];
    }

    let disc = b * b - 4.0 * a * c;
    if disc < -1e-12 {
        return vec![];
    }
    let disc = disc.max(0.0);
    let sqrt_disc = disc.sqrt();

    let mut angles = Vec::new();
    for t in [(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)] {
        if t >= -1e-10 && t <= 1.0 + 1e-10 {
            let t_clamped = t.clamp(0.0, 1.0);
            let p = p0 + d * t_clamped;
            let angle = normalize((p.y - center.y).atan2(p.x - center.x));
            angles.push(angle);
        }
    }

    // Deduplicate tangent case.
    if angles.len() == 2 && (angles[0] - angles[1]).abs() < 1e-10 {
        angles.pop();
    }

    angles
}

/// Intersect a full circle (tool) with a circular arc (boundary segment).
/// Returns angles on the *tool* circle.
fn circle_arc_angles(
    center: Point2,
    radius: f64,
    arc_center: Point2,
    arc_radius: f64,
    arc_start_angle: f64,
    arc_end_angle: f64,
) -> Vec<f64> {
    let d = arc_center - center;
    let dist = d.norm();

    if dist < 1e-15 {
        return vec![]; // concentric
    }
    if dist > radius + arc_radius + 1e-10 {
        return vec![]; // too far
    }
    if dist < (radius - arc_radius).abs() - 1e-10 {
        return vec![]; // one inside the other
    }

    // Standard circle-circle intersection.
    let a_dist = (radius * radius - arc_radius * arc_radius + dist * dist) / (2.0 * dist);
    let h_sq = radius * radius - a_dist * a_dist;
    let h = if h_sq > 0.0 { h_sq.sqrt() } else { 0.0 };

    let unit_d = d / dist;
    let perp = Vec2::new(-unit_d.y, unit_d.x);
    let mid = center + unit_d * a_dist;

    let mut angles = Vec::new();

    let candidates = if h < 1e-12 {
        vec![mid]
    } else {
        vec![mid + perp * h, mid - perp * h]
    };

    // Arc angular span check.
    let arc_span = {
        let mut s = (arc_end_angle - arc_start_angle) % TAU;
        if s <= 0.0 {
            s += TAU;
        }
        s
    };

    for p in candidates {
        // Check if point lies on the boundary arc.
        let angle_on_arc = (p.y - arc_center.y).atan2(p.x - arc_center.x);
        let delta = normalize(angle_on_arc - arc_start_angle);
        if delta <= arc_span + 1e-10 {
            // Point is on the arc — compute angle on the tool circle.
            let angle = normalize((p.y - center.y).atan2(p.x - center.x));
            angles.push(angle);
        }
    }

    angles
}

// ---------------------------------------------------------------------------
// Core: engagement_angle
// ---------------------------------------------------------------------------

/// Compute the engagement angle of a circular cutter against a stock boundary.
///
/// # Arguments
///
/// * `tool_center` — center of the cutter in the XY work plane.
/// * `tool_radius` — radius of the cutter.
/// * `stock` — closed 2D contour representing the current stock boundary.
///   Must be a closed contour with consistent winding (CCW = material inside).
///
/// # Returns
///
/// An [`EngagementResult`] with the total engagement angle and the individual
/// engaged arcs. If the tool is entirely outside the stock, `angle` is 0.
/// If entirely inside (full slotting), `angle` is 2π.
///
/// # Panics
///
/// Panics if `tool_radius <= 0`.
pub fn engagement_angle(
    tool_center: Point2,
    tool_radius: f64,
    stock: &Contour2D,
) -> EngagementResult {
    assert!(tool_radius > 0.0, "tool_radius must be positive");

    let segments = stock.to_segments();
    if segments.is_empty() {
        return EngagementResult {
            angle: 0.0,
            arcs: vec![],
        };
    }

    // Step 1: Collect all intersection angles on the tool circle.
    let mut hit_angles: Vec<f64> = Vec::new();

    for seg in &segments {
        match seg {
            Segment2D::Line(line) => {
                hit_angles.extend(circle_line_angles(
                    tool_center,
                    tool_radius,
                    line.start,
                    line.end,
                ));
            }
            Segment2D::Arc(arc) => {
                hit_angles.extend(circle_arc_angles(
                    tool_center,
                    tool_radius,
                    arc.center,
                    arc.radius.abs(),
                    arc.start_angle,
                    arc.end_angle,
                ));
            }
        }
    }

    // Step 2: Handle degenerate cases.
    if hit_angles.is_empty() {
        // No intersections — tool is entirely inside or entirely outside.
        let wn = stock.winding_number(tool_center);
        if wn != 0 {
            // Tool center inside stock → fully engaged (slotting).
            return EngagementResult {
                angle: TAU,
                arcs: vec![EngagedArc {
                    start_angle: 0.0,
                    end_angle: TAU,
                    span: TAU,
                }],
            };
        } else {
            // Tool entirely outside stock.
            return EngagementResult {
                angle: 0.0,
                arcs: vec![],
            };
        }
    }

    // Deduplicate angles that are very close (multiple segments meeting at a vertex
    // can produce near-duplicate intersection points).
    hit_angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
    hit_angles.dedup_by(|a, b| (*a - *b).abs() < 1e-8);

    // Single intersection point — tangent contact, treat as zero engagement.
    if hit_angles.len() == 1 {
        return EngagementResult {
            angle: 0.0,
            arcs: vec![],
        };
    }

    // Step 3: Walk around the tool circle, testing each arc's midpoint.
    let n = hit_angles.len();
    let mut arcs = Vec::new();
    let mut total_angle = 0.0;

    for i in 0..n {
        let a0 = hit_angles[i];
        let a1 = hit_angles[(i + 1) % n];

        // Angular span of this arc (CCW from a0 to a1).
        let span = if i + 1 < n { a1 - a0 } else { TAU - a0 + a1 };

        // Midpoint angle.
        let mid_angle = a0 + span / 2.0;
        let mid_point = circle_point(tool_center, tool_radius * 0.999, mid_angle);

        // Test if midpoint is inside stock.
        let wn = stock.winding_number(mid_point);
        if wn != 0 {
            arcs.push(EngagedArc {
                start_angle: a0,
                end_angle: if i + 1 < n { a1 } else { a1 + TAU },
                span,
            });
            total_angle += span;
        }
    }

    EngagementResult {
        angle: total_angle.min(TAU),
        arcs,
    }
}

/// Convenience: compute just the engagement angle in radians.
pub fn engagement_angle_simple(
    tool_center: Point2,
    tool_radius: f64,
    stock: &Contour2D,
) -> f64 {
    engagement_angle(tool_center, tool_radius, stock).angle
}

/// Compute the radial depth of cut from tool radius and engagement angle.
///
/// `ae = r * (1 - cos(α/2))` where α is the engagement angle.
pub fn radial_depth_from_engagement(tool_radius: f64, engagement: f64) -> f64 {
    tool_radius * (1.0 - (engagement / 2.0).cos())
}

/// Compute the engagement angle from tool radius and radial depth of cut.
///
/// `α = 2 * acos(1 - ae/r)`
pub fn engagement_from_radial_depth(tool_radius: f64, radial_depth: f64) -> f64 {
    let ratio = (1.0 - radial_depth / tool_radius).clamp(-1.0, 1.0);
    2.0 * ratio.acos()
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// Square stock: 10x10 centered at origin.
    fn square_stock() -> Contour2D {
        Contour2D::from_vertices(
            &[
                (-5.0, -5.0, 0.0),
                (5.0, -5.0, 0.0),
                (5.0, 5.0, 0.0),
                (-5.0, 5.0, 0.0),
            ],
            true,
        )
    }

    #[test]
    fn tool_fully_inside_is_full_slotting() {
        let stock = square_stock();
        let result = engagement_angle(Point2::new(0.0, 0.0), 1.0, &stock);
        assert_relative_eq!(result.angle, TAU, epsilon = 1e-6);
        assert_eq!(result.arcs.len(), 1);
    }

    #[test]
    fn tool_fully_outside_is_zero() {
        let stock = square_stock();
        let result = engagement_angle(Point2::new(20.0, 0.0), 1.0, &stock);
        assert_relative_eq!(result.angle, 0.0, epsilon = 1e-6);
        assert!(result.arcs.is_empty());
    }

    #[test]
    fn tool_half_engaged_on_straight_edge() {
        // Tool centered on the right edge of a 10x10 square.
        // Center at (5, 0), radius 1. Left half of tool is inside stock.
        let stock = square_stock();
        let result = engagement_angle(Point2::new(5.0, 0.0), 1.0, &stock);
        // Engagement should be ~π (half the tool circle is in material).
        assert_relative_eq!(result.angle, PI, epsilon = 0.15);
    }

    #[test]
    fn tool_shallow_engagement() {
        // Tool barely cutting into the right edge.
        // ae = 0.1 (radial depth), r = 1.0
        // engagement = 2 * acos(1 - 0.1/1.0) = 2 * acos(0.9) ≈ 0.902 rad ≈ 51.7°
        let stock = square_stock();
        let expected = engagement_from_radial_depth(1.0, 0.1);
        let result = engagement_angle(Point2::new(5.9, 0.0), 1.0, &stock);
        assert_relative_eq!(result.angle, expected, epsilon = 0.1);
    }

    #[test]
    fn radial_depth_round_trip() {
        let r = 5.0;
        let ae = 1.5;
        let eng = engagement_from_radial_depth(r, ae);
        let ae2 = radial_depth_from_engagement(r, eng);
        assert_relative_eq!(ae, ae2, epsilon = 1e-12);
    }

    #[test]
    fn full_slot_engagement_is_tau() {
        // ae = 2*r (full width) → engagement = 2π
        let eng = engagement_from_radial_depth(5.0, 10.0);
        assert_relative_eq!(eng, TAU, epsilon = 1e-12);
    }

    #[test]
    fn half_diameter_engagement_is_pi() {
        // ae = r → engagement = 2*acos(0) = π
        let eng = engagement_from_radial_depth(5.0, 5.0);
        assert_relative_eq!(eng, PI, epsilon = 1e-12);
    }

    #[test]
    fn concave_corner_increases_engagement() {
        // L-shaped pocket: concave 90° corner at (0, 0).
        // Material is inside the contour (CCW winding).
        //
        //     (0,5)---(-5,5)
        //       |        |
        //     (0,0)---(-5,0)... but we need a concave corner so
        //
        // Use a large square with a rectangular notch cut out.
        // Outer stock: [-10, 10] x [-10, 10]
        // Inner notch: [0, 10] x [0, 10] removed
        // The concave corner is at (0, 0).
        let stock = Contour2D::from_vertices(
            &[
                (-10.0, -10.0, 0.0),
                (10.0, -10.0, 0.0),
                (10.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),   // concave corner
                (0.0, 10.0, 0.0),
                (-10.0, 10.0, 0.0),
            ],
            true,
        );

        // Tool on a straight wall: center at (-1, -9), radius 2.
        // Only the bottom wall contributes. ae = 2 - 1 = 1.
        let straight = engagement_angle(Point2::new(-1.0, -9.0), 2.0, &stock);

        // Tool approaching the concave corner at (0, 0): center at (-1, -1), radius 2.
        // Both the right wall (x=0 for y>0) and top wall (y=0 for x>0) contribute.
        let corner = engagement_angle(Point2::new(-1.0, -1.0), 2.0, &stock);

        assert!(
            corner.angle > straight.angle,
            "concave corner engagement ({:.3} rad, {:.1}°) should exceed \
             straight wall ({:.3} rad, {:.1}°)",
            corner.angle,
            corner.angle.to_degrees(),
            straight.angle,
            straight.angle.to_degrees(),
        );
    }

    #[test]
    fn circular_stock_tool_on_edge() {
        // Circular stock: full-circle arc as contour.
        // cavalier_contours represents a circle as 2 vertices with bulge=1.
        let stock = Contour2D::from_vertices(
            &[
                (5.0, 0.0, 1.0),  // bulge=1 → semicircle CCW
                (-5.0, 0.0, 1.0), // bulge=1 → second semicircle
            ],
            true,
        );
        // Tool centered on the boundary of the circle.
        let result = engagement_angle(Point2::new(5.0, 0.0), 1.0, &stock);
        // Should be roughly π (half engaged) but slightly less because the
        // stock boundary curves away.
        assert!(result.angle > 0.5, "should have some engagement");
        assert!(result.angle < TAU, "should not be full slotting");
    }

    #[test]
    fn engagement_from_depth_zero() {
        let eng = engagement_from_radial_depth(5.0, 0.0);
        assert_relative_eq!(eng, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn simple_convenience_matches_full() {
        let stock = square_stock();
        let center = Point2::new(4.0, 0.0);
        let full = engagement_angle(center, 2.0, &stock);
        let simple = engagement_angle_simple(center, 2.0, &stock);
        assert_relative_eq!(full.angle, simple, epsilon = 1e-12);
    }
}

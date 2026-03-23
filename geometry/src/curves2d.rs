// 2D curve types and curve-curve intersection
// Used by toolpath crate for cutter comp, pocketing, contour clipping

use std::f64::consts::TAU;

use serde::{Deserialize, Serialize};

use crate::types::*;

// ---------------------------------------------------------------------------
// 2D curve types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineSegment2D {
    pub start: Point2,
    pub end: Point2,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Arc2D {
    pub center: Point2,
    pub radius: f64,
    pub start_angle: f64, // radians
    pub end_angle: f64,   // radians, arc goes CCW from start to end
}

// ---------------------------------------------------------------------------
// Angle helpers
// ---------------------------------------------------------------------------

/// Normalize angle to [0, 2pi).
fn normalize_angle(a: f64) -> f64 {
    let mut r = a % TAU;
    if r < 0.0 {
        r += TAU;
    }
    r
}

/// Angular span from start to end going CCW. Returns value in (0, 2pi].
fn arc_span(start: f64, end: f64) -> f64 {
    let mut span = (end - start) % TAU;
    if span <= 0.0 {
        span += TAU;
    }
    span
}

// ---------------------------------------------------------------------------
// Arc2D methods
// ---------------------------------------------------------------------------

impl Arc2D {
    /// Angular span in radians (always positive).
    pub fn span(&self) -> f64 {
        arc_span(self.start_angle, self.end_angle)
    }

    /// Check if an angle (in radians) falls within the arc's angular range.
    pub fn contains_angle(&self, angle: f64) -> bool {
        let a = normalize_angle(angle - self.start_angle);
        a <= self.span() + 1e-10
    }

    /// Evaluate position at parameter t in [0, 1].
    pub fn eval(&self, t: f64) -> Point2 {
        let theta = self.start_angle + t * self.span();
        Point2::new(
            self.center.x + self.radius * theta.cos(),
            self.center.y + self.radius * theta.sin(),
        )
    }

    /// Unit tangent at parameter t in [0, 1] (CCW direction).
    pub fn tangent(&self, t: f64) -> Vec2 {
        let theta = self.start_angle + t * self.span();
        Vec2::new(-theta.sin(), theta.cos())
    }

    /// Arc length.
    pub fn length(&self) -> f64 {
        self.radius.abs() * self.span()
    }
}

impl LineSegment2D {
    /// Evaluate at t in [0, 1].
    pub fn eval(&self, t: f64) -> Point2 {
        let d = self.end - self.start;
        self.start + d * t
    }

    /// Unit tangent.
    pub fn tangent(&self) -> Vec2 {
        (self.end - self.start).normalize()
    }

    /// Length.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }
}

// ---------------------------------------------------------------------------
// 2D line-line intersection
// ---------------------------------------------------------------------------

/// Intersect two 2D line segments. Returns the intersection point if the
/// segments cross (both parameters in [0, 1]).
pub fn intersect_lines_2d(a: &LineSegment2D, b: &LineSegment2D) -> Option<Point2> {
    let da = a.end - a.start;
    let db = b.end - b.start;

    let det = da.x * db.y - da.y * db.x;
    if det.abs() < 1e-15 {
        return None; // parallel or coincident
    }

    let d = b.start - a.start;
    let t = (d.x * db.y - d.y * db.x) / det;
    let s = (d.x * da.y - d.y * da.x) / det;

    if t >= -1e-12 && t <= 1.0 + 1e-12 && s >= -1e-12 && s <= 1.0 + 1e-12 {
        Some(a.start + da * t)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// 2D line-arc intersection
// ---------------------------------------------------------------------------

/// Intersect a 2D line segment with a 2D circular arc.
/// Returns all intersection points (0, 1, or 2).
pub fn intersect_line_arc_2d(line: &LineSegment2D, arc: &Arc2D) -> Vec<Point2> {
    let d = line.end - line.start;
    let f = line.start - arc.center;

    let a = d.dot(&d);
    let b = 2.0 * f.dot(&d);
    let c = f.dot(&f) - arc.radius * arc.radius;

    if a.abs() < 1e-30 {
        return vec![]; // degenerate line
    }

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return vec![];
    }

    let sqrt_disc = disc.sqrt();
    let mut results = Vec::new();

    for t in [(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)] {
        if t >= -1e-12 && t <= 1.0 + 1e-12 {
            let p = line.start + d * t.clamp(0.0, 1.0);
            let angle = (p.y - arc.center.y).atan2(p.x - arc.center.x);
            if arc.contains_angle(angle) {
                results.push(p);
            }
        }
    }

    // Deduplicate near-coincident hits (tangent case)
    if results.len() == 2 {
        if (results[0] - results[1]).norm() < 1e-10 {
            results.pop();
        }
    }

    results
}

// ---------------------------------------------------------------------------
// 2D arc-arc intersection
// ---------------------------------------------------------------------------

/// Intersect two 2D circular arcs.
/// Returns all intersection points (0, 1, or 2).
pub fn intersect_arc_arc_2d(a: &Arc2D, b: &Arc2D) -> Vec<Point2> {
    let d = b.center - a.center;
    let dist = d.norm();

    if dist < 1e-15 {
        return vec![]; // concentric
    }
    if dist > a.radius + b.radius + 1e-10 {
        return vec![]; // too far apart
    }
    if dist < (a.radius - b.radius).abs() - 1e-10 {
        return vec![]; // one fully inside the other
    }

    // Distance from a.center to the radical line along the center-to-center axis
    let a_dist =
        (a.radius * a.radius - b.radius * b.radius + dist * dist) / (2.0 * dist);
    let h_sq = a.radius * a.radius - a_dist * a_dist;
    let h = if h_sq > 0.0 { h_sq.sqrt() } else { 0.0 };

    let unit_d = d / dist;
    let perp = Vec2::new(-unit_d.y, unit_d.x);
    let mid = a.center + unit_d * a_dist;

    let mut results = Vec::new();

    if h < 1e-12 {
        // Tangent — single point
        let p = mid;
        let angle_a = (p.y - a.center.y).atan2(p.x - a.center.x);
        let angle_b = (p.y - b.center.y).atan2(p.x - b.center.x);
        if a.contains_angle(angle_a) && b.contains_angle(angle_b) {
            results.push(p);
        }
    } else {
        for sign in [-1.0_f64, 1.0] {
            let p = mid + perp * (sign * h);
            let angle_a = (p.y - a.center.y).atan2(p.x - a.center.x);
            let angle_b = (p.y - b.center.y).atan2(p.x - b.center.x);
            if a.contains_angle(angle_a) && b.contains_angle(angle_b) {
                results.push(p);
            }
        }
    }

    results
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    // -- Arc2D --

    #[test]
    fn test_arc_span_normal() {
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI,
        };
        assert_relative_eq!(arc.span(), PI, epsilon = 1e-12);
    }

    #[test]
    fn test_arc_span_wrapping() {
        // Arc from 350 to 10 degrees (20 degree span crossing 0)
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 350.0_f64.to_radians(),
            end_angle: 10.0_f64.to_radians(),
        };
        assert_relative_eq!(arc.span(), 20.0_f64.to_radians(), epsilon = 1e-10);
    }

    #[test]
    fn test_arc_contains_angle() {
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI,
        };
        assert!(arc.contains_angle(FRAC_PI_2));  // 90 deg — inside
        assert!(!arc.contains_angle(PI + 0.5));    // 180+28 deg — outside
    }

    #[test]
    fn test_arc_contains_angle_wrapping() {
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 5.5,
            end_angle: 0.5,
        };
        // Wraps through 0. Should contain 0.0 and 6.0 but not 3.0
        assert!(arc.contains_angle(0.0));
        assert!(arc.contains_angle(6.0));
        assert!(!arc.contains_angle(3.0));
    }

    #[test]
    fn test_arc_eval_endpoints() {
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 2.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let p0 = arc.eval(0.0);
        assert_relative_eq!(p0.x, 2.0, epsilon = 1e-12);
        assert_relative_eq!(p0.y, 0.0, epsilon = 1e-12);
        let p1 = arc.eval(1.0);
        assert_relative_eq!(p1.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(p1.y, 2.0, epsilon = 1e-12);
    }

    // -- Line-Line --

    #[test]
    fn test_lines_2d_cross() {
        let a = LineSegment2D {
            start: Point2::new(0.0, 0.0),
            end: Point2::new(2.0, 0.0),
        };
        let b = LineSegment2D {
            start: Point2::new(1.0, -1.0),
            end: Point2::new(1.0, 1.0),
        };
        let p = intersect_lines_2d(&a, &b).unwrap();
        assert_relative_eq!(p.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(p.y, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_lines_2d_parallel() {
        let a = LineSegment2D {
            start: Point2::new(0.0, 0.0),
            end: Point2::new(1.0, 0.0),
        };
        let b = LineSegment2D {
            start: Point2::new(0.0, 1.0),
            end: Point2::new(1.0, 1.0),
        };
        assert!(intersect_lines_2d(&a, &b).is_none());
    }

    #[test]
    fn test_lines_2d_no_overlap() {
        // Lines would intersect if extended, but segments don't overlap
        let a = LineSegment2D {
            start: Point2::new(0.0, 0.0),
            end: Point2::new(1.0, 0.0),
        };
        let b = LineSegment2D {
            start: Point2::new(2.0, -1.0),
            end: Point2::new(2.0, 1.0),
        };
        assert!(intersect_lines_2d(&a, &b).is_none());
    }

    // -- Line-Arc --

    #[test]
    fn test_line_arc_two_hits() {
        let line = LineSegment2D {
            start: Point2::new(-2.0, 0.0),
            end: Point2::new(2.0, 0.0),
        };
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: TAU - 0.01, // nearly full circle
        };
        let hits = intersect_line_arc_2d(&line, &arc);
        assert_eq!(hits.len(), 2);
        // Should hit at (-1, 0) and (1, 0)
        let xs: Vec<f64> = hits.iter().map(|p| p.x).collect();
        assert!(xs.iter().any(|x| (x - 1.0).abs() < 1e-10));
        assert!(xs.iter().any(|x| (x + 1.0).abs() < 1e-10));
    }

    #[test]
    fn test_line_arc_miss() {
        let line = LineSegment2D {
            start: Point2::new(-2.0, 5.0),
            end: Point2::new(2.0, 5.0),
        };
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI,
        };
        assert!(intersect_line_arc_2d(&line, &arc).is_empty());
    }

    #[test]
    fn test_line_arc_one_hit_segment_limit() {
        // Line hits the full circle at two points, but the arc only covers the top half
        let line = LineSegment2D {
            start: Point2::new(-2.0, 0.5),
            end: Point2::new(2.0, 0.5),
        };
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI, // top half only
        };
        let hits = intersect_line_arc_2d(&line, &arc);
        // Both circle hits are in the upper half, so both should be on the arc
        // hit at y=0.5: x = ±sqrt(1 - 0.25) = ±0.866
        // angle for (0.866, 0.5) ≈ 30° — in [0, 180°] ✓
        // angle for (-0.866, 0.5) ≈ 150° — in [0, 180°] ✓
        assert_eq!(hits.len(), 2);
    }

    // -- Arc-Arc --

    #[test]
    fn test_arc_arc_two_hits() {
        // Two overlapping semicircles
        let a = Arc2D {
            center: Point2::new(0.0, 0.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: TAU - 0.01,
        };
        let b = Arc2D {
            center: Point2::new(1.0, 0.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: TAU - 0.01,
        };
        let hits = intersect_arc_arc_2d(&a, &b);
        assert_eq!(hits.len(), 2);
        // Intersection at (0.5, ±sqrt(3)/2)
        for p in &hits {
            assert_relative_eq!(p.x, 0.5, epsilon = 1e-10);
            assert_relative_eq!(p.y.abs(), (3.0_f64).sqrt() / 2.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_arc_arc_no_overlap() {
        let a = Arc2D {
            center: Point2::new(0.0, 0.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI,
        };
        let b = Arc2D {
            center: Point2::new(5.0, 0.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI,
        };
        assert!(intersect_arc_arc_2d(&a, &b).is_empty());
    }

    #[test]
    fn test_arc_arc_tangent() {
        let a = Arc2D {
            center: Point2::new(0.0, 0.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: TAU - 0.01,
        };
        let b = Arc2D {
            center: Point2::new(2.0, 0.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: TAU - 0.01,
        };
        let hits = intersect_arc_arc_2d(&a, &b);
        assert_eq!(hits.len(), 1);
        assert_relative_eq!(hits[0].x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(hits[0].y, 0.0, epsilon = 1e-10);
    }
}

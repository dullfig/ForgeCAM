// 2D contour (chain of line/arc segments) — wraps cavalier_contours
//
// A Contour2D is a connected sequence of 2D line and arc segments, closed or open.
// This is the fundamental representation for toolpath boundaries, pocket profiles,
// cutter comp paths, and stock outlines.

use cavalier_contours::polyline::{
    PlineCreation, PlineOffsetOptions, PlineSource, PlineSourceMut, Polyline,
};
use serde::{Deserialize, Serialize};

use crate::curves2d::{Arc2D, LineSegment2D};
use crate::types::*;

// ---------------------------------------------------------------------------
// Segment2D — a single element in a contour chain
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Segment2D {
    Line(LineSegment2D),
    Arc(Arc2D),
}

// ---------------------------------------------------------------------------
// BooleanOp
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    Union,
    Intersection,
    Difference,
    Xor,
}

// ---------------------------------------------------------------------------
// Contour2D
// ---------------------------------------------------------------------------

/// A connected chain of 2D line and arc segments.
///
/// Wraps `cavalier_contours::Polyline<f64>` for offset and boolean operations
/// while providing conversion to/from ForgeCAM's curve types.
#[derive(Debug, Clone)]
pub struct Contour2D {
    inner: Polyline<f64>,
}

impl Contour2D {
    // -------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------

    /// Create from cavalier Polyline directly.
    pub fn from_polyline(pline: Polyline<f64>) -> Self {
        Contour2D { inner: pline }
    }

    /// Create from raw (x, y, bulge) tuples.
    /// bulge = 0 for lines, tan(sweep/4) for arcs (positive = CCW).
    pub fn from_vertices(vertices: &[(f64, f64, f64)], closed: bool) -> Self {
        let mut pline = Polyline::with_capacity(vertices.len(), closed);
        for &(x, y, bulge) in vertices {
            pline.add(x, y, bulge);
        }
        Contour2D { inner: pline }
    }

    /// Create from a sequence of connected Segment2D.
    pub fn from_segments(segments: &[Segment2D], closed: bool) -> Self {
        if segments.is_empty() {
            return Contour2D {
                inner: Polyline::new_closed(),
            };
        }

        let mut pline = Polyline::with_capacity(segments.len() + 1, closed);

        for seg in segments {
            match seg {
                Segment2D::Line(line) => {
                    pline.add(line.start.x, line.start.y, 0.0);
                }
                Segment2D::Arc(arc) => {
                    let start_pt = arc.eval(0.0);
                    let bulge = bulge_from_arc(arc);
                    pline.add(start_pt.x, start_pt.y, bulge);
                }
            }
        }

        // Add final vertex (end of last segment) if open
        if !closed {
            let last = segments.last().unwrap();
            let end_pt = match last {
                Segment2D::Line(l) => l.end,
                Segment2D::Arc(a) => a.eval(1.0),
            };
            pline.add(end_pt.x, end_pt.y, 0.0);
        }

        Contour2D { inner: pline }
    }

    /// Create a closed contour from a set of points (all line segments).
    pub fn from_points(points: &[Point2], closed: bool) -> Self {
        let mut pline = Polyline::with_capacity(points.len(), closed);
        for p in points {
            pline.add(p.x, p.y, 0.0);
        }
        Contour2D { inner: pline }
    }

    /// Whether this contour is closed (forms a loop).
    pub fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.inner.vertex_count()
    }

    // -------------------------------------------------------------------
    // Conversion back to segments
    // -------------------------------------------------------------------

    /// Convert the contour back to a sequence of Segment2D.
    pub fn to_segments(&self) -> Vec<Segment2D> {
        let n = self.inner.vertex_count();
        if n < 2 {
            return vec![];
        }

        let seg_count = if self.inner.is_closed() { n } else { n - 1 };
        let mut segments = Vec::with_capacity(seg_count);

        for i in 0..seg_count {
            let j = (i + 1) % n;
            let v0 = self.inner.at(i);
            let v1 = self.inner.at(j);

            let start = Point2::new(v0.x, v0.y);
            let end = Point2::new(v1.x, v1.y);

            if v0.bulge.abs() < 1e-15 {
                segments.push(Segment2D::Line(LineSegment2D { start, end }));
            } else {
                segments.push(Segment2D::Arc(arc_from_bulge(start, end, v0.bulge)));
            }
        }

        segments
    }

    /// Get the raw vertices as (x, y, bulge) tuples.
    pub fn vertices(&self) -> Vec<(f64, f64, f64)> {
        (0..self.inner.vertex_count())
            .map(|i| {
                let v = self.inner.at(i);
                (v.x, v.y, v.bulge)
            })
            .collect()
    }

    // -------------------------------------------------------------------
    // Geometric queries
    // -------------------------------------------------------------------

    /// Signed area (positive = CCW winding, negative = CW). Only meaningful for closed contours.
    pub fn area(&self) -> f64 {
        self.inner.area()
    }

    /// Total path length.
    pub fn path_length(&self) -> f64 {
        self.inner.path_length()
    }

    /// Winding number of a point relative to this closed contour.
    pub fn winding_number(&self, point: Point2) -> i32 {
        use cavalier_contours::core::math::Vector2;
        self.inner.winding_number(Vector2::new(point.x, point.y))
    }

    /// Axis-aligned bounding box.
    pub fn extents(&self) -> Option<(Point2, Point2)> {
        self.inner.extents().map(|aabb| {
            (
                Point2::new(aabb.min_x, aabb.min_y),
                Point2::new(aabb.max_x, aabb.max_y),
            )
        })
    }

    // -------------------------------------------------------------------
    // Offset
    // -------------------------------------------------------------------

    /// Parallel offset by `distance`.
    /// Positive = offset to the left of travel direction (outward for CCW contours).
    /// Returns multiple contours (offset may split or create islands).
    pub fn offset(&self, distance: f64) -> Vec<Contour2D> {
        self.inner
            .parallel_offset(distance)
            .into_iter()
            .map(Contour2D::from_polyline)
            .collect()
    }

    /// Parallel offset with custom options.
    pub fn offset_with_options(
        &self,
        distance: f64,
        options: &PlineOffsetOptions<f64>,
    ) -> Vec<Contour2D> {
        self.inner
            .parallel_offset_opt(distance, options)
            .into_iter()
            .map(Contour2D::from_polyline)
            .collect()
    }

    // -------------------------------------------------------------------
    // Boolean operations
    // -------------------------------------------------------------------

    /// Boolean operation between two closed contours.
    /// Returns the resulting contour(s).
    pub fn boolean(&self, other: &Contour2D, op: BooleanOp) -> Vec<Contour2D> {
        use cavalier_contours::polyline::BooleanOp as CavcBoolOp;

        let cavc_op = match op {
            BooleanOp::Union => CavcBoolOp::Or,
            BooleanOp::Intersection => CavcBoolOp::And,
            BooleanOp::Difference => CavcBoolOp::Not,
            BooleanOp::Xor => CavcBoolOp::Xor,
        };

        let result = self.inner.boolean(&other.inner, cavc_op);
        result
            .pos_plines
            .into_iter()
            .chain(result.neg_plines.into_iter())
            .map(|r| Contour2D::from_polyline(r.pline))
            .collect()
    }

    // -------------------------------------------------------------------
    // Access to inner polyline (escape hatch)
    // -------------------------------------------------------------------

    /// Borrow the underlying cavalier Polyline.
    pub fn as_polyline(&self) -> &Polyline<f64> {
        &self.inner
    }

    /// Consume and return the underlying cavalier Polyline.
    pub fn into_polyline(self) -> Polyline<f64> {
        self.inner
    }
}

// ---------------------------------------------------------------------------
// Bulge conversion helpers
// ---------------------------------------------------------------------------

/// Compute bulge value from an Arc2D.
/// bulge = tan(sweep/4), positive for CCW arcs.
fn bulge_from_arc(arc: &Arc2D) -> f64 {
    let sweep = arc.span(); // always positive (CCW)
    (sweep / 4.0).tan()
}

/// Reconstruct an Arc2D from two endpoints and a bulge value.
fn arc_from_bulge(start: Point2, end: Point2, bulge: f64) -> Arc2D {
    let dx = end.x - start.x;
    let dy = end.y - start.y;
    let chord_len = (dx * dx + dy * dy).sqrt();

    if chord_len < 1e-15 {
        // Degenerate: start == end
        return Arc2D {
            center: start,
            radius: 0.0,
            start_angle: 0.0,
            end_angle: 0.0,
        };
    }

    // Sweep angle (signed: positive = CCW, negative = CW)
    let sweep = 4.0 * bulge.atan();

    // Radius from chord length and sweep
    let half_sweep = sweep / 2.0;
    let radius = (chord_len / 2.0) / half_sweep.sin().abs();

    // Center: offset from chord midpoint perpendicular to chord
    let mx = (start.x + end.x) / 2.0;
    let my = (start.y + end.y) / 2.0;

    // Left perpendicular direction (normalized)
    let px = -dy / chord_len;
    let py = dx / chord_len;

    // Signed distance from midpoint to center along perpendicular
    let d = (chord_len / 2.0) / half_sweep.tan();

    let center = Point2::new(mx + d * px, my + d * py);

    // Compute angles
    let start_angle = (start.y - center.y).atan2(start.x - center.x);
    let end_angle = start_angle + sweep;

    Arc2D {
        center,
        radius,
        start_angle,
        end_angle,
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    fn unit_square() -> Contour2D {
        Contour2D::from_vertices(
            &[
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0),
            ],
            true,
        )
    }

    #[test]
    fn test_from_vertices_and_back() {
        let contour = unit_square();
        assert_eq!(contour.vertex_count(), 4);
        assert!(contour.is_closed());
        let segs = contour.to_segments();
        assert_eq!(segs.len(), 4); // 4 segments for closed square
    }

    #[test]
    fn test_area_unit_square() {
        let sq = unit_square();
        let area = sq.area();
        assert_relative_eq!(area.abs(), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_path_length_unit_square() {
        let sq = unit_square();
        assert_relative_eq!(sq.path_length(), 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_offset_square_outward() {
        let sq = unit_square();
        // Negative offset = outward for CCW contours in cavalier convention
        let offsets = sq.offset(-0.1);
        assert!(!offsets.is_empty());
        let outer_area = offsets[0].area().abs();
        assert!(outer_area > 1.0);
    }

    #[test]
    fn test_offset_square_inward() {
        let sq = unit_square();
        // Positive offset = inward for CCW contours
        let offsets = sq.offset(0.1);
        assert!(!offsets.is_empty());
        let inner_area = offsets[0].area().abs();
        assert!(inner_area < 1.0);
    }

    #[test]
    fn test_bulge_roundtrip() {
        let arc = Arc2D {
            center: Point2::origin(),
            radius: 5.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let bulge = bulge_from_arc(&arc);
        assert!(bulge > 0.0); // CCW

        let start = arc.eval(0.0);
        let end = arc.eval(1.0);
        let reconstructed = arc_from_bulge(start, end, bulge);

        assert_relative_eq!(reconstructed.radius, 5.0, epsilon = 1e-8);
        assert_relative_eq!(reconstructed.center.x, 0.0, epsilon = 1e-8);
        assert_relative_eq!(reconstructed.center.y, 0.0, epsilon = 1e-8);
    }

    #[test]
    fn test_from_segments() {
        let segments = vec![
            Segment2D::Line(LineSegment2D {
                start: Point2::new(0.0, 0.0),
                end: Point2::new(1.0, 0.0),
            }),
            Segment2D::Arc(Arc2D {
                center: Point2::new(1.0, 0.5),
                radius: 0.5,
                start_angle: -FRAC_PI_2,
                end_angle: FRAC_PI_2,
            }),
            Segment2D::Line(LineSegment2D {
                start: Point2::new(1.0, 1.0),
                end: Point2::new(0.0, 1.0),
            }),
        ];
        let contour = Contour2D::from_segments(&segments, false);
        assert_eq!(contour.vertex_count(), 4); // 3 segments + final endpoint
    }

    #[test]
    fn test_winding_number() {
        let sq = unit_square();
        // Point inside
        let wn_inside = sq.winding_number(Point2::new(0.5, 0.5));
        assert_ne!(wn_inside, 0);
        // Point outside
        let wn_outside = sq.winding_number(Point2::new(5.0, 5.0));
        assert_eq!(wn_outside, 0);
    }

    #[test]
    fn test_extents() {
        let sq = unit_square();
        let (min, max) = sq.extents().unwrap();
        assert_relative_eq!(min.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(min.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(max.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(max.y, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_boolean_union() {
        // Two overlapping squares
        let a = unit_square();
        let b = Contour2D::from_vertices(
            &[
                (0.5, 0.5, 0.0),
                (1.5, 0.5, 0.0),
                (1.5, 1.5, 0.0),
                (0.5, 1.5, 0.0),
            ],
            true,
        );
        let result = a.boolean(&b, BooleanOp::Union);
        assert!(!result.is_empty());
        // Union area should be > 1.0 (two overlapping unit squares)
        let total_area: f64 = result.iter().map(|c| c.area().abs()).sum();
        assert!(total_area > 1.0);
    }
}

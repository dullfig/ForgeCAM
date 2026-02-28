/// 2D polygon utilities for face splitting and point-in-polygon tests.

/// A 2D point.
pub type Point2 = [f64; 2];

/// A hit from intersecting a line with a polygon edge.
#[derive(Debug, Clone)]
pub struct LineHit {
    /// Parameter along the line: hit_point = line_origin + t_line * line_dir.
    pub t_line: f64,
    /// Which polygon edge was hit (index into vertices, edge i → i+1).
    pub edge_index: usize,
    /// Parameter along the polygon edge [0, 1].
    pub t_edge: f64,
}

/// Classification of a point relative to a polygon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointClassification {
    Inside,
    Outside,
    OnBoundary,
}

/// A 2D polygon as an ordered list of vertices (CCW winding = positive area).
#[derive(Debug, Clone)]
pub struct Polygon2D {
    pub vertices: Vec<Point2>,
}

impl Polygon2D {
    /// Signed area via the shoelace formula. Positive = CCW winding.
    pub fn signed_area(&self) -> f64 {
        let n = self.vertices.len();
        if n < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += self.vertices[i][0] * self.vertices[j][1];
            area -= self.vertices[j][0] * self.vertices[i][1];
        }
        area * 0.5
    }

    /// Classify a point as Inside, Outside, or OnBoundary using the crossing number test.
    pub fn classify_point(&self, p: Point2, tol: f64) -> PointClassification {
        let n = self.vertices.len();
        if n < 3 {
            return PointClassification::Outside;
        }

        // First check if the point is on any edge.
        for i in 0..n {
            let j = (i + 1) % n;
            if point_on_segment(p, self.vertices[i], self.vertices[j], tol) {
                return PointClassification::OnBoundary;
            }
        }

        // Crossing number (ray cast in +x direction).
        let mut crossings = 0i32;
        for i in 0..n {
            let j = (i + 1) % n;
            let yi = self.vertices[i][1];
            let yj = self.vertices[j][1];

            // Does this edge straddle the ray's y-level?
            if (yi <= p[1] && yj > p[1]) || (yj <= p[1] && yi > p[1]) {
                // Compute x-intercept of the edge with y = p[1].
                let t = (p[1] - yi) / (yj - yi);
                let x_intercept = self.vertices[i][0] + t * (self.vertices[j][0] - self.vertices[i][0]);
                if p[0] < x_intercept {
                    crossings += 1;
                }
            }
        }

        if crossings % 2 == 1 {
            PointClassification::Inside
        } else {
            PointClassification::Outside
        }
    }

    /// Intersect an infinite line with the polygon boundary.
    /// Line is defined as: point = line_origin + t * line_dir.
    /// Returns hits sorted by t_line.
    pub fn intersect_line(&self, line_origin: Point2, line_dir: Point2) -> Vec<LineHit> {
        let n = self.vertices.len();
        let mut hits = Vec::new();

        for i in 0..n {
            let j = (i + 1) % n;
            let edge_start = self.vertices[i];
            let edge_end = self.vertices[j];
            let edge_dir = [edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]];

            // Solve: line_origin + t_line * line_dir = edge_start + t_edge * edge_dir
            // [ line_dir.x  -edge_dir.x ] [ t_line ] = [ edge_start.x - line_origin.x ]
            // [ line_dir.y  -edge_dir.y ] [ t_edge ]   [ edge_start.y - line_origin.y ]
            let det = line_dir[0] * (-edge_dir[1]) - line_dir[1] * (-edge_dir[0]);

            if det.abs() < 1e-12 {
                // Parallel — skip.
                continue;
            }

            let dx = edge_start[0] - line_origin[0];
            let dy = edge_start[1] - line_origin[1];

            let t_line = (dx * (-edge_dir[1]) - dy * (-edge_dir[0])) / det;
            let t_edge = (line_dir[0] * dy - line_dir[1] * dx) / det;

            // Accept if hit is within the edge (with small tolerance for endpoints).
            if t_edge >= -1e-10 && t_edge <= 1.0 + 1e-10 {
                let t_edge_clamped = t_edge.clamp(0.0, 1.0);
                hits.push(LineHit {
                    t_line,
                    edge_index: i,
                    t_edge: t_edge_clamped,
                });
            }
        }

        // Sort by t_line.
        hits.sort_by(|a, b| a.t_line.partial_cmp(&b.t_line).unwrap());

        // Deduplicate hits at nearly the same t_line (e.g., line passing through a vertex
        // hits two adjacent edges).
        dedup_hits(&mut hits);

        hits
    }
}

/// Check if point p lies on the segment from a to b within tolerance.
fn point_on_segment(p: Point2, a: Point2, b: Point2, tol: f64) -> bool {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];

    // Cross product (perpendicular distance * |ab|).
    let cross = ab[0] * ap[1] - ab[1] * ap[0];
    let len_sq = ab[0] * ab[0] + ab[1] * ab[1];

    if len_sq < tol * tol {
        // Degenerate edge — check distance to point a.
        let dist_sq = ap[0] * ap[0] + ap[1] * ap[1];
        return dist_sq < tol * tol;
    }

    // Check perpendicular distance.
    if cross * cross / len_sq > tol * tol {
        return false;
    }

    // Check projection parameter.
    let dot = ab[0] * ap[0] + ab[1] * ap[1];
    let t = dot / len_sq;
    t >= -tol && t <= 1.0 + tol
}

/// Remove duplicate hits (same vertex hit from two adjacent edges).
fn dedup_hits(hits: &mut Vec<LineHit>) {
    if hits.len() <= 1 {
        return;
    }
    let mut i = 0;
    while i + 1 < hits.len() {
        if (hits[i].t_line - hits[i + 1].t_line).abs() < 1e-9 {
            hits.remove(i + 1);
        } else {
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_square() -> Polygon2D {
        Polygon2D {
            vertices: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        }
    }

    #[test]
    fn test_signed_area_unit_square() {
        let sq = unit_square();
        let area = sq.signed_area();
        assert!((area - 1.0).abs() < 1e-12, "Area should be 1.0, got {area}");
    }

    #[test]
    fn test_signed_area_cw_is_negative() {
        let sq = Polygon2D {
            vertices: vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]],
        };
        assert!(sq.signed_area() < 0.0);
    }

    #[test]
    fn test_point_inside_square() {
        let sq = unit_square();
        assert_eq!(sq.classify_point([0.5, 0.5], 1e-9), PointClassification::Inside);
    }

    #[test]
    fn test_point_outside_square() {
        let sq = unit_square();
        assert_eq!(sq.classify_point([2.0, 0.5], 1e-9), PointClassification::Outside);
        assert_eq!(sq.classify_point([-1.0, 0.5], 1e-9), PointClassification::Outside);
    }

    #[test]
    fn test_point_on_boundary() {
        let sq = unit_square();
        assert_eq!(sq.classify_point([0.5, 0.0], 1e-9), PointClassification::OnBoundary);
        assert_eq!(sq.classify_point([1.0, 0.5], 1e-9), PointClassification::OnBoundary);
        assert_eq!(sq.classify_point([0.0, 0.0], 1e-9), PointClassification::OnBoundary);
    }

    #[test]
    fn test_line_intersects_square_horizontal() {
        let sq = unit_square();
        // Horizontal line through center: y=0.5, direction = (1, 0)
        let hits = sq.intersect_line([0.0, 0.5], [1.0, 0.0]);
        assert_eq!(hits.len(), 2, "Expected 2 hits, got {}", hits.len());
        // Should hit left edge (x=0) and right edge (x=1)
        assert!((hits[0].t_line - 0.0).abs() < 1e-9, "First hit t={}", hits[0].t_line);
        assert!((hits[1].t_line - 1.0).abs() < 1e-9, "Second hit t={}", hits[1].t_line);
    }

    #[test]
    fn test_line_misses_square() {
        let sq = unit_square();
        // Horizontal line above the square
        let hits = sq.intersect_line([0.0, 2.0], [1.0, 0.0]);
        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn test_line_through_vertex() {
        let sq = unit_square();
        // Diagonal line through origin corner (0,0) to (1,1)
        let hits = sq.intersect_line([0.0, 0.0], [1.0, 1.0]);
        // Should get exactly 2 hits (at corners (0,0) and (1,1)), after dedup
        assert_eq!(hits.len(), 2, "Expected 2 hits through corners, got {}", hits.len());
        assert!((hits[0].t_line - 0.0).abs() < 1e-9);
        assert!((hits[1].t_line - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_line_vertical_through_square() {
        let sq = unit_square();
        // Vertical line at x=0.5, direction = (0, 1)
        let hits = sq.intersect_line([0.5, 0.0], [0.0, 1.0]);
        assert_eq!(hits.len(), 2);
        assert!((hits[0].t_line - 0.0).abs() < 1e-9, "Bottom hit t={}", hits[0].t_line);
        assert!((hits[1].t_line - 1.0).abs() < 1e-9, "Top hit t={}", hits[1].t_line);
    }
}

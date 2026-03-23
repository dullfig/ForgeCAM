// 2D polygon operations: point classification, line intersection, signed area

use serde::{Deserialize, Serialize};

/// Lightweight 2D point — no nalgebra overhead.
pub type Point2D = [f64; 2];

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineHit {
    /// Parameter along the infinite line.
    pub t_line: f64,
    /// Which polygon edge was hit (i → i+1).
    pub edge_index: usize,
    /// Parameter along edge [0, 1].
    pub t_edge: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PointClassification {
    Inside,
    Outside,
    OnBoundary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Polygon2D {
    pub vertices: Vec<Point2D>,
}

impl Polygon2D {
    /// Shoelace formula. Positive = CCW winding, negative = CW.
    pub fn signed_area(&self) -> f64 {
        let n = self.vertices.len();
        if n < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += self.vertices[i][0] * self.vertices[j][1]
                - self.vertices[j][0] * self.vertices[i][1];
        }
        area * 0.5
    }

    /// Point-in-polygon classification via crossing number + boundary snap.
    pub fn classify_point(&self, p: Point2D, tol: f64) -> PointClassification {
        let n = self.vertices.len();
        if n < 3 {
            return PointClassification::Outside;
        }

        // Check boundary proximity first
        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = (self.vertices[i][0], self.vertices[i][1]);
            let (bx, by) = (self.vertices[j][0], self.vertices[j][1]);

            // Distance from point to edge segment
            let dx = bx - ax;
            let dy = by - ay;
            let len_sq = dx * dx + dy * dy;

            if len_sq < tol * tol {
                // Degenerate edge — check distance to vertex
                let dist_sq = (p[0] - ax) * (p[0] - ax) + (p[1] - ay) * (p[1] - ay);
                if dist_sq <= tol * tol {
                    return PointClassification::OnBoundary;
                }
                continue;
            }

            let t = ((p[0] - ax) * dx + (p[1] - ay) * dy) / len_sq;
            let t_clamped = t.clamp(0.0, 1.0);
            let closest_x = ax + t_clamped * dx;
            let closest_y = ay + t_clamped * dy;
            let dist_sq =
                (p[0] - closest_x) * (p[0] - closest_x) + (p[1] - closest_y) * (p[1] - closest_y);

            if dist_sq <= tol * tol {
                return PointClassification::OnBoundary;
            }
        }

        // Ray casting in +x direction
        let mut crossings = 0;
        for i in 0..n {
            let j = (i + 1) % n;
            let (yi, yj) = (self.vertices[i][1], self.vertices[j][1]);

            // Edge must straddle the ray's y-coordinate
            if (yi <= p[1] && yj > p[1]) || (yj <= p[1] && yi > p[1]) {
                // Compute x-intersection of ray with edge
                let t = (p[1] - yi) / (yj - yi);
                let x_intersect = self.vertices[i][0] + t * (self.vertices[j][0] - self.vertices[i][0]);
                if p[0] < x_intersect {
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
    /// Returns hits sorted by t_line parameter.
    pub fn intersect_line(&self, line_origin: Point2D, line_dir: Point2D) -> Vec<LineHit> {
        let n = self.vertices.len();
        let mut hits = Vec::new();

        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = (self.vertices[i][0], self.vertices[i][1]);
            let (bx, by) = (self.vertices[j][0], self.vertices[j][1]);

            let edge_dx = bx - ax;
            let edge_dy = by - ay;

            // Solve: line_origin + t * line_dir = edge_start + s * edge_dir
            // [line_dir.x  -edge_dx] [t]   [ax - line_origin.x]
            // [line_dir.y  -edge_dy] [s] = [ay - line_origin.y]
            let det = line_dir[0] * (-edge_dy) - line_dir[1] * (-edge_dx);

            if det.abs() < 1e-15 {
                continue; // parallel
            }

            let rhs_x = ax - line_origin[0];
            let rhs_y = ay - line_origin[1];

            let t = (rhs_x * (-edge_dy) - rhs_y * (-edge_dx)) / det;
            let s = (line_dir[0] * rhs_y - line_dir[1] * rhs_x) / det;

            if s >= 0.0 && s <= 1.0 {
                hits.push(LineHit {
                    t_line: t,
                    edge_index: i,
                    t_edge: s,
                });
            }
        }

        // Sort by t_line
        hits.sort_by(|a, b| a.t_line.partial_cmp(&b.t_line).unwrap());

        // Deduplicate consecutive hits within 1e-9 (shared vertices)
        if hits.len() > 1 {
            let mut deduped = vec![hits[0].clone()];
            for h in &hits[1..] {
                if (h.t_line - deduped.last().unwrap().t_line).abs() > 1e-9 {
                    deduped.push(h.clone());
                }
            }
            hits = deduped;
        }

        hits
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
    fn test_signed_area_ccw_square() {
        let sq = unit_square();
        let area = sq.signed_area();
        assert!((area - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_signed_area_cw_square() {
        let sq = Polygon2D {
            vertices: vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]],
        };
        let area = sq.signed_area();
        assert!((area + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_classify_point_inside() {
        let sq = unit_square();
        assert_eq!(
            sq.classify_point([0.5, 0.5], 1e-8),
            PointClassification::Inside
        );
    }

    #[test]
    fn test_classify_point_outside() {
        let sq = unit_square();
        assert_eq!(
            sq.classify_point([2.0, 2.0], 1e-8),
            PointClassification::Outside
        );
    }

    #[test]
    fn test_classify_point_on_boundary() {
        let sq = unit_square();
        assert_eq!(
            sq.classify_point([0.5, 0.0], 1e-8),
            PointClassification::OnBoundary
        );
    }

    #[test]
    fn test_classify_point_vertex() {
        let sq = unit_square();
        assert_eq!(
            sq.classify_point([0.0, 0.0], 1e-8),
            PointClassification::OnBoundary
        );
    }

    #[test]
    fn test_intersect_line_horizontal() {
        let sq = unit_square();
        let hits = sq.intersect_line([0.0, 0.5], [1.0, 0.0]);
        assert_eq!(hits.len(), 2);
        // Should hit left edge (x=0) and right edge (x=1)
        assert!((hits[0].t_line - 0.0).abs() < 1e-12);
        assert!((hits[1].t_line - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_intersect_line_diagonal() {
        let sq = unit_square();
        let hits = sq.intersect_line([-1.0, -1.0], [1.0, 1.0]);
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn test_intersect_line_miss() {
        let sq = unit_square();
        let hits = sq.intersect_line([0.0, 2.0], [1.0, 0.0]);
        assert_eq!(hits.len(), 0);
    }
}

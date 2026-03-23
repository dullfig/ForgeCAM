pub mod constraint;
pub mod solver;

use rustkernel_math::{orthonormal_basis, Point3, Vec3};

use constraint::Constraint;
use solver::SolveResult;

/// Error type for sketch operations.
#[derive(Debug, Clone)]
pub enum SketchError {
    /// The sketch lines do not form a single closed loop.
    NotClosed,
    /// The sketch has no lines.
    Empty,
}

impl std::fmt::Display for SketchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SketchError::NotClosed => write!(f, "sketch lines do not form a closed loop"),
            SketchError::Empty => write!(f, "sketch has no lines"),
        }
    }
}

impl std::error::Error for SketchError {}

/// A 2D parametric sketch on a workplane.
///
/// Points have 2 DOF each (x, y in the workplane frame).
/// Lines reference two points. Constraints are equations the solver enforces.
pub struct Sketch {
    points: Vec<(f64, f64)>,
    lines: Vec<(usize, usize)>,  // (start_point_idx, end_point_idx)
    constraints: Vec<Constraint>,
    // Workplane definition for lifting 2D → 3D.
    origin: Point3,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Sketch {
    /// Create a new sketch on the given workplane.
    /// `origin` is the 3D origin; `normal` defines the workplane orientation.
    /// The x_axis and y_axis are computed automatically from the normal.
    pub fn new(origin: Point3, normal: Vec3) -> Self {
        let n = normal.normalize();
        let (x_axis, y_axis) = orthonormal_basis(&n);
        Self {
            points: Vec::new(),
            lines: Vec::new(),
            constraints: Vec::new(),
            origin,
            x_axis,
            y_axis,
        }
    }

    // ── Low-level API ──

    /// Add a point at (x, y) in the workplane. Returns the point index.
    pub fn add_point(&mut self, x: f64, y: f64) -> usize {
        let idx = self.points.len();
        self.points.push((x, y));
        idx
    }

    /// Add a line between two existing points. Returns the line index.
    pub fn add_line(&mut self, start: usize, end: usize) -> usize {
        let idx = self.lines.len();
        self.lines.push((start, end));
        idx
    }

    // ── High-level sketch tools (geometry + constraint in one call) ──

    /// Add a line between two points with a horizontal constraint.
    /// Returns the line index.
    pub fn add_horizontal_line(&mut self, start: usize, end: usize) -> usize {
        let line_idx = self.add_line(start, end);
        self.constraints.push(Constraint::Horizontal { line: line_idx });
        line_idx
    }

    /// Add a line between two points with a vertical constraint.
    /// Returns the line index.
    pub fn add_vertical_line(&mut self, start: usize, end: usize) -> usize {
        let line_idx = self.add_line(start, end);
        self.constraints.push(Constraint::Vertical { line: line_idx });
        line_idx
    }

    // ── Constraint API ──

    pub fn constrain_fixed(&mut self, point: usize, x: f64, y: f64) {
        self.constraints.push(Constraint::Fixed { point, x, y });
    }

    pub fn constrain_coincident(&mut self, p1: usize, p2: usize) {
        self.constraints.push(Constraint::Coincident { p1, p2 });
    }

    pub fn constrain_distance(&mut self, p1: usize, p2: usize, d: f64) {
        self.constraints.push(Constraint::Distance { p1, p2, d });
    }

    pub fn constrain_horizontal(&mut self, line: usize) {
        self.constraints.push(Constraint::Horizontal { line });
    }

    pub fn constrain_vertical(&mut self, line: usize) {
        self.constraints.push(Constraint::Vertical { line });
    }

    pub fn constrain_parallel(&mut self, l1: usize, l2: usize) {
        self.constraints.push(Constraint::Parallel { l1, l2 });
    }

    pub fn constrain_perpendicular(&mut self, l1: usize, l2: usize) {
        self.constraints.push(Constraint::Perpendicular { l1, l2 });
    }

    pub fn constrain_equal_length(&mut self, l1: usize, l2: usize) {
        self.constraints.push(Constraint::EqualLength { l1, l2 });
    }

    pub fn constrain_angle(&mut self, l1: usize, l2: usize, angle_rad: f64) {
        self.constraints.push(Constraint::Angle { l1, l2, angle_rad });
    }

    // ── Solver ──

    /// Run the Newton-Raphson constraint solver.
    /// On success, point positions are updated to satisfy all constraints.
    pub fn solve(&mut self) -> SolveResult {
        solver::solve(&mut self.points, &self.lines, &self.constraints)
    }

    /// Compute the current degrees of freedom without solving.
    pub fn dof(&self) -> usize {
        solver::compute_dof(&self.points, &self.lines, &self.constraints)
    }

    // ── Profile extraction ──

    /// Walk the sketch lines to find a closed loop, then lift each vertex
    /// to 3D via the workplane frame: `P_3d = origin + x * x_axis + y * y_axis`.
    ///
    /// Returns the 3D polygon vertices in loop order (suitable for extrude).
    pub fn to_profile_3d(&self) -> Result<Vec<Point3>, SketchError> {
        if self.lines.is_empty() {
            return Err(SketchError::Empty);
        }

        // Build adjacency: for each point, which lines start or end there?
        // We walk from the first line's start, always following the next line.
        let n_lines = self.lines.len();

        // Build a map from point → outgoing line(s) using (start, end) pairs.
        use std::collections::HashMap;
        let mut outgoing: HashMap<usize, Vec<(usize, usize)>> = HashMap::new(); // point → [(line_idx, dest_point)]
        for (li, &(s, e)) in self.lines.iter().enumerate() {
            outgoing.entry(s).or_default().push((li, e));
            outgoing.entry(e).or_default().push((li, s));
        }

        // Walk from the first line's start point.
        let start_point = self.lines[0].0;
        let mut current = start_point;
        let mut visited_lines = vec![false; n_lines];
        let mut loop_points = Vec::new();

        loop {
            loop_points.push(current);

            // Find an unvisited line from current.
            let next = outgoing.get(&current).and_then(|edges| {
                edges.iter().find(|(li, _)| !visited_lines[*li])
            });

            match next {
                Some(&(li, dest)) => {
                    visited_lines[li] = true;
                    current = dest;
                    if current == start_point {
                        // Closed the loop.
                        break;
                    }
                }
                None => {
                    // No unvisited line — not closed.
                    return Err(SketchError::NotClosed);
                }
            }
        }

        // Check that all lines were visited (single loop).
        if visited_lines.iter().any(|&v| !v) {
            return Err(SketchError::NotClosed);
        }

        // Lift 2D points to 3D.
        let profile_3d: Vec<Point3> = loop_points
            .iter()
            .map(|&pi| {
                let (x, y) = self.points[pi];
                self.origin + x * self.x_axis + y * self.y_axis
            })
            .collect();

        Ok(profile_3d)
    }

    /// Read-only access to point positions (for testing/inspection).
    pub fn point(&self, idx: usize) -> (f64, f64) {
        self.points[idx]
    }

    /// Number of points in the sketch.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Number of lines in the sketch.
    pub fn num_lines(&self) -> usize {
        self.lines.len()
    }

    /// Number of constraints in the sketch.
    pub fn num_constraints(&self) -> usize {
        self.constraints.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    // ── Solver unit tests ──

    #[test]
    fn test_fully_constrained_rectangle() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(2.0, 0.0);
        let p2 = sketch.add_point(2.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_horizontal_line(p0, p1);
        sketch.add_vertical_line(p1, p2);
        sketch.add_horizontal_line(p2, p3);
        sketch.add_vertical_line(p3, p0);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 2.0);
        sketch.constrain_distance(p1, p2, 1.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);
        assert_eq!(sketch.dof(), 0);

        // Verify positions.
        let (x0, y0) = sketch.point(p0);
        let (x1, y1) = sketch.point(p1);
        let (x2, y2) = sketch.point(p2);
        let (x3, y3) = sketch.point(p3);
        assert!((x0 - 0.0).abs() < 1e-8);
        assert!((y0 - 0.0).abs() < 1e-8);
        assert!((x1 - 2.0).abs() < 1e-8);
        assert!((y1 - 0.0).abs() < 1e-8);
        assert!((x2 - 2.0).abs() < 1e-8);
        assert!((y2 - 1.0).abs() < 1e-8);
        assert!((x3 - 0.0).abs() < 1e-8);
        assert!((y3 - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_under_constrained_reports_dof() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);

        sketch.add_horizontal_line(p0, p1);
        // Only horizontal constraint (1 eq) + 4 variables → 3 DOF

        let result = sketch.solve();
        assert_eq!(result, SolveResult::UnderConstrained { dof: 3 });
    }

    #[test]
    fn test_over_constrained_sketch() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);

        sketch.add_line(p0, p1);
        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_fixed(p1, 1.0, 0.0);
        // Already fully constrained (4 eqs, 4 vars). Add conflicting:
        sketch.constrain_distance(p0, p1, 5.0); // conflicts with fixed positions

        let result = sketch.solve();
        assert_eq!(result, SolveResult::OverConstrained);
    }

    #[test]
    fn test_horizontal_constraint_moves_point() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(3.0, 0.5); // slightly off-horizontal

        sketch.add_horizontal_line(p0, p1);
        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 3.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        let (_, y1) = sketch.point(p1);
        assert!(y1.abs() < 1e-8, "p1.y should be 0 (horizontal), got {y1}");
    }

    #[test]
    fn test_vertical_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(0.1, 2.0); // slightly off-vertical

        sketch.add_vertical_line(p0, p1);
        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 2.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        let (x1, _) = sketch.point(p1);
        assert!(x1.abs() < 1e-8, "p1.x should be 0 (vertical), got {x1}");
    }

    #[test]
    fn test_distance_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);

        sketch.add_line(p0, p1);
        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_fixed(p1, 5.0, 0.0);
        sketch.constrain_distance(p0, p1, 5.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        let (x1, y1) = sketch.point(p1);
        let dist = (x1 * x1 + y1 * y1).sqrt();
        assert!((dist - 5.0).abs() < 1e-8);
    }

    #[test]
    fn test_coincident_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(1.0, 2.0);
        let p1 = sketch.add_point(1.1, 2.1); // slightly off

        sketch.constrain_coincident(p0, p1);
        sketch.constrain_fixed(p0, 1.0, 2.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        let (x0, y0) = sketch.point(p0);
        let (x1, y1) = sketch.point(p1);
        assert!((x0 - x1).abs() < 1e-8);
        assert!((y0 - y1).abs() < 1e-8);
    }

    #[test]
    fn test_parallel_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(2.0, 0.0);
        let p2 = sketch.add_point(0.0, 1.0);
        let p3 = sketch.add_point(2.0, 1.2); // slightly non-parallel

        let l0 = sketch.add_line(p0, p1);
        let l1 = sketch.add_line(p2, p3);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_fixed(p1, 2.0, 0.0);
        sketch.constrain_fixed(p2, 0.0, 1.0);
        sketch.constrain_parallel(l0, l1);
        sketch.constrain_distance(p2, p3, 2.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        // Lines should be parallel: same direction vector.
        let (x2, y2) = sketch.point(p2);
        let (x3, y3) = sketch.point(p3);
        let dy = y3 - y2;
        assert!(dy.abs() < 1e-8, "Parallel lines should have same y: dy={dy}");
    }

    #[test]
    fn test_perpendicular_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);
        let p2 = sketch.add_point(0.0, 0.0);
        let p3 = sketch.add_point(0.0, 1.0);

        let l0 = sketch.add_line(p0, p1);
        let l1 = sketch.add_line(p2, p3);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_fixed(p1, 1.0, 0.0);
        sketch.constrain_fixed(p2, 0.0, 0.0);
        sketch.constrain_perpendicular(l0, l1);
        sketch.constrain_distance(p2, p3, 1.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        let (x0, y0) = sketch.point(p0);
        let (x1, y1) = sketch.point(p1);
        let (x2, y2) = sketch.point(p2);
        let (x3, y3) = sketch.point(p3);
        let dot = (x1 - x0) * (x3 - x2) + (y1 - y0) * (y3 - y2);
        assert!(dot.abs() < 1e-8, "Lines should be perpendicular, dot={dot}");
    }

    #[test]
    fn test_equal_length_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(3.0, 0.0);
        let p2 = sketch.add_point(0.0, 1.0);
        let p3 = sketch.add_point(0.0, 4.0); // length 3, will be adjusted

        let l0 = sketch.add_horizontal_line(p0, p1);
        let l1 = sketch.add_vertical_line(p2, p3);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 3.0);
        sketch.constrain_fixed(p2, 0.0, 1.0);
        sketch.constrain_equal_length(l0, l1);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        let (x0, y0) = sketch.point(p0);
        let (x1, y1) = sketch.point(p1);
        let (x2, y2) = sketch.point(p2);
        let (x3, y3) = sketch.point(p3);
        let len1 = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
        let len2 = ((x3 - x2).powi(2) + (y3 - y2).powi(2)).sqrt();
        assert!((len1 - len2).abs() < 1e-8, "Lengths should be equal: {len1} vs {len2}");
    }

    #[test]
    fn test_angle_constraint() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);
        let p2 = sketch.add_point(0.0, 0.0);
        let p3 = sketch.add_point(0.7, 0.7); // ~45 degrees

        let l0 = sketch.add_line(p0, p1);
        let l1 = sketch.add_line(p2, p3);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_fixed(p1, 1.0, 0.0);
        sketch.constrain_fixed(p2, 0.0, 0.0);
        sketch.constrain_angle(l0, l1, FRAC_PI_2); // 90 degrees
        sketch.constrain_distance(p2, p3, 1.0);

        let result = sketch.solve();
        assert_eq!(result, SolveResult::FullyConstrained);

        // l1 should be perpendicular to l0 (which is along +x),
        // so l1 direction should be along +y.
        let (x3, y3) = sketch.point(p3);
        assert!(x3.abs() < 1e-8, "x3 should be ~0, got {x3}");
        assert!((y3 - 1.0).abs() < 1e-8, "y3 should be ~1, got {y3}");
    }

    #[test]
    fn test_dof_unconstrained() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        sketch.add_point(0.0, 0.0);
        sketch.add_point(1.0, 0.0);
        // 2 points × 2 DOF each = 4 DOF
        assert_eq!(sketch.dof(), 4);
    }

    #[test]
    fn test_dof_after_fix() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let _p1 = sketch.add_point(1.0, 0.0);
        sketch.constrain_fixed(p0, 0.0, 0.0);
        // 4 vars - 2 fixed = 2 DOF
        assert_eq!(sketch.dof(), 2);
    }

    #[test]
    fn test_dof_rectangle() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(2.0, 0.0);
        let p2 = sketch.add_point(2.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_horizontal_line(p0, p1);
        sketch.add_vertical_line(p1, p2);
        sketch.add_horizontal_line(p2, p3);
        sketch.add_vertical_line(p3, p0);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 2.0);
        sketch.constrain_distance(p1, p2, 1.0);

        // 8 vars, 4 H/V constraints + 2 fixed + 2 distance = 8 eq → 0 DOF
        assert_eq!(sketch.dof(), 0);
    }

    // ── Profile extraction tests ──

    #[test]
    fn test_closed_rectangle_profile() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(2.0, 0.0);
        let p2 = sketch.add_point(2.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p0);

        let profile = sketch.to_profile_3d().unwrap();
        assert_eq!(profile.len(), 4);
    }

    #[test]
    fn test_l_shaped_profile() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(2.0, 0.0);
        let p2 = sketch.add_point(2.0, 1.0);
        let p3 = sketch.add_point(1.0, 1.0);
        let p4 = sketch.add_point(1.0, 2.0);
        let p5 = sketch.add_point(0.0, 2.0);

        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p4);
        sketch.add_line(p4, p5);
        sketch.add_line(p5, p0);

        let profile = sketch.to_profile_3d().unwrap();
        assert_eq!(profile.len(), 6);
    }

    #[test]
    fn test_open_sketch_error() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);
        let p2 = sketch.add_point(1.0, 1.0);

        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        // Not closed: p2 ≠ p0

        let result = sketch.to_profile_3d();
        assert!(result.is_err());
    }

    #[test]
    fn test_profile_3d_on_xy_plane() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);
        let p2 = sketch.add_point(1.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p0);

        let profile = sketch.to_profile_3d().unwrap();
        // All z should be 0 (on XY plane).
        for p in &profile {
            assert!(p.z.abs() < 1e-10, "Expected z=0, got {}", p.z);
        }
    }

    #[test]
    fn test_profile_3d_on_offset_plane() {
        let origin = Point3::new(0.0, 0.0, 5.0);
        let mut sketch = Sketch::new(origin, Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);
        let p2 = sketch.add_point(1.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p0);

        let profile = sketch.to_profile_3d().unwrap();
        for p in &profile {
            assert!((p.z - 5.0).abs() < 1e-10, "Expected z=5, got {}", p.z);
        }
    }

    #[test]
    fn test_solve_then_profile() {
        let mut sketch = Sketch::new(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));
        let p0 = sketch.add_point(0.1, -0.1); // slightly off
        let p1 = sketch.add_point(2.1, 0.1);
        let p2 = sketch.add_point(1.9, 0.9);
        let p3 = sketch.add_point(-0.1, 1.1);

        sketch.add_horizontal_line(p0, p1);
        sketch.add_vertical_line(p1, p2);
        sketch.add_horizontal_line(p2, p3);
        sketch.add_vertical_line(p3, p0);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 2.0);
        sketch.constrain_distance(p1, p2, 1.0);

        assert_eq!(sketch.solve(), SolveResult::FullyConstrained);

        // Verify solved 2D positions.
        let tol = 1e-6;
        let (x0, y0) = sketch.point(p0);
        let (x1, y1) = sketch.point(p1);
        let (x2, y2) = sketch.point(p2);
        let (x3, y3) = sketch.point(p3);
        assert!((x0 - 0.0).abs() < tol && (y0 - 0.0).abs() < tol);
        assert!((x1 - 2.0).abs() < tol && (y1 - 0.0).abs() < tol);
        assert!((x2 - 2.0).abs() < tol && (y2 - 1.0).abs() < tol);
        assert!((x3 - 0.0).abs() < tol && (y3 - 1.0).abs() < tol);

        let profile = sketch.to_profile_3d().unwrap();
        assert_eq!(profile.len(), 4);

        // All z should be 0 (on the XY workplane).
        for p in &profile {
            assert!(p.z.abs() < tol, "Expected z=0, got {}", p.z);
        }

        // Edge lengths in 3D should match the 2D rectangle (2×1).
        let mut edge_lengths: Vec<f64> = (0..4)
            .map(|i| (profile[(i + 1) % 4] - profile[i]).norm())
            .collect();
        edge_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((edge_lengths[0] - 1.0).abs() < tol);
        assert!((edge_lengths[1] - 1.0).abs() < tol);
        assert!((edge_lengths[2] - 2.0).abs() < tol);
        assert!((edge_lengths[3] - 2.0).abs() < tol);
    }
}

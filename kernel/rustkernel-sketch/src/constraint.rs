/// Constraint types for 2D parametric sketching.
///
/// Each constraint removes one or more degrees of freedom from the sketch.
/// All constraints provide analytical Jacobian entries for the Newton-Raphson solver.
#[derive(Debug, Clone)]
pub enum Constraint {
    /// Fix a point at (x, y). Removes 2 DOF.
    Fixed { point: usize, x: f64, y: f64 },
    /// Two points coincide. Removes 2 DOF.
    Coincident { p1: usize, p2: usize },
    /// Distance between two points equals d. Removes 1 DOF.
    /// Uses squared-distance formulation to avoid sqrt.
    Distance { p1: usize, p2: usize, d: f64 },
    /// Line is horizontal (start.y == end.y). Removes 1 DOF.
    Horizontal { line: usize },
    /// Line is vertical (start.x == end.x). Removes 1 DOF.
    Vertical { line: usize },
    /// Two lines are parallel (cross product of directions = 0). Removes 1 DOF.
    Parallel { l1: usize, l2: usize },
    /// Two lines are perpendicular (dot product of directions = 0). Removes 1 DOF.
    Perpendicular { l1: usize, l2: usize },
    /// Two lines have equal length. Removes 1 DOF.
    /// Uses squared-length formulation.
    EqualLength { l1: usize, l2: usize },
    /// Angle between two lines equals theta. Removes 1 DOF.
    /// Uses cross*cos(theta) - dot*sin(theta) = 0 to avoid atan2.
    Angle { l1: usize, l2: usize, angle_rad: f64 },
}

/// A (row, col, value) entry in the Jacobian matrix.
pub struct JacobianEntry {
    pub row: usize,
    pub col: usize,
    pub value: f64,
}

impl Constraint {
    /// Number of scalar equations this constraint contributes.
    pub fn num_equations(&self) -> usize {
        match self {
            Constraint::Fixed { .. } | Constraint::Coincident { .. } => 2,
            _ => 1,
        }
    }

    /// Evaluate the constraint residual(s) given the variable vector.
    /// `vars` is laid out as [x0, y0, x1, y1, ...] for points.
    /// `lines` maps line index to (start_point_idx, end_point_idx).
    pub fn evaluate(&self, vars: &[f64], lines: &[(usize, usize)]) -> Vec<f64> {
        match *self {
            Constraint::Fixed { point, x, y } => {
                let px = vars[point * 2];
                let py = vars[point * 2 + 1];
                vec![px - x, py - y]
            }
            Constraint::Coincident { p1, p2 } => {
                let x1 = vars[p1 * 2];
                let y1 = vars[p1 * 2 + 1];
                let x2 = vars[p2 * 2];
                let y2 = vars[p2 * 2 + 1];
                vec![x2 - x1, y2 - y1]
            }
            Constraint::Distance { p1, p2, d } => {
                let dx = vars[p2 * 2] - vars[p1 * 2];
                let dy = vars[p2 * 2 + 1] - vars[p1 * 2 + 1];
                vec![dx * dx + dy * dy - d * d]
            }
            Constraint::Horizontal { line } => {
                let (s, e) = lines[line];
                let ys = vars[s * 2 + 1];
                let ye = vars[e * 2 + 1];
                vec![ye - ys]
            }
            Constraint::Vertical { line } => {
                let (s, e) = lines[line];
                let xs = vars[s * 2];
                let xe = vars[e * 2];
                vec![xe - xs]
            }
            Constraint::Parallel { l1, l2 } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                // cross product: dx1*dy2 - dy1*dx2 = 0
                vec![dx1 * dy2 - dy1 * dx2]
            }
            Constraint::Perpendicular { l1, l2 } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                // dot product: dx1*dx2 + dy1*dy2 = 0
                vec![dx1 * dx2 + dy1 * dy2]
            }
            Constraint::EqualLength { l1, l2 } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                vec![dx1 * dx1 + dy1 * dy1 - dx2 * dx2 - dy2 * dy2]
            }
            Constraint::Angle { l1, l2, angle_rad } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                let cross = dx1 * dy2 - dy1 * dx2;
                let dot = dx1 * dx2 + dy1 * dy2;
                let cos_t = angle_rad.cos();
                let sin_t = angle_rad.sin();
                vec![cross * cos_t - dot * sin_t]
            }
        }
    }

    /// Compute analytical Jacobian entries for this constraint.
    /// `row_offset` is the starting row in the global Jacobian.
    pub fn jacobian_entries(
        &self,
        vars: &[f64],
        lines: &[(usize, usize)],
        row_offset: usize,
    ) -> Vec<JacobianEntry> {
        let mut entries = Vec::new();
        match *self {
            Constraint::Fixed { point, .. } => {
                let cx = point * 2;
                let cy = point * 2 + 1;
                // d(px - x)/dpx = 1
                entries.push(JacobianEntry { row: row_offset, col: cx, value: 1.0 });
                // d(py - y)/dpy = 1
                entries.push(JacobianEntry { row: row_offset + 1, col: cy, value: 1.0 });
            }
            Constraint::Coincident { p1, p2 } => {
                let x1c = p1 * 2;
                let y1c = p1 * 2 + 1;
                let x2c = p2 * 2;
                let y2c = p2 * 2 + 1;
                // d(x2-x1)/dx1 = -1, d(x2-x1)/dx2 = 1
                entries.push(JacobianEntry { row: row_offset, col: x1c, value: -1.0 });
                entries.push(JacobianEntry { row: row_offset, col: x2c, value: 1.0 });
                // d(y2-y1)/dy1 = -1, d(y2-y1)/dy2 = 1
                entries.push(JacobianEntry { row: row_offset + 1, col: y1c, value: -1.0 });
                entries.push(JacobianEntry { row: row_offset + 1, col: y2c, value: 1.0 });
            }
            Constraint::Distance { p1, p2, .. } => {
                let dx = vars[p2 * 2] - vars[p1 * 2];
                let dy = vars[p2 * 2 + 1] - vars[p1 * 2 + 1];
                // f = dx²+dy²-d²  →  df/dx1=-2dx, df/dy1=-2dy, df/dx2=2dx, df/dy2=2dy
                entries.push(JacobianEntry { row: row_offset, col: p1 * 2, value: -2.0 * dx });
                entries.push(JacobianEntry { row: row_offset, col: p1 * 2 + 1, value: -2.0 * dy });
                entries.push(JacobianEntry { row: row_offset, col: p2 * 2, value: 2.0 * dx });
                entries.push(JacobianEntry { row: row_offset, col: p2 * 2 + 1, value: 2.0 * dy });
            }
            Constraint::Horizontal { line } => {
                let (s, e) = lines[line];
                // f = ye - ys  →  df/dys = -1, df/dye = 1
                entries.push(JacobianEntry { row: row_offset, col: s * 2 + 1, value: -1.0 });
                entries.push(JacobianEntry { row: row_offset, col: e * 2 + 1, value: 1.0 });
            }
            Constraint::Vertical { line } => {
                let (s, e) = lines[line];
                // f = xe - xs  →  df/dxs = -1, df/dxe = 1
                entries.push(JacobianEntry { row: row_offset, col: s * 2, value: -1.0 });
                entries.push(JacobianEntry { row: row_offset, col: e * 2, value: 1.0 });
            }
            Constraint::Parallel { l1, l2 } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                // f = dx1*dy2 - dy1*dx2
                // df/d(xs1) = -dy2, df/d(ys1) = dx2
                // df/d(xe1) = dy2,  df/d(ye1) = -dx2
                // df/d(xs2) = dy1,  df/d(ys2) = -dx1
                // df/d(xe2) = -dy1, df/d(ye2) = dx1
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2, value: -dy2 });
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2 + 1, value: dx2 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2, value: dy2 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2 + 1, value: -dx2 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2, value: dy1 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2 + 1, value: -dx1 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2, value: -dy1 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2 + 1, value: dx1 });
            }
            Constraint::Perpendicular { l1, l2 } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                // f = dx1*dx2 + dy1*dy2
                // df/d(xs1) = -dx2, df/d(ys1) = -dy2
                // df/d(xe1) = dx2,  df/d(ye1) = dy2
                // df/d(xs2) = -dx1, df/d(ys2) = -dy1
                // df/d(xe2) = dx1,  df/d(ye2) = dy1
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2, value: -dx2 });
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2 + 1, value: -dy2 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2, value: dx2 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2 + 1, value: dy2 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2, value: -dx1 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2 + 1, value: -dy1 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2, value: dx1 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2 + 1, value: dy1 });
            }
            Constraint::EqualLength { l1, l2 } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                // f = dx1²+dy1² - dx2²-dy2²
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2, value: -2.0 * dx1 });
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2 + 1, value: -2.0 * dy1 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2, value: 2.0 * dx1 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2 + 1, value: 2.0 * dy1 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2, value: 2.0 * dx2 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2 + 1, value: 2.0 * dy2 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2, value: -2.0 * dx2 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2 + 1, value: -2.0 * dy2 });
            }
            Constraint::Angle { l1, l2, angle_rad } => {
                let (s1, e1) = lines[l1];
                let (s2, e2) = lines[l2];
                let dx1 = vars[e1 * 2] - vars[s1 * 2];
                let dy1 = vars[e1 * 2 + 1] - vars[s1 * 2 + 1];
                let dx2 = vars[e2 * 2] - vars[s2 * 2];
                let dy2 = vars[e2 * 2 + 1] - vars[s2 * 2 + 1];
                let cos_t = angle_rad.cos();
                let sin_t = angle_rad.sin();
                // f = (dx1*dy2 - dy1*dx2)*cos_t - (dx1*dx2 + dy1*dy2)*sin_t
                // df/d(xs1): d_cross/d(xs1)*cos - d_dot/d(xs1)*sin
                //   d_cross/d(xs1) = -dy2, d_dot/d(xs1) = -dx2
                //   → -dy2*cos + dx2*sin
                // df/d(ys1): dx2*cos + dy2*sin → but let me be careful
                //   d_cross/d(ys1) = dx2, d_dot/d(ys1) = -dy2
                //   → dx2*cos + dy2*sin
                // ... and so on for all 8 partials.
                // Let's compute systematically:
                // cross = dx1*dy2 - dy1*dx2
                // dot = dx1*dx2 + dy1*dy2
                // f = cross*C - dot*S where C=cos_t, S=sin_t
                //
                // d(dx1)/d(xs1) = -1, d(dx1)/d(xe1) = 1
                // d(dy1)/d(ys1) = -1, d(dy1)/d(ye1) = 1
                // d(dx2)/d(xs2) = -1, d(dx2)/d(xe2) = 1
                // d(dy2)/d(ys2) = -1, d(dy2)/d(ye2) = 1
                //
                // df/d(xs1) = (d_cross/d_dx1 * C - d_dot/d_dx1 * S) * d_dx1/d_xs1
                //           = (dy2*C - dx2*S) * (-1)
                let dc_xs1 = -(dy2 * cos_t - dx2 * sin_t);
                let dc_ys1 = -(-dx2 * cos_t - dy2 * sin_t);
                let dc_xe1 = dy2 * cos_t - dx2 * sin_t;
                let dc_ye1 = -dx2 * cos_t - dy2 * sin_t;
                let dc_xs2 = -(- dy1 * cos_t - dx1 * sin_t);
                let dc_ys2 = -(dx1 * cos_t - dy1 * sin_t);
                let dc_xe2 = -dy1 * cos_t - dx1 * sin_t;
                let dc_ye2 = dx1 * cos_t - dy1 * sin_t;

                entries.push(JacobianEntry { row: row_offset, col: s1 * 2, value: dc_xs1 });
                entries.push(JacobianEntry { row: row_offset, col: s1 * 2 + 1, value: dc_ys1 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2, value: dc_xe1 });
                entries.push(JacobianEntry { row: row_offset, col: e1 * 2 + 1, value: dc_ye1 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2, value: dc_xs2 });
                entries.push(JacobianEntry { row: row_offset, col: s2 * 2 + 1, value: dc_ys2 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2, value: dc_xe2 });
                entries.push(JacobianEntry { row: row_offset, col: e2 * 2 + 1, value: dc_ye2 });
            }
        }
        entries
    }
}

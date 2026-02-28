use rustkernel_math::nalgebra::{DMatrix, DVector};

use crate::constraint::Constraint;

const MAX_ITERATIONS: usize = 50;
const TOLERANCE: f64 = 1e-10;
/// Singular values below this fraction of the largest are treated as zero for rank computation.
const SVD_RANK_TOL: f64 = 1e-8;

/// Result of the constraint solver.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolveResult {
    FullyConstrained,
    UnderConstrained { dof: usize },
    OverConstrained,
    DidNotConverge,
}

/// Newton-Raphson constraint solver with SVD-based rank analysis.
pub fn solve(
    points: &mut [(f64, f64)],
    lines: &[(usize, usize)],
    constraints: &[Constraint],
) -> SolveResult {
    let n_vars = points.len() * 2;
    let n_eqs: usize = constraints.iter().map(|c| c.num_equations()).sum();

    if n_vars == 0 {
        if n_eqs == 0 {
            return SolveResult::FullyConstrained;
        } else {
            return SolveResult::OverConstrained;
        }
    }

    if n_eqs == 0 {
        return SolveResult::UnderConstrained { dof: n_vars };
    }

    // Flatten point coordinates into variable vector.
    let mut vars = DVector::zeros(n_vars);
    for (i, &(x, y)) in points.iter().enumerate() {
        vars[i * 2] = x;
        vars[i * 2 + 1] = y;
    }

    let vars_slice = |v: &DVector<f64>| -> Vec<f64> { v.iter().copied().collect() };

    for _iter in 0..MAX_ITERATIONS {
        let vs = vars_slice(&vars);

        // Evaluate residuals.
        let mut residuals = DVector::zeros(n_eqs);
        let mut row = 0;
        for c in constraints {
            let vals = c.evaluate(&vs, lines);
            for &v in &vals {
                residuals[row] = v;
                row += 1;
            }
        }

        // Check convergence.
        if residuals.norm() < TOLERANCE {
            // Converged — determine DOF via Jacobian rank.
            let jacobian = build_jacobian(constraints, &vs, lines, n_eqs, n_vars);
            let rank = compute_rank(&jacobian);
            unflatten(points, &vars);

            return if rank > n_vars {
                // More independent equations than variables — shouldn't happen with correct rank
                SolveResult::OverConstrained
            } else {
                let dof = n_vars - rank;
                if dof == 0 {
                    SolveResult::FullyConstrained
                } else {
                    SolveResult::UnderConstrained { dof }
                }
            };
        }

        // Build Jacobian.
        let jacobian = build_jacobian(constraints, &vs, lines, n_eqs, n_vars);

        // Solve J * dx = -r via SVD pseudoinverse.
        let svd = jacobian.svd(true, true);
        let dx = svd.solve(&(-&residuals), TOLERANCE).unwrap_or_else(|_| DVector::zeros(n_vars));

        vars += dx;
    }

    // Did not converge — check if we're close enough and over-constrained.
    let vs = vars_slice(&vars);
    let mut residuals = DVector::zeros(n_eqs);
    let mut row = 0;
    for c in constraints {
        let vals = c.evaluate(&vs, lines);
        for &v in &vals {
            residuals[row] = v;
            row += 1;
        }
    }

    // If residual is large and n_eqs > n_vars, likely over-constrained.
    if residuals.norm() > 1e-4 && n_eqs > n_vars {
        return SolveResult::OverConstrained;
    }

    // Check if it's over-constrained by examining Jacobian rank vs equation count.
    let jacobian = build_jacobian(constraints, &vs, lines, n_eqs, n_vars);
    let rank = compute_rank(&jacobian);
    if residuals.norm() > 1e-4 && rank >= n_vars {
        return SolveResult::OverConstrained;
    }

    SolveResult::DidNotConverge
}

/// Compute DOF without modifying point positions — just analyzes the Jacobian.
pub fn compute_dof(
    points: &[(f64, f64)],
    lines: &[(usize, usize)],
    constraints: &[Constraint],
) -> usize {
    let n_vars = points.len() * 2;
    let n_eqs: usize = constraints.iter().map(|c| c.num_equations()).sum();

    if n_vars == 0 || n_eqs == 0 {
        return n_vars;
    }

    let mut vars = Vec::with_capacity(n_vars);
    for &(x, y) in points {
        vars.push(x);
        vars.push(y);
    }

    let jacobian = build_jacobian(constraints, &vars, lines, n_eqs, n_vars);
    let rank = compute_rank(&jacobian);
    n_vars.saturating_sub(rank)
}

fn build_jacobian(
    constraints: &[Constraint],
    vars: &[f64],
    lines: &[(usize, usize)],
    n_rows: usize,
    n_cols: usize,
) -> DMatrix<f64> {
    let mut jac = DMatrix::zeros(n_rows, n_cols);
    let mut row = 0;
    for c in constraints {
        let entries = c.jacobian_entries(vars, lines, row);
        for e in entries {
            if e.col < n_cols && e.row < n_rows {
                jac[(e.row, e.col)] += e.value;
            }
        }
        row += c.num_equations();
    }
    jac
}

fn compute_rank(m: &DMatrix<f64>) -> usize {
    let svd = m.clone().svd(false, false);
    let svals = &svd.singular_values;
    if svals.is_empty() {
        return 0;
    }
    let max_sv = svals[0];
    if max_sv < 1e-15 {
        return 0;
    }
    let threshold = max_sv * SVD_RANK_TOL;
    svals.iter().filter(|&&s| s > threshold).count()
}

fn unflatten(points: &mut [(f64, f64)], vars: &DVector<f64>) {
    for (i, pt) in points.iter_mut().enumerate() {
        pt.0 = vars[i * 2];
        pt.1 = vars[i * 2 + 1];
    }
}

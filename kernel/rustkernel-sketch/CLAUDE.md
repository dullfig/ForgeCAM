# rustkernel-sketch — Sub-Agent

You own the 2D parametric sketch constraint solver.

## Your files

| File | Lines | What it does |
|------|-------|--------------|
| `lib.rs` | ~670 | Sketch struct, workplane, profile extraction, all tests |
| `constraint.rs` | — | Constraint enum (9 types) with analytical Jacobians |
| `solver.rs` | — | Newton-Raphson solver with SVD pseudoinverse |

## Architecture

### Sketch
- Points have 2 DOF each (x, y in workplane frame)
- Lines reference two point indices
- Constraints are equations the solver enforces
- Workplane: origin + normal → auto-computed x_axis/y_axis

### Constraint types
Fixed, Coincident, Distance, Horizontal, Vertical, Parallel, Perpendicular, EqualLength, Angle

### Solver (Newton-Raphson + SVD)
- Builds residual vector + Jacobian matrix each iteration
- SVD pseudoinverse handles rank-deficient systems
- Reports: FullyConstrained, UnderConstrained{dof}, OverConstrained

### Profile extraction
`to_profile_3d()`: walk lines to find single closed loop, lift 2D→3D via workplane.
Output feeds directly into `extrude_builder` or `revolve_builder`.

## Key design rule

**All constraint solving is 2D.** Never do 2D operations in 3D space. Project to
workplane first, solve, then lift to 3D only via `to_profile_3d()`.

## Dependencies

- `rustkernel_math` — Point3, Vec3, orthonormal_basis
- `nalgebra` — DMatrix, DVector, SVD

## Future work

- Arc/circle sketch elements (not just lines)
- Tangent constraints
- Symmetry constraints
- Mirror/pattern within sketch

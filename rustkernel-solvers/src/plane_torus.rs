use rustkernel_math::{orthonormal_basis, Point3, Vec3};
use rustkernel_topology::geom_store::SurfaceKind;
use rustkernel_topology::intersection::{
    IntersectionCircle, IntersectionCurve, IntersectionError, SurfaceSurfaceResult,
    SurfaceSurfaceSolver, INTERSECTION_TOLERANCE,
};

/// Plane-torus intersection solver (special cases only).
///
/// Handles:
/// - Plane perpendicular to torus axis: two concentric circles (outer R+r, inner R-r)
/// - General case: quartic curve — returns NoSolverAvailable (deferred to Phase 4)
pub struct PlaneTorusSolver;

impl PlaneTorusSolver {
    fn solve_impl(
        plane_origin: Point3,
        plane_normal: Vec3,
        torus_center: Point3,
        torus_axis: Vec3,
        major_r: f64,
        minor_r: f64,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let n = plane_normal.normalize();
        let a = torus_axis.normalize();
        let r = minor_r.abs();
        let big_r = major_r;

        let cos_na = n.dot(&a).abs();

        if cos_na > 1.0 - INTERSECTION_TOLERANCE {
            // Plane perpendicular to torus axis.
            // Distance from torus center to plane along axis.
            let d = n.dot(&(plane_origin - torus_center));

            if d.abs() > r + INTERSECTION_TOLERANCE {
                return Ok(SurfaceSurfaceResult::Empty);
            }

            if d.abs() < INTERSECTION_TOLERANCE {
                // Plane through torus center perpendicular to axis.
                // Two concentric circles: R+r and R-r.
                let (ref_dir, _) = orthonormal_basis(&n);

                let mut curves = Vec::new();
                curves.push(IntersectionCurve::Circle(IntersectionCircle {
                    center: torus_center,
                    axis: n,
                    radius: big_r + r,
                    ref_dir,
                }));
                if big_r > r + INTERSECTION_TOLERANCE {
                    curves.push(IntersectionCurve::Circle(IntersectionCircle {
                        center: torus_center,
                        axis: n,
                        radius: big_r - r,
                        ref_dir,
                    }));
                }
                return Ok(SurfaceSurfaceResult::Curves(curves));
            }

            // Offset perpendicular plane: still two circles but with adjusted radii.
            let rho = (r * r - d * d).sqrt(); // cross-section radius at height d
            let circle_center = torus_center + d * a;
            let (ref_dir, _) = orthonormal_basis(&n);

            let mut curves = Vec::new();
            curves.push(IntersectionCurve::Circle(IntersectionCircle {
                center: circle_center,
                axis: n,
                radius: big_r + rho,
                ref_dir,
            }));
            if big_r > rho + INTERSECTION_TOLERANCE {
                curves.push(IntersectionCurve::Circle(IntersectionCircle {
                    center: circle_center,
                    axis: n,
                    radius: big_r - rho,
                    ref_dir,
                }));
            }
            return Ok(SurfaceSurfaceResult::Curves(curves));
        }

        // General case: quartic Villarceau or spiric curve — defer to Phase 4.
        Err(IntersectionError::NoSolverAvailable)
    }
}

impl SurfaceSurfaceSolver for PlaneTorusSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane { .. }, SurfaceKind::Torus { .. })
                | (SurfaceKind::Torus { .. }, SurfaceKind::Plane { .. })
        )
    }

    fn solve(
        &self,
        a: &SurfaceKind,
        b: &SurfaceKind,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        match (a, b) {
            (
                SurfaceKind::Plane { origin, normal },
                SurfaceKind::Torus {
                    center,
                    axis,
                    major_radius,
                    minor_radius,
                },
            ) => Self::solve_impl(*origin, *normal, *center, *axis, *major_radius, *minor_radius),
            (
                SurfaceKind::Torus {
                    center,
                    axis,
                    major_radius,
                    minor_radius,
                },
                SurfaceKind::Plane { origin, normal },
            ) => Self::solve_impl(*origin, *normal, *center, *axis, *major_radius, *minor_radius),
            _ => Err(IntersectionError::InvalidInput(
                "expected plane and torus".into(),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    #[test]
    fn equatorial_cut() {
        // Plane through torus center perpendicular to axis.
        // R=5, r=2 → circles at radii 7 and 3.
        let result = PlaneTorusSolver::solve_impl(
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            5.0,
            2.0,
        )
        .unwrap();

        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 2);
                let IntersectionCurve::Circle(c0) = &curves[0] else {
                    panic!("expected Circle");
                };
                let IntersectionCurve::Circle(c1) = &curves[1] else {
                    panic!("expected Circle");
                };
                let r0 = c0.radius;
                let r1 = c1.radius;
                let (big, small) = if r0 > r1 { (r0, r1) } else { (r1, r0) };
                assert!((big - 7.0).abs() < TOL);
                assert!((small - 3.0).abs() < TOL);
            }
            other => panic!("expected 2 circles, got {other:?}"),
        }
    }

    #[test]
    fn miss() {
        let result = PlaneTorusSolver::solve_impl(
            Point3::new(0.0, 0.0, 10.0),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            5.0,
            2.0,
        )
        .unwrap();
        assert!(matches!(result, SurfaceSurfaceResult::Empty));
    }

    #[test]
    fn general_case_deferred() {
        // Non-perpendicular cut — should return NoSolverAvailable.
        let result = PlaneTorusSolver::solve_impl(
            Point3::origin(),
            Vec3::new(1.0, 0.0, 0.0), // perpendicular to torus axis → general case
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            5.0,
            2.0,
        );
        assert!(matches!(result, Err(IntersectionError::NoSolverAvailable)));
    }

    #[test]
    fn offset_perpendicular_cut() {
        // Plane at z=1 perpendicular to axis. R=5, r=2.
        // Cross-section at height 1: rho = sqrt(4-1) = sqrt(3).
        // Circles at radii 5+sqrt(3) and 5-sqrt(3).
        let result = PlaneTorusSolver::solve_impl(
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            5.0,
            2.0,
        )
        .unwrap();

        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 2);
                let IntersectionCurve::Circle(c0) = &curves[0] else {
                    panic!("expected Circle");
                };
                let IntersectionCurve::Circle(c1) = &curves[1] else {
                    panic!("expected Circle");
                };
                let rho = (3.0_f64).sqrt();
                let expected_big = 5.0 + rho;
                let expected_small = 5.0 - rho;
                let r0 = c0.radius;
                let r1 = c1.radius;
                let (big, small) = if r0 > r1 { (r0, r1) } else { (r1, r0) };
                assert!((big - expected_big).abs() < TOL);
                assert!((small - expected_small).abs() < TOL);
            }
            other => panic!("expected 2 circles, got {other:?}"),
        }
    }
}

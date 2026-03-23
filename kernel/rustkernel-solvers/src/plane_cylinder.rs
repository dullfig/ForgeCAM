use rustkernel_math::{orthonormal_basis, Point3, Vec3};
use rustkernel_topology::geom_store::SurfaceKind;
use rustkernel_topology::intersection::{
    IntersectionCircle, IntersectionCurve, IntersectionEllipse, IntersectionError,
    SurfaceSurfaceResult, SurfaceSurfaceSolver, INTERSECTION_TOLERANCE,
};
use tracing::debug;

/// Analytical plane-cylinder intersection solver.
///
/// Depending on the angle between the plane normal and the cylinder axis:
/// - Perpendicular to axis (cos_theta ≈ 1): Circle
/// - General angle: Ellipse
/// - Parallel to axis (cos_theta ≈ 0): 0 or 2 lines (not implemented — returns Empty or NoSolverAvailable)
pub struct PlaneCylinderSolver;

impl PlaneCylinderSolver {
    fn solve_impl(
        plane_origin: Point3,
        plane_normal: Vec3,
        cyl_origin: Point3,
        cyl_axis: Vec3,
        cyl_radius: f64,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let r = cyl_radius.abs();
        let n = plane_normal.normalize();
        let a = cyl_axis.normalize();

        let cos_theta = n.dot(&a).abs();

        if cos_theta > 1.0 - INTERSECTION_TOLERANCE {
            // Plane perpendicular to cylinder axis → Circle.
            // Project cylinder axis onto plane to find circle center.
            let d = n.dot(&(plane_origin - cyl_origin));
            let t = d / n.dot(&a);
            let circle_center = cyl_origin + t * a;

            let (ref_dir, _) = orthonormal_basis(&n);

            return Ok(SurfaceSurfaceResult::Curves(vec![
                IntersectionCurve::Circle(IntersectionCircle {
                    center: circle_center,
                    axis: n,
                    radius: r,
                    ref_dir,
                }),
            ]));
        }

        if cos_theta < INTERSECTION_TOLERANCE {
            // Plane parallel to cylinder axis.
            // Distance from cylinder axis to plane determines the result.
            // This produces 0 or 2 lines — defer to Phase 4.
            let axis_to_plane = plane_origin - cyl_origin;
            let axis_component = axis_to_plane.dot(&a) * a;
            let perp = axis_to_plane - axis_component;
            let dist = n.dot(&perp);

            if dist.abs() > r + INTERSECTION_TOLERANCE {
                return Ok(SurfaceSurfaceResult::Empty);
            }
            // Tangent or two lines — not implemented yet.
            return Err(IntersectionError::NoSolverAvailable);
        }

        // General case: Ellipse.
        // The plane cuts the cylinder at an oblique angle.
        // semi_minor = r, semi_major = r / cos_theta.
        let semi_minor = r;
        let semi_major = r / cos_theta;

        // Find the center: project cylinder axis onto the plane.
        // The center is where the cylinder axis pierces the plane.
        let denom = n.dot(&a);
        let t = n.dot(&(plane_origin - cyl_origin)) / denom;
        let ellipse_center = cyl_origin + t * a;

        // major_dir: projection of cylinder axis onto the plane, normalized.
        let a_proj = a - n.dot(&a) * n;
        let major_dir = a_proj.normalize();

        Ok(SurfaceSurfaceResult::Curves(vec![
            IntersectionCurve::Ellipse(IntersectionEllipse {
                center: ellipse_center,
                axis: n,
                semi_major,
                semi_minor,
                major_dir,
            }),
        ]))
    }
}

impl SurfaceSurfaceSolver for PlaneCylinderSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane { .. }, SurfaceKind::Cylinder { .. })
                | (SurfaceKind::Cylinder { .. }, SurfaceKind::Plane { .. })
        )
    }

    fn solve(
        &self,
        a: &SurfaceKind,
        b: &SurfaceKind,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        debug!("plane_cylinder solve");
        match (a, b) {
            (
                SurfaceKind::Plane { origin, normal },
                SurfaceKind::Cylinder {
                    origin: co,
                    axis,
                    radius,
                },
            ) => Self::solve_impl(*origin, *normal, *co, *axis, *radius),
            (
                SurfaceKind::Cylinder {
                    origin: co,
                    axis,
                    radius,
                },
                SurfaceKind::Plane { origin, normal },
            ) => Self::solve_impl(*origin, *normal, *co, *axis, *radius),
            _ => Err(IntersectionError::InvalidInput(
                "expected plane and cylinder".into(),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    #[test]
    fn perpendicular_cut_circle() {
        // Plane z=2 cutting vertical cylinder of radius 3.
        let result = PlaneCylinderSolver::solve_impl(
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            3.0,
        )
        .unwrap();

        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                let IntersectionCurve::Circle(c) = &curves[0] else {
                    panic!("expected Circle");
                };
                assert!((c.radius - 3.0).abs() < TOL);
                assert!((c.center.z - 2.0).abs() < TOL);
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    #[test]
    fn oblique_cut_ellipse() {
        // Plane with normal (0, 0.5, 0.866) ≈ 30° from Z axis, cutting vertical cylinder r=1.
        let normal = Vec3::new(0.0, 0.5, (3.0_f64).sqrt() / 2.0).normalize();
        let result = PlaneCylinderSolver::solve_impl(
            Point3::origin(),
            normal,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
        )
        .unwrap();

        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                let IntersectionCurve::Ellipse(e) = &curves[0] else {
                    panic!("expected Ellipse");
                };
                assert!((e.semi_minor - 1.0).abs() < TOL);
                assert!(e.semi_major > 1.0 + TOL, "semi_major should be > 1");
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    #[test]
    fn parallel_far_empty() {
        // Plane parallel to cylinder axis, far away.
        let result = PlaneCylinderSolver::solve_impl(
            Point3::new(10.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
        );
        match result {
            Ok(SurfaceSurfaceResult::Empty) => {}
            other => panic!("expected Empty, got {other:?}"),
        }
    }

    #[test]
    fn circle_points_satisfy_both_surfaces() {
        let plane_origin = Point3::new(0.0, 0.0, 5.0);
        let plane_normal = Vec3::new(0.0, 0.0, 1.0);
        let cyl_origin = Point3::origin();
        let cyl_axis = Vec3::new(0.0, 0.0, 1.0);
        let cyl_radius = 3.0;

        let result = PlaneCylinderSolver::solve_impl(
            plane_origin, plane_normal, cyl_origin, cyl_axis, cyl_radius,
        )
        .unwrap();

        if let SurfaceSurfaceResult::Curves(curves) = result {
            let IntersectionCurve::Circle(c) = &curves[0] else {
                panic!("expected Circle");
            };
            let ref_y = c.axis.cross(&c.ref_dir);
            for i in 0..8 {
                let angle = (i as f64 / 8.0) * std::f64::consts::TAU;
                let p = c.center + c.radius * (angle.cos() * c.ref_dir + angle.sin() * ref_y);
                // On plane?
                let dist_plane = plane_normal.dot(&(p - plane_origin));
                assert!(dist_plane.abs() < TOL);
                // On cylinder? distance from axis = radius.
                let to_axis = p - cyl_origin;
                let along = to_axis.dot(&cyl_axis.normalize());
                let perp = (to_axis - along * cyl_axis.normalize()).norm();
                assert!((perp - cyl_radius).abs() < TOL);
            }
        }
    }
}

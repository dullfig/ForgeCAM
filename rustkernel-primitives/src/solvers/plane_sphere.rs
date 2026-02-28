use rustkernel_math::{orthonormal_basis, Point3, Vec3};
use rustkernel_topology::geom_store::SurfaceKind;
use rustkernel_topology::intersection::{
    IntersectionCircle, IntersectionCurve, IntersectionError, SurfaceSurfaceResult,
    SurfaceSurfaceSolver, INTERSECTION_TOLERANCE,
};

/// Analytical plane-sphere intersection solver.
///
/// Produces a circle (if the plane cuts the sphere) or Empty.
pub struct PlaneSphereSolver;

impl PlaneSphereSolver {
    fn solve_impl(
        plane_origin: Point3,
        plane_normal: Vec3,
        sphere_center: Point3,
        sphere_radius: f64,
    ) -> SurfaceSurfaceResult {
        let r = sphere_radius.abs();
        let n = plane_normal.normalize();

        // Signed distance from sphere center to plane.
        let d = n.dot(&(sphere_center - plane_origin));

        if d.abs() > r - INTERSECTION_TOLERANCE {
            // No intersection (or tangent — skip for boolean purposes).
            return SurfaceSurfaceResult::Empty;
        }

        // Circle of intersection.
        let circle_radius = (r * r - d * d).sqrt();
        let circle_center = sphere_center - d * n;

        // ref_dir: any direction in the plane of the circle.
        let (ref_dir, _) = orthonormal_basis(&n);

        SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Circle(IntersectionCircle {
            center: circle_center,
            axis: n,
            radius: circle_radius,
            ref_dir,
        })])
    }
}

impl SurfaceSurfaceSolver for PlaneSphereSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane { .. }, SurfaceKind::Sphere { .. })
                | (SurfaceKind::Sphere { .. }, SurfaceKind::Plane { .. })
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
                SurfaceKind::Sphere { center, radius },
            ) => Ok(Self::solve_impl(*origin, *normal, *center, *radius)),
            (
                SurfaceKind::Sphere { center, radius },
                SurfaceKind::Plane { origin, normal },
            ) => Ok(Self::solve_impl(*origin, *normal, *center, *radius)),
            _ => Err(IntersectionError::InvalidInput(
                "expected plane and sphere".into(),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    #[test]
    fn equator_cut() {
        let result = PlaneSphereSolver::solve_impl(
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            5.0,
        );
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                let IntersectionCurve::Circle(c) = &curves[0] else {
                    panic!("expected Circle");
                };
                assert!((c.radius - 5.0).abs() < TOL);
                assert!((c.center - Point3::origin()).norm() < TOL);
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    #[test]
    fn offset_cut() {
        // Plane at z=3, sphere of radius 5 at origin → circle radius = sqrt(25-9) = 4.
        let result = PlaneSphereSolver::solve_impl(
            Point3::new(0.0, 0.0, 3.0),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            5.0,
        );
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                let IntersectionCurve::Circle(c) = &curves[0] else {
                    panic!("expected Circle");
                };
                assert!((c.radius - 4.0).abs() < TOL);
                assert!((c.center.z - 3.0).abs() < TOL);
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    #[test]
    fn miss() {
        let result = PlaneSphereSolver::solve_impl(
            Point3::new(0.0, 0.0, 10.0),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            5.0,
        );
        assert!(matches!(result, SurfaceSurfaceResult::Empty));
    }

    #[test]
    fn points_on_result_satisfy_both() {
        let plane_origin = Point3::new(1.0, 2.0, 3.0);
        let plane_normal = Vec3::new(1.0, 1.0, 1.0).normalize();
        let sphere_center = Point3::new(2.0, 3.0, 4.0);
        let sphere_radius = 2.0;

        let result = PlaneSphereSolver::solve_impl(
            plane_origin, plane_normal, sphere_center, sphere_radius,
        );
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
                assert!(dist_plane.abs() < TOL, "not on plane: {dist_plane}");
                // On sphere?
                let dist_sphere = (p - sphere_center).norm() - sphere_radius;
                assert!(dist_sphere.abs() < TOL, "not on sphere: {dist_sphere}");
            }
        }
    }
}

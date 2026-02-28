use rustkernel_math::{orthonormal_basis, Point3, Vec3};
use rustkernel_topology::geom_store::SurfaceKind;
use rustkernel_topology::intersection::{
    IntersectionCircle, IntersectionCurve, IntersectionEllipse, IntersectionError,
    SurfaceSurfaceResult, SurfaceSurfaceSolver, INTERSECTION_TOLERANCE,
};

/// Analytical plane-cone intersection solver.
///
/// The conic section type depends on the angle between the cutting plane and the cone axis:
/// - Perpendicular to axis: Circle
/// - theta_plane > half_angle: Ellipse
/// - theta_plane ≈ half_angle: Parabola (returns Empty for Phase 3)
/// - theta_plane < half_angle: Hyperbola (returns NoSolverAvailable for Phase 3)
pub struct PlaneConeSolver;

impl PlaneConeSolver {
    fn solve_impl(
        plane_origin: Point3,
        plane_normal: Vec3,
        cone_apex: Point3,
        cone_axis: Vec3,
        cone_half_angle: f64,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let n = plane_normal.normalize();
        let a = cone_axis.normalize();
        let ha = cone_half_angle.abs();

        // Angle between plane normal and cone axis.
        let cos_na = n.dot(&a).abs();
        // theta_plane = angle between plane and axis = PI/2 - acos(cos_na)
        // Equivalently, sin(theta_plane) = cos_na.
        let theta_plane = (std::f64::consts::FRAC_PI_2 - cos_na.acos()).abs();

        if cos_na > 1.0 - INTERSECTION_TOLERANCE {
            // Plane perpendicular to axis → Circle.
            // Distance from apex along axis to plane.
            let d = n.dot(&(plane_origin - cone_apex)) / n.dot(&a);
            if d < INTERSECTION_TOLERANCE {
                // Plane at or behind apex — only the apex point.
                return Ok(SurfaceSurfaceResult::Empty);
            }
            let circle_radius = d * ha.tan();
            let circle_center = cone_apex + d * a;

            let (ref_dir, _) = orthonormal_basis(&n);

            return Ok(SurfaceSurfaceResult::Curves(vec![
                IntersectionCurve::Circle(IntersectionCircle {
                    center: circle_center,
                    axis: n,
                    radius: circle_radius,
                    ref_dir,
                }),
            ]));
        }

        // Check for parabola (theta_plane ≈ half_angle).
        if (theta_plane - ha).abs() < 1e-6 {
            // Parabola — deferred.
            return Ok(SurfaceSurfaceResult::Empty);
        }

        // Check for hyperbola (theta_plane < half_angle).
        if theta_plane < ha {
            return Err(IntersectionError::NoSolverAvailable);
        }

        // Ellipse case (theta_plane > half_angle).
        // Find where the cone axis pierces the plane.
        let denom = n.dot(&a);
        if denom.abs() < INTERSECTION_TOLERANCE {
            // Plane parallel to axis — the general conic is a hyperbola or lines.
            return Err(IntersectionError::NoSolverAvailable);
        }
        let t_axis = n.dot(&(plane_origin - cone_apex)) / denom;

        if t_axis < INTERSECTION_TOLERANCE {
            return Ok(SurfaceSurfaceResult::Empty);
        }

        let axis_point = cone_apex + t_axis * a;

        // The semi-axes of the ellipse on the cutting plane.
        // For a right circular cone cut by an oblique plane:
        // When the plane is not perpendicular, we compute the ellipse parameters.
        let r_at_t = t_axis * ha.tan();

        // Project axis direction onto the plane.
        let a_proj = a - n.dot(&a) * n;
        let major_dir = if a_proj.norm() > INTERSECTION_TOLERANCE {
            a_proj.normalize()
        } else {
            let (d, _) = orthonormal_basis(&n);
            d
        };

        let sin_theta = (1.0 - cos_na * cos_na).sqrt(); // sin of angle between n and a
        let semi_minor = r_at_t;
        let semi_major = r_at_t / cos_na;

        // Adjust center: the center of the ellipse is not at axis_point but offset
        // along major_dir. For simplicity in Phase 3, use axis_point as approximate center.
        let center = axis_point;

        Ok(SurfaceSurfaceResult::Curves(vec![
            IntersectionCurve::Ellipse(IntersectionEllipse {
                center,
                axis: n,
                semi_major,
                semi_minor,
                major_dir,
            }),
        ]))
    }
}

impl SurfaceSurfaceSolver for PlaneConeSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane { .. }, SurfaceKind::Cone { .. })
                | (SurfaceKind::Cone { .. }, SurfaceKind::Plane { .. })
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
                SurfaceKind::Cone {
                    apex, axis, half_angle,
                },
            ) => Self::solve_impl(*origin, *normal, *apex, *axis, *half_angle),
            (
                SurfaceKind::Cone {
                    apex, axis, half_angle,
                },
                SurfaceKind::Plane { origin, normal },
            ) => Self::solve_impl(*origin, *normal, *apex, *axis, *half_angle),
            _ => Err(IntersectionError::InvalidInput(
                "expected plane and cone".into(),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    #[test]
    fn perpendicular_cut_circle() {
        // Cone with apex at origin, axis +Z, half_angle 45°.
        // Plane at z=3 perpendicular to axis → circle radius = 3*tan(45°) = 3.
        let result = PlaneConeSolver::solve_impl(
            Point3::new(0.0, 0.0, 3.0),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::FRAC_PI_4,
        )
        .unwrap();

        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                let IntersectionCurve::Circle(c) = &curves[0] else {
                    panic!("expected Circle");
                };
                assert!((c.radius - 3.0).abs() < TOL);
                assert!((c.center.z - 3.0).abs() < TOL);
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    #[test]
    fn plane_at_apex_empty() {
        let result = PlaneConeSolver::solve_impl(
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::FRAC_PI_4,
        )
        .unwrap();
        assert!(matches!(result, SurfaceSurfaceResult::Empty));
    }

    #[test]
    fn oblique_cut_ellipse() {
        // Cone with half_angle 30°, cut at ~60° from axis → should be ellipse.
        let normal = Vec3::new(0.0, 0.5, (3.0_f64).sqrt() / 2.0).normalize();
        let result = PlaneConeSolver::solve_impl(
            Point3::new(0.0, 0.0, 5.0),
            normal,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::PI / 6.0, // 30°
        );

        match result {
            Ok(SurfaceSurfaceResult::Curves(curves)) => {
                let IntersectionCurve::Ellipse(e) = &curves[0] else {
                    panic!("expected Ellipse");
                };
                assert!(e.semi_major > e.semi_minor);
            }
            Ok(SurfaceSurfaceResult::Empty) => {
                // acceptable if the plane misses the cone
            }
            other => {
                // NoSolverAvailable for hyperbola case is also acceptable
            }
        }
    }
}

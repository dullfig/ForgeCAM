use rustkernel_math::Point3;
use rustkernel_topology::geom_store::GeomAccess;
use tracing::{debug, warn};

/// Refine approximate intersection points onto both exact surfaces using
/// alternating projection.
///
/// For each input point:
/// 1. Project onto surface A → get closest point on A.
/// 2. Project onto surface B → get closest point on B.
/// 3. Repeat until the two projections converge.
/// 4. Return the midpoint.
///
/// Uses `GeomAccess::surface_inverse_uv` (Gauss-Newton for NURBS, analytical
/// for others) followed by `surface_eval` to snap points to exact geometry.
pub fn refine_intersection_points(
    geom: &dyn GeomAccess,
    surface_a: u32,
    surface_b: u32,
    points: &[Point3],
    max_iter: usize,
    tolerance: f64,
) -> Vec<Point3> {
    let mut refined = Vec::with_capacity(points.len());

    for (i, p) in points.iter().enumerate() {
        let mut current = *p;

        for _iter in 0..max_iter {
            // Project onto surface A
            let (u_a, v_a) = geom.surface_inverse_uv(surface_a, &current);
            let p_a = geom.surface_eval(surface_a, u_a, v_a);

            // Project onto surface B
            let (u_b, v_b) = geom.surface_inverse_uv(surface_b, &p_a);
            let p_b = geom.surface_eval(surface_b, u_b, v_b);

            let gap = (p_a - p_b).norm();
            current = Point3::from((p_a.coords + p_b.coords) * 0.5);

            if gap < tolerance {
                break;
            }
        }

        // Final check
        let (u_a, v_a) = geom.surface_inverse_uv(surface_a, &current);
        let p_a = geom.surface_eval(surface_a, u_a, v_a);
        let (u_b, v_b) = geom.surface_inverse_uv(surface_b, &current);
        let p_b = geom.surface_eval(surface_b, u_b, v_b);
        let final_gap = (p_a - p_b).norm();

        if final_gap > 1e-4 {
            warn!(
                point_idx = i,
                gap = final_gap,
                "alternating projection did not converge"
            );
        }

        refined.push(Point3::from((p_a.coords + p_b.coords) * 0.5));
    }

    debug!(
        count = refined.len(),
        "refine_intersection_points complete"
    );
    refined
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_math::{Point3, Vec3};
    use rustkernel_topology::geom_store::SurfaceKind;

    /// Simple geometry store with a plane and a sphere for testing refinement.
    struct TestGeom;

    impl GeomAccess for TestGeom {
        fn point(&self, _id: u32) -> Point3 { Point3::origin() }
        fn curve_eval(&self, _id: u32, _t: f64) -> Point3 { Point3::origin() }
        fn curve_kind(&self, _id: u32) -> rustkernel_topology::geom_store::CurveKind {
            rustkernel_topology::geom_store::CurveKind::Unknown
        }
        fn curve_tangent(&self, _id: u32, _t: f64) -> Vec3 { Vec3::x() }

        fn surface_eval(&self, surface_id: u32, u: f64, v: f64) -> Point3 {
            match surface_id {
                0 => {
                    // Plane z=0: point is (u, v, 0)
                    Point3::new(u, v, 0.0)
                }
                1 => {
                    // Sphere radius 2 at origin: (u, v) = (theta, phi)
                    let r = 2.0;
                    Point3::new(
                        r * v.cos() * u.cos(),
                        r * v.cos() * u.sin(),
                        r * v.sin(),
                    )
                }
                _ => Point3::origin(),
            }
        }

        fn surface_normal(&self, surface_id: u32, _u: f64, _v: f64) -> Vec3 {
            match surface_id {
                0 => Vec3::new(0.0, 0.0, 1.0),
                _ => Vec3::new(0.0, 0.0, 1.0),
            }
        }

        fn surface_kind(&self, surface_id: u32) -> SurfaceKind {
            match surface_id {
                0 => SurfaceKind::Plane {
                    origin: Point3::origin(),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                },
                1 => SurfaceKind::Sphere {
                    center: Point3::origin(),
                    radius: 2.0,
                },
                _ => SurfaceKind::Unknown,
            }
        }

        fn surface_inverse_uv(&self, surface_id: u32, point: &Point3) -> (f64, f64) {
            match surface_id {
                0 => {
                    // Plane z=0: (u,v) = (x, y)
                    (point.x, point.y)
                }
                1 => {
                    // Sphere: inverse mapping
                    let d = point.coords;
                    let r = d.norm();
                    if r < 1e-15 { return (0.0, 0.0); }
                    let phi = (d.z / r).asin();
                    let theta = d.y.atan2(d.x);
                    (theta, phi)
                }
                _ => (0.0, 0.0),
            }
        }
    }

    #[test]
    fn test_refine_plane_sphere() {
        let geom = TestGeom;

        // Plane z=0 intersects sphere r=2 in circle of radius 2 at z=0.
        // An approximate point near the intersection.
        let approx = vec![
            Point3::new(2.0, 0.0, 0.1),  // slightly off z=0
            Point3::new(0.0, 1.95, 0.05), // slightly off sphere
        ];

        let refined = refine_intersection_points(&geom, 0, 1, &approx, 10, 1e-10);
        assert_eq!(refined.len(), 2);

        // After refinement, points should be near z=0 and on the sphere
        for p in &refined {
            // Should be near z=0 (plane)
            assert!(p.z.abs() < 0.1, "refined point z={} not near plane", p.z);
            // Should be near sphere radius
            let r = p.coords.norm();
            assert!((r - 2.0).abs() < 0.1, "refined point r={} not near sphere", r);
        }
    }
}

use rustkernel_math::Point3;
use rustkernel_topology::geom_store::SurfaceKind;
use rustkernel_topology::intersection::{
    IntersectionCurve, IntersectionError, IntersectionLine, SurfaceSurfaceResult,
    SurfaceSurfaceSolver, INTERSECTION_TOLERANCE,
};

/// Analytical plane-plane intersection solver.
///
/// Given two planes (origin + normal each), computes the intersection line,
/// or detects parallel/coincident cases. Results are unbounded — trimming
/// to face boundaries is a separate operation.
pub struct PlanePlaneSolver;

impl SurfaceSurfaceSolver for PlanePlaneSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane { .. }, SurfaceKind::Plane { .. })
        )
    }

    fn solve(
        &self,
        a: &SurfaceKind,
        b: &SurfaceKind,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let (o1, n1) = match a {
            SurfaceKind::Plane { origin, normal } => (*origin, *normal),
            _ => return Err(IntersectionError::InvalidInput("expected plane".into())),
        };
        let (o2, n2) = match b {
            SurfaceKind::Plane { origin, normal } => (*origin, *normal),
            _ => return Err(IntersectionError::InvalidInput("expected plane".into())),
        };

        let d = n1.cross(&n2);
        let d_sq = d.dot(&d);

        if d_sq < INTERSECTION_TOLERANCE * INTERSECTION_TOLERANCE {
            // Parallel planes — check if coincident via signed distance.
            let signed_dist = n1.dot(&(o2 - o1));
            if signed_dist.abs() < INTERSECTION_TOLERANCE {
                return Ok(SurfaceSurfaceResult::Coincident);
            } else {
                return Ok(SurfaceSurfaceResult::Empty);
            }
        }

        // Not parallel — compute intersection line.
        // Point on line closest to world origin:
        //   p = (d1 * (n2 × d) + d2 * (d × n1)) / |d|²
        // where d1 = n1·o1, d2 = n2·o2
        let d1 = n1.dot(&o1.coords);
        let d2 = n2.dot(&o2.coords);

        let origin = Point3::from((d1 * n2.cross(&d) + d2 * d.cross(&n1)) / d_sq);
        let direction = d.normalize();

        Ok(SurfaceSurfaceResult::Curves(vec![
            IntersectionCurve::Line(IntersectionLine { origin, direction }),
        ]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{AnalyticalGeomStore, Plane};
    use crate::solvers::default_pipeline;
    use rustkernel_math::Vec3;
    use rustkernel_topology::geom_store::GeomAccess;

    const TOL: f64 = 1e-9;

    fn assert_near(a: f64, b: f64, msg: &str) {
        assert!((a - b).abs() < TOL, "{msg}: {a} vs {b}");
    }

    // 1. Orthogonal planes → line along correct axis
    #[test]
    fn orthogonal_xy_xz() {
        let solver = PlanePlaneSolver;
        let a = SurfaceKind::Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0), // XY plane
        };
        let b = SurfaceKind::Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 1.0, 0.0), // XZ plane
        };
        let result = solver.solve(&a, &b).unwrap();
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                let IntersectionCurve::Line(line) = &curves[0];
                // Intersection of XY and XZ planes is the X axis.
                // Direction should be ±X.
                assert_near(line.direction.y.abs(), 0.0, "dir.y");
                assert_near(line.direction.z.abs(), 0.0, "dir.z");
                assert_near(line.direction.x.abs(), 1.0, "dir.x");
                // Origin should be at world origin.
                assert_near(line.origin.x, 0.0, "origin.x");
                assert_near(line.origin.y, 0.0, "origin.y");
                assert_near(line.origin.z, 0.0, "origin.z");
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    // 2. Parallel disjoint → Empty
    #[test]
    fn parallel_disjoint() {
        let solver = PlanePlaneSolver;
        let a = SurfaceKind::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };
        let b = SurfaceKind::Plane {
            origin: Point3::new(0.0, 5.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = solver.solve(&a, &b).unwrap();
        assert!(matches!(result, SurfaceSurfaceResult::Empty));
    }

    // 3. Coincident planes → Coincident
    #[test]
    fn coincident() {
        let solver = PlanePlaneSolver;
        let a = SurfaceKind::Plane {
            origin: Point3::new(1.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };
        let b = SurfaceKind::Plane {
            origin: Point3::new(3.0, 0.0, 7.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = solver.solve(&a, &b).unwrap();
        assert!(matches!(result, SurfaceSurfaceResult::Coincident));
    }

    // 4. Anti-parallel coincident (opposite normals, same plane) → Coincident
    #[test]
    fn anti_parallel_coincident() {
        let solver = PlanePlaneSolver;
        let a = SurfaceKind::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let b = SurfaceKind::Plane {
            origin: Point3::new(5.0, 3.0, 0.0),
            normal: Vec3::new(0.0, 0.0, -1.0),
        };
        let result = solver.solve(&a, &b).unwrap();
        assert!(matches!(result, SurfaceSurfaceResult::Coincident));
    }

    // 5. Anti-parallel disjoint → Empty
    #[test]
    fn anti_parallel_disjoint() {
        let solver = PlanePlaneSolver;
        let a = SurfaceKind::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
        };
        let b = SurfaceKind::Plane {
            origin: Point3::new(2.0, 0.0, 0.0),
            normal: Vec3::new(-1.0, 0.0, 0.0),
        };
        let result = solver.solve(&a, &b).unwrap();
        assert!(matches!(result, SurfaceSurfaceResult::Empty));
    }

    // 6. Offset non-orthogonal planes → correct intersection point
    #[test]
    fn offset_non_orthogonal() {
        let solver = PlanePlaneSolver;
        // Plane A: y = 0 (XZ plane)
        let a = SurfaceKind::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };
        // Plane B: y = x (rotated 45° about Z, through origin)
        let n = Vec3::new(-1.0, 1.0, 0.0).normalize();
        let b = SurfaceKind::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: n,
        };
        let result = solver.solve(&a, &b).unwrap();
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                let IntersectionCurve::Line(line) = &curves[0];
                // Intersection should be along Z axis (both planes contain it).
                assert_near(line.direction.x.abs(), 0.0, "dir.x");
                assert_near(line.direction.y.abs(), 0.0, "dir.y");
                assert_near(line.direction.z.abs(), 1.0, "dir.z");
                // Origin should be at world origin (both pass through it).
                assert_near(line.origin.x, 0.0, "origin.x");
                assert_near(line.origin.y, 0.0, "origin.y");
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    // 7. Arbitrary planes → points on result line satisfy both plane equations
    #[test]
    fn arbitrary_planes_satisfy_both() {
        let solver = PlanePlaneSolver;
        let a = SurfaceKind::Plane {
            origin: Point3::new(1.0, 2.0, 3.0),
            normal: Vec3::new(1.0, 1.0, 0.0).normalize(),
        };
        let b = SurfaceKind::Plane {
            origin: Point3::new(-1.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 1.0, 1.0).normalize(),
        };
        let result = solver.solve(&a, &b).unwrap();
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                let IntersectionCurve::Line(line) = &curves[0];
                let n1 = Vec3::new(1.0, 1.0, 0.0).normalize();
                let n2 = Vec3::new(0.0, 1.0, 1.0).normalize();
                let o1 = Point3::new(1.0, 2.0, 3.0);
                let o2 = Point3::new(-1.0, 0.0, 5.0);
                // Check several points along the line.
                for t in &[-10.0, 0.0, 1.0, 42.0] {
                    let p = line.origin + *t * line.direction;
                    let dist_a = n1.dot(&(p - o1));
                    let dist_b = n2.dot(&(p - o2));
                    assert_near(dist_a, 0.0, &format!("plane A at t={t}"));
                    assert_near(dist_b, 0.0, &format!("plane B at t={t}"));
                }
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    // 8. Pipeline dispatches to PlanePlaneSolver via AnalyticalGeomStore
    #[test]
    fn pipeline_dispatches() {
        let mut geom = AnalyticalGeomStore::new();
        let s0 = geom.add_surface(Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        });
        let s1 = geom.add_surface(Plane {
            origin: Point3::origin(),
            normal: Vec3::new(1.0, 0.0, 0.0),
        });
        let pipeline = default_pipeline();
        let result = pipeline.solve(&geom, s0, s1).unwrap();
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                let IntersectionCurve::Line(line) = &curves[0];
                // XY ∩ YZ = Y axis
                assert_near(line.direction.y.abs(), 1.0, "dir.y");
                assert_near(line.direction.x.abs(), 0.0, "dir.x");
                assert_near(line.direction.z.abs(), 0.0, "dir.z");
            }
            other => panic!("expected Curves, got {other:?}"),
        }
    }

    // 9. Pipeline returns NoSolverAvailable for Unknown surface kinds
    #[test]
    fn pipeline_no_solver_for_unknown() {
        // Build a mock geom store that returns Unknown.
        struct UnknownStore;
        impl GeomAccess for UnknownStore {
            fn point(&self, _: u32) -> Point3 {
                unimplemented!()
            }
            fn curve_eval(&self, _: u32, _: f64) -> Point3 {
                unimplemented!()
            }
            fn surface_eval(&self, _: u32, _: f64, _: f64) -> Point3 {
                unimplemented!()
            }
            fn surface_normal(&self, _: u32, _: f64, _: f64) -> Vec3 {
                unimplemented!()
            }
            fn surface_kind(&self, _: u32) -> SurfaceKind {
                SurfaceKind::Unknown
            }
        }
        let pipeline = default_pipeline();
        let result = pipeline.solve(&UnknownStore, 0, 1);
        assert!(matches!(
            result,
            Err(rustkernel_topology::intersection::IntersectionError::NoSolverAvailable)
        ));
    }

    // 10. accepts() rejects non-plane kinds
    #[test]
    fn accepts_rejects_non_plane() {
        let solver = PlanePlaneSolver;
        let plane = SurfaceKind::Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 1.0, 0.0),
        };
        let unknown = SurfaceKind::Unknown;
        assert!(!solver.accepts(&plane, &unknown));
        assert!(!solver.accepts(&unknown, &plane));
        assert!(!solver.accepts(&unknown, &unknown));
        assert!(solver.accepts(&plane, &plane));
    }
}

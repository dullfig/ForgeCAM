use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
use rustkernel_topology::intersection::{
    IntersectionCurve, IntersectionError, IntersectionPolyline, SurfaceSurfaceResult,
    SurfaceSurfaceSolver,
};
use tracing::{debug, info_span, warn};

use crate::mesh_intersect::{chain_segments, mesh_mesh_intersect, mesh_plane_intersect};
use crate::refine::refine_intersection_points;

/// Default tessellation resolution for NURBS surfaces during SSI.
const NURBS_SSI_DIVS: usize = 16;

/// Tolerance for chaining segments.
const CHAIN_TOLERANCE: f64 = 0.05;

/// Solver for Plane × NURBS surface intersection.
///
/// Uses the fast path: slices the NURBS mesh against the plane, then refines
/// the resulting polyline points onto both exact surfaces.
pub struct PlaneNurbsSolver;

impl SurfaceSurfaceSolver for PlaneNurbsSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane { .. }, SurfaceKind::Nurbs)
                | (SurfaceKind::Nurbs, SurfaceKind::Plane { .. })
        )
    }

    fn solve(
        &self,
        _a: &SurfaceKind,
        _b: &SurfaceKind,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        // PlaneNurbsSolver requires geometry access; should not be called directly.
        Err(IntersectionError::InvalidInput(
            "PlaneNurbsSolver requires solve_with_geom".into(),
        ))
    }

    fn solve_with_geom(
        &self,
        geom: &dyn GeomAccess,
        surface_a: u32,
        surface_b: u32,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let _span = info_span!("plane_nurbs_solve", surface_a, surface_b).entered();

        let kind_a = geom.surface_kind(surface_a);
        let kind_b = geom.surface_kind(surface_b);

        // Determine which is the plane and which is the NURBS surface.
        let (plane_origin, plane_normal, nurbs_sid) = match (&kind_a, &kind_b) {
            (SurfaceKind::Plane { origin, normal }, SurfaceKind::Nurbs) => {
                (*origin, *normal, surface_b)
            }
            (SurfaceKind::Nurbs, SurfaceKind::Plane { origin, normal }) => {
                (*origin, *normal, surface_a)
            }
            _ => {
                return Err(IntersectionError::InvalidInput(
                    "expected one Plane and one Nurbs".into(),
                ));
            }
        };

        // Tessellate the NURBS surface.
        let mesh = match geom.tessellate_surface(nurbs_sid, NURBS_SSI_DIVS, NURBS_SSI_DIVS) {
            Some(m) => m,
            None => {
                warn!(nurbs_sid, "failed to tessellate NURBS surface for SSI");
                return Ok(SurfaceSurfaceResult::Empty);
            }
        };

        // Slice the mesh against the plane.
        let raw_segments = mesh_plane_intersect(&mesh, &plane_origin, &plane_normal);
        if raw_segments.is_empty() {
            debug!("no mesh-plane intersection segments");
            return Ok(SurfaceSurfaceResult::Empty);
        }

        // Chain segments into polylines.
        let chains = chain_segments(raw_segments, CHAIN_TOLERANCE);
        if chains.is_empty() {
            return Ok(SurfaceSurfaceResult::Empty);
        }

        // Refine each polyline onto both exact surfaces.
        let mut curves = Vec::new();
        for chain in chains {
            if chain.len() < 2 {
                continue;
            }
            let refined =
                refine_intersection_points(geom, surface_a, surface_b, &chain, 5, 1e-8);
            if refined.len() >= 2 {
                curves.push(IntersectionCurve::Polyline(IntersectionPolyline {
                    points: refined,
                }));
            }
        }

        if curves.is_empty() {
            Ok(SurfaceSurfaceResult::Empty)
        } else {
            debug!(curve_count = curves.len(), "plane-nurbs SSI complete");
            Ok(SurfaceSurfaceResult::Curves(curves))
        }
    }
}

/// Solver for NURBS × NURBS surface intersection.
///
/// Uses full mesh-mesh intersection: tessellates both surfaces, finds
/// triangle-triangle intersection segments, chains them, and refines.
pub struct NurbsNurbsSolver;

impl SurfaceSurfaceSolver for NurbsNurbsSolver {
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool {
        // Accept if at least one is NURBS and the other isn't handled by a more
        // specific solver. This catches NURBS×NURBS and also NURBS×Analytical
        // pairs that don't have dedicated solvers (e.g., NURBS×Cylinder).
        matches!(a, SurfaceKind::Nurbs) || matches!(b, SurfaceKind::Nurbs)
    }

    fn solve(
        &self,
        _a: &SurfaceKind,
        _b: &SurfaceKind,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        Err(IntersectionError::InvalidInput(
            "NurbsNurbsSolver requires solve_with_geom".into(),
        ))
    }

    fn solve_with_geom(
        &self,
        geom: &dyn GeomAccess,
        surface_a: u32,
        surface_b: u32,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let _span = info_span!("nurbs_nurbs_solve", surface_a, surface_b).entered();

        // Tessellate both surfaces.
        let mesh_a = tessellate_for_ssi(geom, surface_a)?;
        let mesh_b = tessellate_for_ssi(geom, surface_b)?;

        // Find all triangle-triangle intersection segments.
        let raw_segments = mesh_mesh_intersect(&mesh_a, &mesh_b);
        if raw_segments.is_empty() {
            debug!("no mesh-mesh intersection segments");
            return Ok(SurfaceSurfaceResult::Empty);
        }

        // Chain into polylines.
        let chains = chain_segments(raw_segments, CHAIN_TOLERANCE);
        if chains.is_empty() {
            return Ok(SurfaceSurfaceResult::Empty);
        }

        // Refine onto both exact surfaces.
        let mut curves = Vec::new();
        for chain in chains {
            if chain.len() < 2 {
                continue;
            }
            let refined =
                refine_intersection_points(geom, surface_a, surface_b, &chain, 5, 1e-8);
            if refined.len() >= 2 {
                curves.push(IntersectionCurve::Polyline(IntersectionPolyline {
                    points: refined,
                }));
            }
        }

        if curves.is_empty() {
            Ok(SurfaceSurfaceResult::Empty)
        } else {
            debug!(curve_count = curves.len(), "nurbs-nurbs SSI complete");
            Ok(SurfaceSurfaceResult::Curves(curves))
        }
    }
}

/// Tessellate a surface for SSI. For NURBS, uses `GeomAccess::tessellate_surface`.
/// For analytical surfaces, builds a mesh from parametric evaluation.
fn tessellate_for_ssi(
    geom: &dyn GeomAccess,
    surface_id: u32,
) -> Result<rustkernel_topology::mesh_cache::FaceMesh, IntersectionError> {
    // Try the GeomAccess tessellation first (works for NURBS).
    if let Some(mesh) = geom.tessellate_surface(surface_id, NURBS_SSI_DIVS, NURBS_SSI_DIVS) {
        return Ok(mesh);
    }

    // For analytical surfaces, build a mesh from parametric sampling.
    let ((u_min, u_max), (v_min, v_max)) = geom.surface_domain(surface_id);
    let divs = NURBS_SSI_DIVS;
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();

    for iv in 0..=divs {
        let v = u_min + (v_max - v_min) * (iv as f64 / divs as f64);
        for iu in 0..=divs {
            let u = u_min + (u_max - u_min) * (iu as f64 / divs as f64);
            positions.push(geom.surface_eval(surface_id, u, v));
            normals.push(geom.surface_normal(surface_id, u, v));
            uvs.push([u, v]);
        }
    }

    let row = divs + 1;
    for iv in 0..divs {
        for iu in 0..divs {
            let i00 = (iv * row + iu) as u32;
            let i10 = (iv * row + iu + 1) as u32;
            let i01 = ((iv + 1) * row + iu) as u32;
            let i11 = ((iv + 1) * row + iu + 1) as u32;
            indices.push(i00);
            indices.push(i10);
            indices.push(i11);
            indices.push(i00);
            indices.push(i11);
            indices.push(i01);
        }
    }

    Ok(rustkernel_topology::mesh_cache::FaceMesh {
        positions,
        normals,
        indices,
        uvs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_math::{Point3, Vec3};
    use rustkernel_topology::geom_store::CurveKind;
    use rustkernel_topology::mesh_cache::FaceMesh;

    /// Test geometry: surface 0 = plane z=0, surface 1 = plane z=x (tilted).
    struct TwoPlaneGeom;

    impl GeomAccess for TwoPlaneGeom {
        fn point(&self, _id: u32) -> Point3 { Point3::origin() }
        fn curve_eval(&self, _id: u32, _t: f64) -> Point3 { Point3::origin() }
        fn curve_kind(&self, _id: u32) -> CurveKind { CurveKind::Unknown }
        fn curve_tangent(&self, _id: u32, _t: f64) -> Vec3 { Vec3::x() }

        fn surface_eval(&self, sid: u32, u: f64, v: f64) -> Point3 {
            match sid {
                0 => Point3::new(u, v, 0.0), // z=0 plane
                1 => Point3::new(u, v, u),   // z=x tilted plane
                _ => Point3::origin(),
            }
        }

        fn surface_normal(&self, sid: u32, _u: f64, _v: f64) -> Vec3 {
            match sid {
                0 => Vec3::new(0.0, 0.0, 1.0),
                1 => Vec3::new(-1.0, 0.0, 1.0).normalize(),
                _ => Vec3::z(),
            }
        }

        fn surface_kind(&self, sid: u32) -> SurfaceKind {
            match sid {
                0 => SurfaceKind::Plane {
                    origin: Point3::origin(),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                },
                // Pretend the tilted plane is NURBS for testing
                1 => SurfaceKind::Nurbs,
                _ => SurfaceKind::Unknown,
            }
        }

        fn surface_domain(&self, _sid: u32) -> ((f64, f64), (f64, f64)) {
            ((-1.0, 1.0), (-1.0, 1.0))
        }

        fn tessellate_surface(&self, sid: u32, divs_u: usize, divs_v: usize) -> Option<FaceMesh> {
            if sid != 1 { return None; }
            // Build a mesh for the tilted plane
            let mut positions = Vec::new();
            let mut normals = Vec::new();
            let mut uvs = Vec::new();
            let mut indices = Vec::new();
            let n = Vec3::new(-1.0, 0.0, 1.0).normalize();

            for iv in 0..=divs_v {
                let v = -1.0 + 2.0 * (iv as f64 / divs_v as f64);
                for iu in 0..=divs_u {
                    let u = -1.0 + 2.0 * (iu as f64 / divs_u as f64);
                    positions.push(Point3::new(u, v, u));
                    normals.push(n);
                    uvs.push([u, v]);
                }
            }
            let row = divs_u + 1;
            for iv in 0..divs_v {
                for iu in 0..divs_u {
                    let i00 = (iv * row + iu) as u32;
                    let i10 = (iv * row + iu + 1) as u32;
                    let i01 = ((iv + 1) * row + iu) as u32;
                    let i11 = ((iv + 1) * row + iu + 1) as u32;
                    indices.push(i00);
                    indices.push(i10);
                    indices.push(i11);
                    indices.push(i00);
                    indices.push(i11);
                    indices.push(i01);
                }
            }

            Some(FaceMesh { positions, normals, indices, uvs })
        }

        fn surface_inverse_uv(&self, sid: u32, point: &Point3) -> (f64, f64) {
            match sid {
                0 => (point.x, point.y),
                1 => (point.x, point.y),
                _ => (0.0, 0.0),
            }
        }
    }

    #[test]
    fn test_plane_nurbs_solver() {
        let geom = TwoPlaneGeom;
        let solver = PlaneNurbsSolver;

        assert!(solver.accepts(
            &SurfaceKind::Plane {
                origin: Point3::origin(),
                normal: Vec3::z(),
            },
            &SurfaceKind::Nurbs,
        ));

        let result = solver.solve_with_geom(&geom, 0, 1).unwrap();
        match result {
            SurfaceSurfaceResult::Curves(curves) => {
                assert!(!curves.is_empty(), "should produce intersection curves");
                for curve in &curves {
                    match curve {
                        IntersectionCurve::Polyline(poly) => {
                            assert!(poly.points.len() >= 2, "polyline too short");
                            // All points should be near z=0 (on the plane)
                            for p in &poly.points {
                                assert!(
                                    p.z.abs() < 0.2,
                                    "point z={} too far from z=0 plane",
                                    p.z
                                );
                            }
                        }
                        _ => panic!("expected Polyline curve"),
                    }
                }
            }
            SurfaceSurfaceResult::Empty => {
                panic!("expected intersection, got empty");
            }
            SurfaceSurfaceResult::Coincident => {
                panic!("unexpected coincident result");
            }
        }
    }

    #[test]
    fn test_nurbs_nurbs_solver_accepts() {
        let solver = NurbsNurbsSolver;
        assert!(solver.accepts(&SurfaceKind::Nurbs, &SurfaceKind::Nurbs));
        assert!(solver.accepts(
            &SurfaceKind::Nurbs,
            &SurfaceKind::Cylinder {
                origin: Point3::origin(),
                axis: Vec3::z(),
                radius: 1.0,
            },
        ));
        assert!(!solver.accepts(
            &SurfaceKind::Plane {
                origin: Point3::origin(),
                normal: Vec3::z(),
            },
            &SurfaceKind::Sphere {
                center: Point3::origin(),
                radius: 1.0,
            },
        ));
    }
}

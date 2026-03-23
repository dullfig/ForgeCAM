// rustkernel-geom: geometry types + AnalyticalGeomStore + GeomAccess impl
//
// Type definitions (structs, enums, methods) are provided by forgecam-geometry.
// This crate re-exports them and provides the kernel-specific AnalyticalGeomStore
// and its GeomAccess implementation.

use rustkernel_math::{orthonormal_basis, Point3, Vec3};
use rustkernel_topology::geom_store::{inverse_map_from_kind, CurveKind, GeomAccess, SurfaceKind};
use rustkernel_topology::mesh_cache::FaceMesh;
use serde::{Serialize, Deserialize};

// ── Re-exports from forgecam-geometry ──

pub use forgecam_geometry::{
    // Curve structs
    LineSegment, CircleCurve, CircularArcCurve, EllipseCurve,
    // Surface structs
    Plane, CylinderSurface, SphereSurface, ConeSurface, TorusSurface,
    // NURBS wrapper
    FlippableNurbs,
    // Enums
    SurfaceDef, CurveDef,
    // Traits
    Curve, Surface,
};
pub use curvo::prelude::{NurbsCurve3D, NurbsSurface3D};

// ── Geometry Store ──

/// Concrete geometry store holding analytical geometry.
#[derive(Serialize, Deserialize)]
pub struct AnalyticalGeomStore {
    pub points: Vec<Point3>,
    pub curves: Vec<CurveDef>,
    pub surfaces: Vec<SurfaceDef>,
}

impl AnalyticalGeomStore {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            curves: Vec::new(),
            surfaces: Vec::new(),
        }
    }

    pub fn add_point(&mut self, p: Point3) -> u32 {
        let id = self.points.len() as u32;
        self.points.push(p);
        id
    }

    /// Add a generic curve definition.
    pub fn add_curve(&mut self, curve: CurveDef) -> u32 {
        let id = self.curves.len() as u32;
        self.curves.push(curve);
        id
    }

    /// Add a generic surface definition.
    pub fn add_surface(&mut self, surface: SurfaceDef) -> u32 {
        let id = self.surfaces.len() as u32;
        self.surfaces.push(surface);
        id
    }

    /// Convenience: add a line segment (wraps in CurveDef::LineSegment).
    pub fn add_line_segment(&mut self, seg: LineSegment) -> u32 {
        self.add_curve(CurveDef::LineSegment(seg))
    }

    /// Convenience: add a plane (wraps in SurfaceDef::Plane).
    pub fn add_plane(&mut self, plane: Plane) -> u32 {
        self.add_surface(SurfaceDef::Plane(plane))
    }

    /// Add a NURBS curve.
    pub fn add_nurbs_curve(&mut self, curve: NurbsCurve3D<f64>) -> u32 {
        self.add_curve(CurveDef::Nurbs(curve))
    }

    /// Add a NURBS surface.
    pub fn add_nurbs_surface(&mut self, surface: NurbsSurface3D<f64>) -> u32 {
        self.add_surface(SurfaceDef::Nurbs(FlippableNurbs { surface, flipped: false }))
    }
}

impl Default for AnalyticalGeomStore {
    fn default() -> Self {
        Self::new()
    }
}

// ── GeomAccess implementation ──
//
// Note: surface_eval / surface_normal use kernel-specific parameterization
// (e.g., sphere uses orthonormal_basis-derived frame) that differs from
// geometry's Surface trait. These are kept as-is for behavioral compatibility.
// Curve methods delegate to geometry's Curve trait (parameterizations match).

impl GeomAccess for AnalyticalGeomStore {
    fn point(&self, point_id: u32) -> Point3 {
        self.points[point_id as usize]
    }

    fn curve_eval(&self, curve_id: u32, t: f64) -> Point3 {
        self.curves[curve_id as usize].eval(t)
    }

    fn curve_kind(&self, curve_id: u32) -> CurveKind {
        match &self.curves[curve_id as usize] {
            CurveDef::LineSegment(seg) => CurveKind::LineSegment {
                start: seg.start,
                end: seg.end,
            },
            CurveDef::Circle(c) => CurveKind::Circle {
                center: c.center,
                axis: c.axis,
                radius: c.radius,
            },
            CurveDef::CircularArc(a) => CurveKind::CircularArc {
                center: a.center,
                axis: a.axis,
                radius: a.radius,
                ref_dir: a.ref_dir,
                start_angle: a.start_angle,
                end_angle: a.end_angle,
            },
            CurveDef::Ellipse(e) => CurveKind::Ellipse {
                center: e.center,
                axis: e.axis,
                semi_major: e.semi_major,
                semi_minor: e.semi_minor,
                major_dir: e.major_dir,
            },
            CurveDef::Nurbs(_) => CurveKind::Nurbs,
        }
    }

    fn curve_tangent(&self, curve_id: u32, t: f64) -> Vec3 {
        self.curves[curve_id as usize].tangent(t)
    }

    fn surface_eval(&self, surface_id: u32, u: f64, v: f64) -> Point3 {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Plane(plane) => {
                // For a plane, (u, v) are ignored (vertex positions are read directly).
                let _ = (u, v);
                plane.origin
            }
            SurfaceDef::Cylinder(cyl) => {
                let (ref_x, ref_y) = orthonormal_basis(&cyl.axis);
                cyl.origin
                    + v * cyl.axis
                    + cyl.radius.abs() * (u.cos() * ref_x + u.sin() * ref_y)
            }
            SurfaceDef::Sphere(sph) => {
                let (ref_x, ref_y) = orthonormal_basis(&Vec3::new(0.0, 0.0, 1.0));
                let axis = Vec3::new(0.0, 0.0, 1.0);
                let r = sph.radius.abs();
                sph.center
                    + r * (v.cos() * u.cos() * ref_x + v.cos() * u.sin() * ref_y + v.sin() * axis)
            }
            SurfaceDef::Cone(cone) => {
                let (ref_x, ref_y) = orthonormal_basis(&cone.axis);
                let r = v * cone.half_angle.abs().tan();
                cone.apex + v * cone.axis + r * (u.cos() * ref_x + u.sin() * ref_y)
            }
            SurfaceDef::Torus(tor) => {
                let (ref_x, ref_y) = orthonormal_basis(&tor.axis);
                let big_r = tor.major_radius;
                let small_r = tor.minor_radius.abs();
                tor.center
                    + (big_r + small_r * v.cos()) * (u.cos() * ref_x + u.sin() * ref_y)
                    + small_r * v.sin() * tor.axis
            }
            SurfaceDef::Nurbs(ns) => ns.surface.point_at(u, v),
        }
    }

    fn surface_normal(&self, surface_id: u32, u: f64, v: f64) -> Vec3 {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Plane(plane) => plane.normal,
            SurfaceDef::Cylinder(cyl) => {
                let (ref_x, ref_y) = orthonormal_basis(&cyl.axis);
                let n = u.cos() * ref_x + u.sin() * ref_y;
                if cyl.radius < 0.0 { -n } else { n }
            }
            SurfaceDef::Sphere(sph) => {
                let (ref_x, ref_y) = orthonormal_basis(&Vec3::new(0.0, 0.0, 1.0));
                let axis = Vec3::new(0.0, 0.0, 1.0);
                let n = v.cos() * u.cos() * ref_x + v.cos() * u.sin() * ref_y + v.sin() * axis;
                if sph.radius < 0.0 { -n } else { n }
            }
            SurfaceDef::Cone(cone) => {
                let (ref_x, ref_y) = orthonormal_basis(&cone.axis);
                let ha = cone.half_angle.abs();
                let n = ha.cos() * (u.cos() * ref_x + u.sin() * ref_y) - ha.sin() * cone.axis;
                if cone.half_angle < 0.0 { -n } else { n }
            }
            SurfaceDef::Torus(tor) => {
                let (ref_x, ref_y) = orthonormal_basis(&tor.axis);
                let n = v.cos() * (u.cos() * ref_x + u.sin() * ref_y) + v.sin() * tor.axis;
                if tor.minor_radius < 0.0 { -n } else { n }
            }
            SurfaceDef::Nurbs(ns) => {
                let n = ns.surface.normal_at(u, v);
                let len = n.norm();
                let normalized = if len > 1e-15 { n / len } else { Vec3::new(0.0, 0.0, 1.0) };
                if ns.flipped { -normalized } else { normalized }
            }
        }
    }

    fn surface_kind(&self, surface_id: u32) -> SurfaceKind {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Plane(plane) => SurfaceKind::Plane {
                origin: plane.origin,
                normal: plane.normal,
            },
            SurfaceDef::Cylinder(cyl) => SurfaceKind::Cylinder {
                origin: cyl.origin,
                axis: cyl.axis,
                radius: cyl.radius,
            },
            SurfaceDef::Sphere(sph) => SurfaceKind::Sphere {
                center: sph.center,
                radius: sph.radius,
            },
            SurfaceDef::Cone(cone) => SurfaceKind::Cone {
                apex: cone.apex,
                axis: cone.axis,
                half_angle: cone.half_angle,
            },
            SurfaceDef::Torus(tor) => SurfaceKind::Torus {
                center: tor.center,
                axis: tor.axis,
                major_radius: tor.major_radius,
                minor_radius: tor.minor_radius,
            },
            SurfaceDef::Nurbs(_) => SurfaceKind::Nurbs,
        }
    }

    fn curve_domain(&self, curve_id: u32) -> (f64, f64) {
        self.curves[curve_id as usize].domain()
    }

    fn surface_domain(&self, surface_id: u32) -> ((f64, f64), (f64, f64)) {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Nurbs(ns) => ns.surface.knots_domain(),
            _ => ((0.0, 1.0), (0.0, 1.0)),
        }
    }

    fn surface_inverse_uv(&self, surface_id: u32, point: &Point3) -> (f64, f64) {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Nurbs(ns) => {
                // Delegate to geometry's Surface::inverse_uv for NURBS
                use forgecam_geometry::Surface;
                ns.inverse_uv(point)
            }
            _ => {
                // Analytical surfaces use kernel's parameterization
                let kind = self.surface_kind(surface_id);
                inverse_map_from_kind(&kind, point)
            }
        }
    }

    fn tessellate_surface(&self, surface_id: u32, divs_u: usize, divs_v: usize) -> Option<FaceMesh> {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Nurbs(ns) => {
                let ((u_min, u_max), (v_min, v_max)) = ns.surface.knots_domain();
                let flip = ns.flipped;
                let mut positions = Vec::with_capacity((divs_u + 1) * (divs_v + 1));
                let mut normals = Vec::with_capacity((divs_u + 1) * (divs_v + 1));
                let mut uvs = Vec::with_capacity((divs_u + 1) * (divs_v + 1));
                for iv in 0..=divs_v {
                    let v = v_min + (v_max - v_min) * iv as f64 / divs_v as f64;
                    for iu in 0..=divs_u {
                        let u = u_min + (u_max - u_min) * iu as f64 / divs_u as f64;
                        positions.push(ns.surface.point_at(u, v));
                        let n = ns.surface.normal_at(u, v);
                        let len = n.norm();
                        let normalized = if len > 1e-15 { n / len } else { Vec3::new(0.0, 0.0, 1.0) };
                        normals.push(if flip { -normalized } else { normalized });
                        uvs.push([u, v]);
                    }
                }
                let mut indices = Vec::with_capacity(divs_u * divs_v * 6);
                let w = (divs_u + 1) as u32;
                for iv in 0..divs_v as u32 {
                    for iu in 0..divs_u as u32 {
                        let i00 = iv * w + iu;
                        let i10 = i00 + 1;
                        let i01 = i00 + w;
                        let i11 = i01 + 1;
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
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use curvo::prelude::Interpolation;

    // Helper: create a NURBS curve through 4 points
    fn make_test_curve() -> NurbsCurve3D<f64> {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
        ];
        NurbsCurve3D::<f64>::interpolate(&pts, 3).unwrap()
    }

    // Helper: create a NURBS surface by extruding a line along Z
    fn make_test_surface() -> NurbsSurface3D<f64> {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let curve = NurbsCurve3D::<f64>::interpolate(&pts, 1).unwrap();
        NurbsSurface3D::<f64>::extrude(&curve, &Vec3::new(0.0, 0.0, 2.0))
    }

    #[test]
    fn test_nurbs_curve_eval() {
        let mut geom = AnalyticalGeomStore::new();
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
        ];
        let curve = NurbsCurve3D::<f64>::interpolate(&pts, 3).unwrap();
        let (t0, t1) = curve.knots_domain();
        let cid = geom.add_nurbs_curve(curve);
        let p_start = geom.curve_eval(cid, t0);
        let p_end = geom.curve_eval(cid, t1);
        assert!((p_start - pts[0]).norm() < 1e-8, "Curve should start at first point");
        assert!((p_end - pts[3]).norm() < 1e-8, "Curve should end at last point");
    }

    #[test]
    fn test_nurbs_curve_tangent() {
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_test_curve();
        let (t0, t1) = curve.knots_domain();
        let cid = geom.add_nurbs_curve(curve);
        let t_mid = (t0 + t1) / 2.0;
        let tangent = geom.curve_tangent(cid, t_mid);
        assert!(tangent.norm() > 0.99, "Tangent should be unit-length");
    }

    #[test]
    fn test_nurbs_curve_domain() {
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_test_curve();
        let expected = curve.knots_domain();
        let cid = geom.add_nurbs_curve(curve);
        let domain = geom.curve_domain(cid);
        assert_eq!(domain, expected);
    }

    #[test]
    fn test_nurbs_surface_eval() {
        let mut geom = AnalyticalGeomStore::new();
        let surface = make_test_surface();
        let ((u0, u1), (v0, v1)) = surface.knots_domain();
        let sid = geom.add_nurbs_surface(surface);
        let corners = [
            geom.surface_eval(sid, u0, v0),
            geom.surface_eval(sid, u1, v0),
            geom.surface_eval(sid, u0, v1),
            geom.surface_eval(sid, u1, v1),
        ];
        let expected = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(1.0, 0.0, 2.0),
        ];
        for exp in &expected {
            let min_dist = corners.iter().map(|c| (c - exp).norm()).fold(f64::MAX, f64::min);
            assert!(min_dist < 1e-6, "Expected corner {:?} not found, corners: {:?}", exp, corners);
        }
    }

    #[test]
    fn test_nurbs_surface_normal() {
        let mut geom = AnalyticalGeomStore::new();
        let surface = make_test_surface();
        let ((u0, u1), (v0, v1)) = surface.knots_domain();
        let sid = geom.add_nurbs_surface(surface);
        let u_mid = (u0 + u1) / 2.0;
        let v_mid = (v0 + v1) / 2.0;
        let n = geom.surface_normal(sid, u_mid, v_mid);
        assert!(n.norm() > 0.99, "Normal should be unit-length");
    }

    #[test]
    fn test_nurbs_surface_tessellate() {
        let mut geom = AnalyticalGeomStore::new();
        let sid = geom.add_nurbs_surface(make_test_surface());
        let mesh = geom.tessellate_surface(sid, 4, 4);
        assert!(mesh.is_some(), "NURBS surface should tessellate");
        let mesh = mesh.unwrap();
        assert_eq!(mesh.positions.len(), 25);
        assert_eq!(mesh.triangle_count(), 32);
    }

    #[test]
    fn test_cylinder_surface_eval() {
        let mut geom = AnalyticalGeomStore::new();
        let sid = geom.add_surface(SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        }));

        let p = geom.surface_eval(sid, 0.0, 0.0);
        let dist_from_axis = (p.x * p.x + p.y * p.y).sqrt();
        assert!((dist_from_axis - 2.0).abs() < 1e-10, "Point should be at radius from axis");

        let p2 = geom.surface_eval(sid, std::f64::consts::FRAC_PI_2, 3.0);
        let dist2 = (p2.x * p2.x + p2.y * p2.y).sqrt();
        assert!((dist2 - 2.0).abs() < 1e-10);
        assert!((p2.z - 3.0).abs() < 1e-10);

        let n = geom.surface_normal(sid, 0.0, 0.0);
        assert!(n.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() < 1e-10, "Normal perpendicular to axis");
        assert!(n.norm() > 0.99);
    }

    #[test]
    fn test_nurbs_inverse_uv_extrude() {
        let mut geom = AnalyticalGeomStore::new();
        let surface = make_test_surface();
        let ((u0, u1), (v0, v1)) = surface.knots_domain();
        let sid = geom.add_nurbs_surface(surface);

        let test_params = [
            (u0, v0),
            (u1, v1),
            ((u0 + u1) / 2.0, (v0 + v1) / 2.0),
            (u0 + (u1 - u0) * 0.25, v0 + (v1 - v0) * 0.75),
        ];
        for &(u, v) in &test_params {
            let pt = geom.surface_eval(sid, u, v);
            let (u_inv, v_inv) = geom.surface_inverse_uv(sid, &pt);
            let pt_recon = geom.surface_eval(sid, u_inv, v_inv);
            let err = (pt - pt_recon).norm();
            assert!(
                err < 1e-4,
                "Inverse UV error too large: {} at ({}, {}) -> ({}, {})",
                err, u, v, u_inv, v_inv
            );
        }
    }

    #[test]
    fn test_nurbs_inverse_uv_analytical_unchanged() {
        let mut geom = AnalyticalGeomStore::new();
        let sid = geom.add_surface(SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        }));

        let p = geom.surface_eval(sid, 1.0, 3.0);
        let (u, v) = geom.surface_inverse_uv(sid, &p);
        assert!((u - 1.0).abs() < 1e-6, "u mismatch: {u}");
        assert!((v - 3.0).abs() < 1e-6, "v mismatch: {v}");
    }

    // ── Surface offset tests (now provided by geometry crate) ──

    #[test]
    fn test_offset_plane() {
        let plane = SurfaceDef::Plane(Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 2.0),
        });
        let offset = plane.offset(3.0).unwrap();
        match offset {
            SurfaceDef::Plane(p) => {
                assert!((p.origin.z - 8.0).abs() < 1e-10, "origin.z = {}", p.origin.z);
                assert!((p.normal - Vec3::new(0.0, 0.0, 2.0)).norm() < 1e-10, "normal preserved");
            }
            _ => panic!("Expected Plane"),
        }
    }

    #[test]
    fn test_offset_cylinder_outward() {
        let cyl = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        });
        let offset = cyl.offset(2.0).unwrap();
        match offset {
            SurfaceDef::Cylinder(c) => {
                assert!((c.radius - 7.0).abs() < 1e-10);
            }
            _ => panic!("Expected Cylinder"),
        }
    }

    #[test]
    fn test_offset_cylinder_degenerate() {
        let cyl = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        });
        assert!(cyl.offset(-3.0).is_none(), "Should degenerate when radius → 0");
    }

    #[test]
    fn test_offset_sphere() {
        let sph = SurfaceDef::Sphere(SphereSurface {
            center: Point3::new(1.0, 2.0, 3.0),
            radius: 4.0,
        });
        let offset = sph.offset(1.5).unwrap();
        match offset {
            SurfaceDef::Sphere(s) => {
                assert!((s.radius - 5.5).abs() < 1e-10);
                assert!((s.center - Point3::new(1.0, 2.0, 3.0)).norm() < 1e-10, "center unchanged");
            }
            _ => panic!("Expected Sphere"),
        }
    }

    #[test]
    fn test_offset_nurbs_returns_some() {
        let surface = make_test_surface();
        let nurbs = SurfaceDef::Nurbs(FlippableNurbs { surface, flipped: false });
        assert!(nurbs.offset(1.0).is_some(), "NURBS offset should succeed");
    }

    #[test]
    fn test_offset_preserves_sign_convention() {
        let cyl = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        });

        let mut a = cyl.offset(2.0).unwrap();
        a.flip_normal();

        let mut b = cyl.clone();
        b.flip_normal();
        let b = b.offset(-2.0).unwrap();

        match (&a, &b) {
            (SurfaceDef::Cylinder(ca), SurfaceDef::Cylinder(cb)) => {
                assert!((ca.radius - cb.radius).abs() < 1e-10,
                    "radii should match: {} vs {}", ca.radius, cb.radius);
            }
            _ => panic!("Expected Cylinder"),
        }
    }
}

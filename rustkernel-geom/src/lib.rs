use curvo::prelude::{NurbsCurve3D, NurbsSurface3D};
use rustkernel_math::{orthonormal_basis, Mat4, Point3, Vec3};
use rustkernel_topology::geom_store::{inverse_map_from_kind, CurveKind, GeomAccess, SurfaceKind};
use rustkernel_topology::mesh_cache::FaceMesh;
use tracing::warn;

// ── Primitive geometry types ──

/// A line segment between two points.
#[derive(Debug, Clone)]
pub struct LineSegment {
    pub start: Point3,
    pub end: Point3,
}

/// A plane defined by an origin point and a normal vector.
#[derive(Debug, Clone)]
pub struct Plane {
    pub origin: Point3,
    pub normal: Vec3,
}

/// A cylindrical surface defined by center of bottom circle, axis direction, and radius.
#[derive(Debug, Clone)]
pub struct CylinderSurface {
    pub origin: Point3,
    pub axis: Vec3,
    pub radius: f64,
}

/// A spherical surface defined by center and radius.
#[derive(Debug, Clone)]
pub struct SphereSurface {
    pub center: Point3,
    pub radius: f64,
}

/// A conical surface defined by apex, axis (from apex toward base), and half-angle.
#[derive(Debug, Clone)]
pub struct ConeSurface {
    pub apex: Point3,
    pub axis: Vec3,
    pub half_angle: f64,
}

/// A toroidal surface defined by center, axis, major radius, and minor radius.
#[derive(Debug, Clone)]
pub struct TorusSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// A full circle curve in 3D.
#[derive(Debug, Clone)]
pub struct CircleCurve {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub ref_dir: Vec3,
}

/// An arc of a circle in 3D.
#[derive(Debug, Clone)]
pub struct CircularArcCurve {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub ref_dir: Vec3,
    pub start_angle: f64,
    pub end_angle: f64,
}

/// An ellipse curve in 3D.
#[derive(Debug, Clone)]
pub struct EllipseCurve {
    pub center: Point3,
    pub axis: Vec3,
    pub semi_major: f64,
    pub semi_minor: f64,
    pub major_dir: Vec3,
}

// ── Storage enums ──

/// Polymorphic surface definition stored in the geometry store.
#[derive(Debug, Clone)]
pub enum SurfaceDef {
    Plane(Plane),
    Cylinder(CylinderSurface),
    Sphere(SphereSurface),
    Cone(ConeSurface),
    Torus(TorusSurface),
    Nurbs(NurbsSurface3D<f64>),
}

impl SurfaceDef {
    /// Flip the surface normal orientation.
    pub fn flip_normal(&mut self) {
        match self {
            SurfaceDef::Plane(p) => p.normal = -p.normal,
            SurfaceDef::Cylinder(c) => c.radius = -c.radius,
            SurfaceDef::Sphere(s) => s.radius = -s.radius,
            SurfaceDef::Cone(c) => c.half_angle = -c.half_angle,
            SurfaceDef::Torus(t) => t.minor_radius = -t.minor_radius,
            SurfaceDef::Nurbs(_) => {} // no-op: NURBS booleans not yet supported
        }
    }

    /// Translate the surface by a delta vector.
    pub fn translate(&mut self, delta: &Vec3) {
        match self {
            SurfaceDef::Plane(p) => p.origin += delta,
            SurfaceDef::Cylinder(c) => c.origin += delta,
            SurfaceDef::Sphere(s) => s.center += delta,
            SurfaceDef::Cone(c) => c.apex += delta,
            SurfaceDef::Torus(t) => t.center += delta,
            SurfaceDef::Nurbs(ns) => {
                use curvo::prelude::Transformable;
                let m = Mat4::new_translation(delta);
                ns.transform(&m);
            }
        }
    }
}

/// Polymorphic curve definition stored in the geometry store.
#[derive(Debug, Clone)]
pub enum CurveDef {
    LineSegment(LineSegment),
    Circle(CircleCurve),
    CircularArc(CircularArcCurve),
    Ellipse(EllipseCurve),
    Nurbs(NurbsCurve3D<f64>),
}

impl CurveDef {
    /// Translate the curve by a delta vector.
    pub fn translate(&mut self, delta: &Vec3) {
        match self {
            CurveDef::LineSegment(s) => {
                s.start += delta;
                s.end += delta;
            }
            CurveDef::Circle(c) => c.center += delta,
            CurveDef::CircularArc(a) => a.center += delta,
            CurveDef::Ellipse(e) => e.center += delta,
            CurveDef::Nurbs(nc) => {
                use curvo::prelude::Transformable;
                let m = Mat4::new_translation(delta);
                nc.transform(&m);
            }
        }
    }
}

// ── Geometry Store ──

/// Concrete geometry store holding analytical geometry.
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
        self.add_surface(SurfaceDef::Nurbs(surface))
    }
}

impl Default for AnalyticalGeomStore {
    fn default() -> Self {
        Self::new()
    }
}

/// Find the closest (u,v) parameters on a NURBS surface to a 3D point.
/// Uses coarse grid search followed by Gauss-Newton refinement.
fn nurbs_closest_uv(ns: &NurbsSurface3D<f64>, point: &Point3) -> (f64, f64) {
    let ((u_min, u_max), (v_min, v_max)) = ns.knots_domain();
    let u_span = u_max - u_min;
    let v_span = v_max - v_min;

    // Phase 1: Grid search (8×8)
    let grid_n = 8;
    let mut best_u = u_min;
    let mut best_v = v_min;
    let mut best_dist2 = f64::MAX;

    for iv in 0..=grid_n {
        let v = v_min + v_span * iv as f64 / grid_n as f64;
        for iu in 0..=grid_n {
            let u = u_min + u_span * iu as f64 / grid_n as f64;
            let s = ns.point_at(u, v);
            let d2 = (s - point).norm_squared();
            if d2 < best_dist2 {
                best_dist2 = d2;
                best_u = u;
                best_v = v;
            }
        }
    }

    // Phase 2: Gauss-Newton refinement (5 iterations)
    let h = 1e-6 * (u_span + v_span);
    for _ in 0..5 {
        let s = ns.point_at(best_u, best_v);
        let r = point - s;

        if r.norm_squared() < 1e-24 {
            break;
        }

        // Finite-difference tangent vectors
        let s_u = (ns.point_at((best_u + h).min(u_max), best_v)
            - ns.point_at((best_u - h).max(u_min), best_v))
            / (2.0 * h);
        let s_v = (ns.point_at(best_u, (best_v + h).min(v_max))
            - ns.point_at(best_u, (best_v - h).max(v_min)))
            / (2.0 * h);

        // 2×2 system: J^T J [du, dv] = J^T r
        let a00 = s_u.dot(&s_u);
        let a01 = s_u.dot(&s_v);
        let a11 = s_v.dot(&s_v);
        let b0 = r.dot(&s_u);
        let b1 = r.dot(&s_v);

        let det = a00 * a11 - a01 * a01;
        if det.abs() < 1e-30 {
            break;
        }

        let du = (a11 * b0 - a01 * b1) / det;
        let dv = (a00 * b1 - a01 * b0) / det;

        best_u = (best_u + du).clamp(u_min, u_max);
        best_v = (best_v + dv).clamp(v_min, v_max);
    }

    let final_err = (ns.point_at(best_u, best_v) - point).norm();
    if final_err > 1e-4 {
        warn!(error = final_err, u = best_u, v = best_v, "NURBS inverse UV: poor convergence");
    }

    (best_u, best_v)
}

impl GeomAccess for AnalyticalGeomStore {
    fn point(&self, point_id: u32) -> Point3 {
        self.points[point_id as usize]
    }

    fn curve_eval(&self, curve_id: u32, t: f64) -> Point3 {
        match &self.curves[curve_id as usize] {
            CurveDef::LineSegment(seg) => {
                Point3::from(seg.start.coords.lerp(&seg.end.coords, t))
            }
            CurveDef::Circle(c) => {
                let (ref_x, _) = (c.ref_dir, c.axis.cross(&c.ref_dir));
                let ref_y = c.axis.cross(&c.ref_dir);
                let angle = t * std::f64::consts::TAU;
                c.center + c.radius * (angle.cos() * ref_x + angle.sin() * ref_y)
            }
            CurveDef::CircularArc(a) => {
                let ref_y = a.axis.cross(&a.ref_dir);
                let angle = a.start_angle + t * (a.end_angle - a.start_angle);
                a.center + a.radius * (angle.cos() * a.ref_dir + angle.sin() * ref_y)
            }
            CurveDef::Ellipse(e) => {
                let minor_dir = e.axis.cross(&e.major_dir);
                let angle = t * std::f64::consts::TAU;
                e.center + e.semi_major * angle.cos() * e.major_dir
                    + e.semi_minor * angle.sin() * minor_dir
            }
            CurveDef::Nurbs(nc) => nc.point_at(t),
        }
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
        match &self.curves[curve_id as usize] {
            CurveDef::LineSegment(seg) => (seg.end - seg.start).normalize(),
            CurveDef::Circle(c) => {
                let ref_y = c.axis.cross(&c.ref_dir);
                let angle = t * std::f64::consts::TAU;
                (-angle.sin() * c.ref_dir + angle.cos() * ref_y).normalize()
            }
            CurveDef::CircularArc(a) => {
                let ref_y = a.axis.cross(&a.ref_dir);
                let angle = a.start_angle + t * (a.end_angle - a.start_angle);
                (-angle.sin() * a.ref_dir + angle.cos() * ref_y).normalize()
            }
            CurveDef::Ellipse(e) => {
                let minor_dir = e.axis.cross(&e.major_dir);
                let angle = t * std::f64::consts::TAU;
                (-e.semi_major * angle.sin() * e.major_dir
                    + e.semi_minor * angle.cos() * minor_dir)
                    .normalize()
            }
            CurveDef::Nurbs(nc) => nc.tangent_at(t).normalize(),
        }
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
            SurfaceDef::Nurbs(ns) => ns.point_at(u, v),
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
                let n = ns.normal_at(u, v);
                let len = n.norm();
                if len > 1e-15 { n / len } else { Vec3::new(0.0, 0.0, 1.0) }
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
        match &self.curves[curve_id as usize] {
            CurveDef::Nurbs(nc) => nc.knots_domain(),
            _ => (0.0, 1.0),
        }
    }

    fn surface_domain(&self, surface_id: u32) -> ((f64, f64), (f64, f64)) {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Nurbs(ns) => ns.knots_domain(),
            _ => ((0.0, 1.0), (0.0, 1.0)),
        }
    }

    fn surface_inverse_uv(&self, surface_id: u32, point: &Point3) -> (f64, f64) {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Nurbs(ns) => nurbs_closest_uv(ns, point),
            _ => {
                let kind = self.surface_kind(surface_id);
                inverse_map_from_kind(&kind, point)
            }
        }
    }

    fn tessellate_surface(&self, surface_id: u32, divs_u: usize, divs_v: usize) -> Option<FaceMesh> {
        match &self.surfaces[surface_id as usize] {
            SurfaceDef::Nurbs(ns) => {
                let ((u_min, u_max), (v_min, v_max)) = ns.knots_domain();
                let mut positions = Vec::with_capacity((divs_u + 1) * (divs_v + 1));
                let mut normals = Vec::with_capacity((divs_u + 1) * (divs_v + 1));
                let mut uvs = Vec::with_capacity((divs_u + 1) * (divs_v + 1));
                for iv in 0..=divs_v {
                    let v = v_min + (v_max - v_min) * iv as f64 / divs_v as f64;
                    for iu in 0..=divs_u {
                        let u = u_min + (u_max - u_min) * iu as f64 / divs_u as f64;
                        positions.push(ns.point_at(u, v));
                        let n = ns.normal_at(u, v);
                        let len = n.norm();
                        normals.push(if len > 1e-15 { n / len } else { Vec3::new(0.0, 0.0, 1.0) });
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
        // Curve should pass through start and end points
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
        // Collect all 4 corners — extrude may map u/v differently
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
        // Each expected corner should be close to some actual corner
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
        // 5x5 grid = 25 vertices, 4x4 = 16 quads = 32 triangles
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

        // At u=0, v=0: point should be at distance radius from axis
        let p = geom.surface_eval(sid, 0.0, 0.0);
        let dist_from_axis = (p.x * p.x + p.y * p.y).sqrt();
        assert!((dist_from_axis - 2.0).abs() < 1e-10, "Point should be at radius from axis");

        // At u=PI/2, v=3.0: should still be at radius, z=3
        let p2 = geom.surface_eval(sid, std::f64::consts::FRAC_PI_2, 3.0);
        let dist2 = (p2.x * p2.x + p2.y * p2.y).sqrt();
        assert!((dist2 - 2.0).abs() < 1e-10);
        assert!((p2.z - 3.0).abs() < 1e-10);

        // Normal should be perpendicular to axis and outward
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

        // Pick several points on the surface, verify inverse_uv recovers them
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
        // Verify that analytical surfaces still work through the new surface_inverse_uv path
        let mut geom = AnalyticalGeomStore::new();
        let sid = geom.add_surface(SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        }));

        let p = geom.surface_eval(sid, 1.0, 3.0);
        let (u, v) = geom.surface_inverse_uv(sid, &p);
        // For cylinder: u should be ~1.0, v should be ~3.0
        assert!((u - 1.0).abs() < 1e-6, "u mismatch: {u}");
        assert!((v - 3.0).abs() < 1e-6, "v mismatch: {v}");
    }
}

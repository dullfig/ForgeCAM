use rustkernel_math::{orthonormal_basis, Point3, Vec3};
use rustkernel_topology::geom_store::{CurveKind, GeomAccess, SurfaceKind};

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
}

impl Default for AnalyticalGeomStore {
    fn default() -> Self {
        Self::new()
    }
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
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}

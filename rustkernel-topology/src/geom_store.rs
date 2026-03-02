use crate::mesh_cache::FaceMesh;
use rustkernel_math::{Point3, Vec3};

/// Classification of a surface with its geometric parameters.
/// Solvers pattern-match on this to avoid a second round-trip to the geometry store.
#[derive(Debug, Clone)]
pub enum SurfaceKind {
    Plane { origin: Point3, normal: Vec3 },
    Cylinder { origin: Point3, axis: Vec3, radius: f64 },
    Sphere { center: Point3, radius: f64 },
    Cone { apex: Point3, axis: Vec3, half_angle: f64 },
    Torus { center: Point3, axis: Vec3, major_radius: f64, minor_radius: f64 },
    Nurbs,
    Unknown,
}

/// Classification of a curve with its geometric parameters.
#[derive(Debug, Clone)]
pub enum CurveKind {
    Line { origin: Point3, direction: Vec3 },
    LineSegment { start: Point3, end: Point3 },
    Circle { center: Point3, axis: Vec3, radius: f64 },
    Ellipse { center: Point3, axis: Vec3, semi_major: f64, semi_minor: f64, major_dir: Vec3 },
    CircularArc { center: Point3, axis: Vec3, radius: f64, ref_dir: Vec3, start_angle: f64, end_angle: f64 },
    Nurbs,
    Unknown,
}

/// Trait that topology queries to access geometry without depending on concrete geometry types.
/// Concrete implementations live in rustkernel-primitives.
pub trait GeomAccess {
    /// Get the 3D position of a point by its geometry index.
    fn point(&self, point_id: u32) -> Point3;

    /// Evaluate a curve at parameter t, returning position.
    fn curve_eval(&self, curve_id: u32, t: f64) -> Point3;

    /// Classify a curve by kind, returning its geometric parameters.
    fn curve_kind(&self, curve_id: u32) -> CurveKind;

    /// Get the tangent vector of a curve at parameter t.
    fn curve_tangent(&self, curve_id: u32, t: f64) -> Vec3;

    /// Evaluate a surface at parameters (u, v), returning position.
    fn surface_eval(&self, surface_id: u32, u: f64, v: f64) -> Point3;

    /// Get the surface normal at parameters (u, v).
    fn surface_normal(&self, surface_id: u32, u: f64, v: f64) -> Vec3;

    /// Classify a surface by kind, returning its geometric parameters.
    fn surface_kind(&self, surface_id: u32) -> SurfaceKind;

    /// Get the parameter domain of a curve. Default: (0.0, 1.0).
    fn curve_domain(&self, _curve_id: u32) -> (f64, f64) {
        (0.0, 1.0)
    }

    /// Get the parameter domain of a surface. Default: ((0.0, 1.0), (0.0, 1.0)).
    fn surface_domain(&self, _surface_id: u32) -> ((f64, f64), (f64, f64)) {
        ((0.0, 1.0), (0.0, 1.0))
    }

    /// Tessellate a surface into a mesh (used for NURBS; analytical surfaces return None).
    fn tessellate_surface(&self, _surface_id: u32, _divs_u: usize, _divs_v: usize) -> Option<FaceMesh> {
        None
    }
}

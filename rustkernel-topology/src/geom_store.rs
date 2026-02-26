use rustkernel_math::{Point3, Vec3};

/// Classification of a surface with its geometric parameters.
/// Solvers pattern-match on this to avoid a second round-trip to the geometry store.
#[derive(Debug, Clone)]
pub enum SurfaceKind {
    Plane { origin: Point3, normal: Vec3 },
    Unknown,
}

/// Trait that topology queries to access geometry without depending on concrete geometry types.
/// Concrete implementations live in rustkernel-primitives.
pub trait GeomAccess {
    /// Get the 3D position of a point by its geometry index.
    fn point(&self, point_id: u32) -> Point3;

    /// Evaluate a curve at parameter t, returning position.
    fn curve_eval(&self, curve_id: u32, t: f64) -> Point3;

    /// Evaluate a surface at parameters (u, v), returning position.
    fn surface_eval(&self, surface_id: u32, u: f64, v: f64) -> Point3;

    /// Get the surface normal at parameters (u, v).
    fn surface_normal(&self, surface_id: u32, u: f64, v: f64) -> Vec3;

    /// Classify a surface by kind, returning its geometric parameters.
    fn surface_kind(&self, surface_id: u32) -> SurfaceKind;
}

use rustkernel_math::{Point3, Vec3};

use crate::geom_store::{GeomAccess, SurfaceKind};

/// Tolerance for intersection computations (1e-10).
pub const INTERSECTION_TOLERANCE: f64 = 1e-10;

/// An infinite line in 3D space, defined by a point and a unit direction.
#[derive(Debug, Clone)]
pub struct IntersectionLine {
    pub origin: Point3,
    pub direction: Vec3,
}

/// An intersection curve resulting from surface-surface intersection.
/// Enum rather than trait object — the set of analytical curve types is small and compile-time known.
#[derive(Debug, Clone)]
pub enum IntersectionCurve {
    Line(IntersectionLine),
}

/// Result of intersecting two surfaces.
#[derive(Debug, Clone)]
pub enum SurfaceSurfaceResult {
    /// The surfaces do not intersect.
    Empty,
    /// The surfaces intersect in one or more curves.
    Curves(Vec<IntersectionCurve>),
    /// The surfaces are coincident (overlap everywhere).
    Coincident,
}

/// Errors that can occur during intersection.
#[derive(Debug, Clone)]
pub enum IntersectionError {
    /// No registered solver can handle the given surface pair.
    NoSolverAvailable,
    /// The solver received invalid input.
    InvalidInput(String),
}

impl std::fmt::Display for IntersectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IntersectionError::NoSolverAvailable => write!(f, "no solver available for this surface pair"),
            IntersectionError::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
        }
    }
}

impl std::error::Error for IntersectionError {}

/// Trait for a surface-surface intersection solver.
/// Solvers work with pure geometry data (`SurfaceKind`), no geom store needed.
pub trait SurfaceSurfaceSolver {
    /// Returns true if this solver can handle the given surface pair.
    fn accepts(&self, a: &SurfaceKind, b: &SurfaceKind) -> bool;

    /// Compute the intersection of two surfaces.
    fn solve(&self, a: &SurfaceKind, b: &SurfaceKind) -> Result<SurfaceSurfaceResult, IntersectionError>;
}

/// Chain-of-responsibility pipeline: register solvers in order of specificity,
/// first one that accepts handles it, failures fall through.
pub struct IntersectionPipeline {
    solvers: Vec<Box<dyn SurfaceSurfaceSolver>>,
}

impl IntersectionPipeline {
    pub fn new() -> Self {
        Self {
            solvers: Vec::new(),
        }
    }

    /// Register a solver. Solvers are tried in registration order.
    pub fn register(&mut self, solver: Box<dyn SurfaceSurfaceSolver>) {
        self.solvers.push(solver);
    }

    /// Intersect two surfaces identified by their geometry indices.
    /// Extracts `SurfaceKind` from the geom store, then dispatches to solvers.
    pub fn solve(
        &self,
        geom: &dyn GeomAccess,
        surface_a: u32,
        surface_b: u32,
    ) -> Result<SurfaceSurfaceResult, IntersectionError> {
        let a = geom.surface_kind(surface_a);
        let b = geom.surface_kind(surface_b);

        for solver in &self.solvers {
            if solver.accepts(&a, &b) {
                return solver.solve(&a, &b);
            }
        }

        Err(IntersectionError::NoSolverAvailable)
    }
}

impl Default for IntersectionPipeline {
    fn default() -> Self {
        Self::new()
    }
}

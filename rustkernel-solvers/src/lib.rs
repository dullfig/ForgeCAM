pub mod mesh_intersect;
pub mod nurbs_solver;
pub mod plane_cone;
pub mod plane_cylinder;
pub mod plane_plane;
pub mod plane_sphere;
pub mod plane_torus;
pub mod refine;

use rustkernel_topology::intersection::IntersectionPipeline;

use crate::nurbs_solver::{NurbsNurbsSolver, PlaneNurbsSolver};
use crate::plane_cone::PlaneConeSolver;
use crate::plane_cylinder::PlaneCylinderSolver;
use crate::plane_plane::PlanePlaneSolver;
use crate::plane_sphere::PlaneSphereSolver;
use crate::plane_torus::PlaneTorusSolver;

/// Create an intersection pipeline pre-loaded with all solvers.
///
/// Analytical solvers are registered first (fast paths), followed by
/// NURBS solvers as fallbacks for surface pairs involving NURBS geometry.
pub fn default_pipeline() -> IntersectionPipeline {
    let mut pipeline = IntersectionPipeline::new();

    // Analytical solvers (fast, exact)
    pipeline.register(Box::new(PlanePlaneSolver));
    pipeline.register(Box::new(PlaneSphereSolver));
    pipeline.register(Box::new(PlaneCylinderSolver));
    pipeline.register(Box::new(PlaneConeSolver));
    pipeline.register(Box::new(PlaneTorusSolver));

    // NURBS solvers (mesh-based, fallback)
    pipeline.register(Box::new(PlaneNurbsSolver));
    pipeline.register(Box::new(NurbsNurbsSolver));

    pipeline
}

pub mod plane_cone;
pub mod plane_cylinder;
pub mod plane_plane;
pub mod plane_sphere;
pub mod plane_torus;

use rustkernel_topology::intersection::IntersectionPipeline;

use crate::solvers::plane_cone::PlaneConeSolver;
use crate::solvers::plane_cylinder::PlaneCylinderSolver;
use crate::solvers::plane_plane::PlanePlaneSolver;
use crate::solvers::plane_sphere::PlaneSphereSolver;
use crate::solvers::plane_torus::PlaneTorusSolver;

/// Create an intersection pipeline pre-loaded with all analytical solvers.
pub fn default_pipeline() -> IntersectionPipeline {
    let mut pipeline = IntersectionPipeline::new();
    pipeline.register(Box::new(PlanePlaneSolver));
    pipeline.register(Box::new(PlaneSphereSolver));
    pipeline.register(Box::new(PlaneCylinderSolver));
    pipeline.register(Box::new(PlaneConeSolver));
    pipeline.register(Box::new(PlaneTorusSolver));
    pipeline
}

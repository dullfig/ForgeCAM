pub mod plane_plane;

use rustkernel_topology::intersection::IntersectionPipeline;

use crate::solvers::plane_plane::PlanePlaneSolver;

/// Create an intersection pipeline pre-loaded with all analytical solvers.
pub fn default_pipeline() -> IntersectionPipeline {
    let mut pipeline = IntersectionPipeline::new();
    pipeline.register(Box::new(PlanePlaneSolver));
    pipeline
}

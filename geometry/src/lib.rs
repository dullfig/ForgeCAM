// forgecam-geometry: pure computational geometry, no topology dependency
//
// See API.md at crate root for the full specification.

pub mod types;
pub mod utils;
pub mod curves;
pub mod curves2d;
pub mod contour2d;
pub mod polygon2d;
pub mod mesh;
pub mod surfaces;
pub mod intersection;
pub mod entity_ref;
pub mod engagement;

// Re-export core types at crate root for convenience
pub use types::{Mat4, Point2, Point3, Vec2, Vec3};
pub use utils::Aabb;
pub use curves::{
    CircleCurve, CircularArcCurve, Curve, CurveDef, CurveEnd, CurveKind, EllipseCurve,
    LineSegment, NurbsCurve3D,
};
pub use curves2d::{Arc2D, LineSegment2D};
pub use contour2d::{BooleanOp, Contour2D, Segment2D};
pub use surfaces::{
    ConeSurface, CylinderSurface, FlippableNurbs, NurbsSurface3D, Plane, SphereSurface, Surface,
    SurfaceDef, SurfaceKind, TorusSurface,
};
pub use mesh::TriMesh;
pub use polygon2d::Polygon2D;
pub use intersection::{SsiPipeline, SurfaceSurfaceResult, SsiSolver};
pub use entity_ref::{EntityKind, EntityRef};
pub use engagement::{
    engagement_angle, engagement_angle_simple, engagement_from_radial_depth,
    radial_depth_from_engagement, EngagedArc, EngagementResult,
};

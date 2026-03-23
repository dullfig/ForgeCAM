// Re-exports from forgecam-geometry::polygon2d

pub use forgecam_geometry::polygon2d::{Polygon2D, PointClassification, LineHit};

/// Compatibility alias: kernel code uses `Point2`, geometry crate uses `Point2D`.
pub type Point2 = forgecam_geometry::polygon2d::Point2D;

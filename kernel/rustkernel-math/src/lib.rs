// rustkernel-math: re-exports from forgecam-geometry
//
// This crate is a thin compatibility shim. All types and functions are
// provided by forgecam-geometry; we re-export them here so that existing
// `use rustkernel_math::*` throughout the kernel continues to work.

pub mod polygon2d;
pub mod tri_tri;

// Core types
pub use forgecam_geometry::{Point3, Vec3, Mat4};

// GPU helpers
pub use forgecam_geometry::types::{point3_to_f32, vec3_to_f32};

// nalgebra re-export (used by sketch solver, tests)
pub use forgecam_geometry::types::nalgebra;

// Utility functions
pub use forgecam_geometry::utils::orthonormal_basis;

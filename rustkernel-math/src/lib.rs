pub mod polygon2d;

pub use nalgebra;

pub type Point3 = nalgebra::Point3<f64>;
pub type Vec3 = nalgebra::Vector3<f64>;
pub type Mat4 = nalgebra::Matrix4<f64>;

/// Convert an f64 Point3 to an [f32; 3] array for GPU upload.
#[inline]
pub fn point3_to_f32(p: &Point3) -> [f32; 3] {
    [p.x as f32, p.y as f32, p.z as f32]
}

/// Convert an f64 Vec3 to an [f32; 3] array for GPU upload.
#[inline]
pub fn vec3_to_f32(v: &Vec3) -> [f32; 3] {
    [v.x as f32, v.y as f32, v.z as f32]
}

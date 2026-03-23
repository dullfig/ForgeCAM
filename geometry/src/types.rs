// Core geometric types — re-exports from nalgebra + GPU helpers

pub use nalgebra;

pub type Point2 = nalgebra::Point2<f64>;
pub type Point3 = nalgebra::Point3<f64>;
pub type Vec2 = nalgebra::Vector2<f64>;
pub type Vec3 = nalgebra::Vector3<f64>;
pub type Mat4 = nalgebra::Matrix4<f64>;

#[inline]
pub fn point3_to_f32(p: &Point3) -> [f32; 3] {
    [p.x as f32, p.y as f32, p.z as f32]
}

#[inline]
pub fn vec3_to_f32(v: &Vec3) -> [f32; 3] {
    [v.x as f32, v.y as f32, v.z as f32]
}

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

/// Given a unit vector `axis`, returns two perpendicular unit vectors `(u, v)`
/// forming a right-handed orthonormal frame `(u, v, axis)`.
pub fn orthonormal_basis(axis: &Vec3) -> (Vec3, Vec3) {
    let n = axis.normalize();
    let arbitrary = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = n.cross(&arbitrary).normalize();
    let v = n.cross(&u);
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orthonormal_basis_z_axis() {
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let (u, v) = orthonormal_basis(&axis);
        assert!((u.norm() - 1.0).abs() < 1e-12);
        assert!((v.norm() - 1.0).abs() < 1e-12);
        assert!(u.dot(&v).abs() < 1e-12);
        assert!(u.dot(&axis).abs() < 1e-12);
        assert!(v.dot(&axis).abs() < 1e-12);
    }

    #[test]
    fn test_orthonormal_basis_x_axis() {
        let axis = Vec3::new(1.0, 0.0, 0.0);
        let (u, v) = orthonormal_basis(&axis);
        assert!((u.norm() - 1.0).abs() < 1e-12);
        assert!((v.norm() - 1.0).abs() < 1e-12);
        assert!(u.dot(&v).abs() < 1e-12);
        assert!(u.dot(&axis).abs() < 1e-12);
        assert!(v.dot(&axis).abs() < 1e-12);
    }

    #[test]
    fn test_orthonormal_basis_arbitrary() {
        let axis = Vec3::new(1.0, 2.0, 3.0).normalize();
        let (u, v) = orthonormal_basis(&axis);
        assert!((u.norm() - 1.0).abs() < 1e-12);
        assert!((v.norm() - 1.0).abs() < 1e-12);
        assert!(u.dot(&v).abs() < 1e-12);
        assert!(u.dot(&axis).abs() < 1e-12);
        assert!(v.dot(&axis).abs() < 1e-12);
        // Right-handed: u x v = axis
        let cross = u.cross(&v);
        assert!((cross - axis).norm() < 1e-12);
    }
}

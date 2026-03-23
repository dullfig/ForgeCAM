// Geometric utilities: orthonormal_basis, rodrigues_rotate, AABB, transform helpers

use serde::{Deserialize, Serialize};

use crate::types::*;

// ---------------------------------------------------------------------------
// Orthonormal basis
// ---------------------------------------------------------------------------

/// Given a unit axis vector, return two perpendicular unit vectors (u, v)
/// forming a right-handed frame (u, v, axis).
pub fn orthonormal_basis(axis: &Vec3) -> (Vec3, Vec3) {
    let arbitrary = if axis.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = axis.cross(&arbitrary).normalize();
    let v = axis.cross(&u);
    (u, v)
}

// ---------------------------------------------------------------------------
// Rodrigues rotation
// ---------------------------------------------------------------------------

/// Rotate a vector around a unit axis by angle (radians) using Rodrigues' formula.
pub fn rodrigues_rotate(v: &Vec3, axis: &Vec3, angle: f64) -> Vec3 {
    let (sin_a, cos_a) = angle.sin_cos();
    *v * cos_a + axis.cross(v) * sin_a + axis * axis.dot(v) * (1.0 - cos_a)
}

/// Rotate a point around an axis through `axis_origin` by angle.
pub fn rodrigues_rotate_point(
    point: &Point3,
    axis_origin: &Point3,
    axis_dir: &Vec3,
    angle: f64,
) -> Point3 {
    let v = point - axis_origin;
    let rotated = rodrigues_rotate(&v, axis_dir, angle);
    axis_origin + rotated
}

// ---------------------------------------------------------------------------
// Arc generation
// ---------------------------------------------------------------------------

/// Generate n_segments+1 points along a circular arc via Rodrigues rotation.
pub fn generate_arc_points(
    center: &Point3,
    start: &Point3,
    axis: &Vec3,
    subtended_angle: f64,
    n_segments: usize,
) -> Vec<Point3> {
    let v = start - center;
    (0..=n_segments)
        .map(|i| {
            let angle = subtended_angle * i as f64 / n_segments as f64;
            let rotated = rodrigues_rotate(&v, axis, angle);
            center + rotated
        })
        .collect()
}

// ---------------------------------------------------------------------------
// AABB
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Aabb {
    pub min: Point3,
    pub max: Point3,
}

impl Aabb {
    /// Build from a set of points. Panics if empty.
    pub fn from_points(pts: &[Point3]) -> Self {
        assert!(!pts.is_empty(), "Aabb::from_points called with empty slice");
        let mut min = pts[0];
        let mut max = pts[0];
        for p in &pts[1..] {
            for i in 0..3 {
                if p[i] < min[i] {
                    min[i] = p[i];
                }
                if p[i] > max[i] {
                    max[i] = p[i];
                }
            }
        }
        Aabb { min, max }
    }

    /// Test overlap with another AABB (all three axes).
    pub fn overlaps(&self, other: &Aabb) -> bool {
        (0..3).all(|i| self.min[i] <= other.max[i] && self.max[i] >= other.min[i])
    }

    /// Expand box uniformly by `amount` in all directions.
    pub fn inflate(&mut self, amount: f64) {
        for i in 0..3 {
            self.min[i] -= amount;
            self.max[i] += amount;
        }
    }

    /// Merge two AABBs.
    pub fn union(&self, other: &Aabb) -> Aabb {
        let mut min = self.min;
        let mut max = self.max;
        for i in 0..3 {
            if other.min[i] < min[i] {
                min[i] = other.min[i];
            }
            if other.max[i] > max[i] {
                max[i] = other.max[i];
            }
        }
        Aabb { min, max }
    }

    /// Test if a point is inside (inclusive).
    pub fn contains_point(&self, p: &Point3) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] <= self.max[i])
    }
}

// ---------------------------------------------------------------------------
// Transform helpers
// ---------------------------------------------------------------------------

/// Apply the 3x3 rotation/scale sub-matrix of a Mat4 to a direction vector.
/// Normalizes result; returns original if degenerate.
pub fn transform_dir(m: &Mat4, v: &Vec3) -> Vec3 {
    let transformed = Vec3::new(
        m[(0, 0)] * v.x + m[(0, 1)] * v.y + m[(0, 2)] * v.z,
        m[(1, 0)] * v.x + m[(1, 1)] * v.y + m[(1, 2)] * v.z,
        m[(2, 0)] * v.x + m[(2, 1)] * v.y + m[(2, 2)] * v.z,
    );
    let norm = transformed.norm();
    if norm > 1e-15 {
        transformed / norm
    } else {
        *v
    }
}

/// Extract uniform scale factor from a Mat4 (column 0 norm).
pub fn scale_factor(m: &Mat4) -> f64 {
    Vec3::new(m[(0, 0)], m[(1, 0)], m[(2, 0)]).norm()
}

/// Transform a Point3 by a Mat4 (applies full affine transform).
pub fn transform_point(m: &Mat4, p: &Point3) -> Point3 {
    Point3::new(
        m[(0, 0)] * p.x + m[(0, 1)] * p.y + m[(0, 2)] * p.z + m[(0, 3)],
        m[(1, 0)] * p.x + m[(1, 1)] * p.y + m[(1, 2)] * p.z + m[(1, 3)],
        m[(2, 0)] * p.x + m[(2, 1)] * p.y + m[(2, 2)] * p.z + m[(2, 3)],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_orthonormal_basis_z_axis() {
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let (u, v) = orthonormal_basis(&axis);
        assert_relative_eq!(u.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(v.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(u.dot(&axis), 0.0, epsilon = 1e-12);
        assert_relative_eq!(v.dot(&axis), 0.0, epsilon = 1e-12);
        assert_relative_eq!(u.dot(&v), 0.0, epsilon = 1e-12);
        // Right-handed: u × v = axis
        let cross = u.cross(&v);
        assert_relative_eq!(cross.x, axis.x, epsilon = 1e-12);
        assert_relative_eq!(cross.y, axis.y, epsilon = 1e-12);
        assert_relative_eq!(cross.z, axis.z, epsilon = 1e-12);
    }

    #[test]
    fn test_orthonormal_basis_x_axis() {
        let axis = Vec3::new(1.0, 0.0, 0.0);
        let (u, v) = orthonormal_basis(&axis);
        assert_relative_eq!(u.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(v.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(u.dot(&axis), 0.0, epsilon = 1e-12);
        assert_relative_eq!(v.dot(&axis), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_rodrigues_90_degrees() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let rotated = rodrigues_rotate(&v, &axis, std::f64::consts::FRAC_PI_2);
        assert_relative_eq!(rotated.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(rotated.y, 1.0, epsilon = 1e-12);
        assert_relative_eq!(rotated.z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_rodrigues_360_returns_original() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let rotated = rodrigues_rotate(&v, &axis, std::f64::consts::TAU);
        assert_relative_eq!(rotated.x, v.x, epsilon = 1e-10);
        assert_relative_eq!(rotated.y, v.y, epsilon = 1e-10);
        assert_relative_eq!(rotated.z, v.z, epsilon = 1e-10);
    }

    #[test]
    fn test_rodrigues_rotate_point() {
        let point = Point3::new(1.0, 0.0, 0.0);
        let origin = Point3::new(0.0, 0.0, 0.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let rotated = rodrigues_rotate_point(&point, &origin, &axis, std::f64::consts::FRAC_PI_2);
        assert_relative_eq!(rotated.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(rotated.y, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_generate_arc_points_semicircle() {
        let center = Point3::new(0.0, 0.0, 0.0);
        let start = Point3::new(1.0, 0.0, 0.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let pts = generate_arc_points(&center, &start, &axis, std::f64::consts::PI, 4);
        assert_eq!(pts.len(), 5);
        // First point should be start
        assert_relative_eq!(pts[0].x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(pts[0].y, 0.0, epsilon = 1e-12);
        // Last point should be (-1, 0, 0)
        assert_relative_eq!(pts[4].x, -1.0, epsilon = 1e-12);
        assert_relative_eq!(pts[4].y, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_aabb_from_points() {
        let pts = vec![
            Point3::new(1.0, -2.0, 3.0),
            Point3::new(-1.0, 5.0, 0.0),
            Point3::new(0.0, 0.0, 10.0),
        ];
        let bb = Aabb::from_points(&pts);
        assert_eq!(bb.min, Point3::new(-1.0, -2.0, 0.0));
        assert_eq!(bb.max, Point3::new(1.0, 5.0, 10.0));
    }

    #[test]
    fn test_aabb_overlaps() {
        let a = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(2.0, 2.0, 2.0),
        };
        let b = Aabb {
            min: Point3::new(1.0, 1.0, 1.0),
            max: Point3::new(3.0, 3.0, 3.0),
        };
        assert!(a.overlaps(&b));

        let c = Aabb {
            min: Point3::new(5.0, 5.0, 5.0),
            max: Point3::new(6.0, 6.0, 6.0),
        };
        assert!(!a.overlaps(&c));
    }

    #[test]
    fn test_aabb_inflate() {
        let mut bb = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        bb.inflate(0.5);
        assert_eq!(bb.min, Point3::new(-0.5, -0.5, -0.5));
        assert_eq!(bb.max, Point3::new(1.5, 1.5, 1.5));
    }

    #[test]
    fn test_aabb_union() {
        let a = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Point3::new(-1.0, 2.0, 0.5),
            max: Point3::new(0.5, 3.0, 0.7),
        };
        let u = a.union(&b);
        assert_eq!(u.min, Point3::new(-1.0, 0.0, 0.0));
        assert_eq!(u.max, Point3::new(1.0, 3.0, 1.0));
    }

    #[test]
    fn test_aabb_contains_point() {
        let bb = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        assert!(bb.contains_point(&Point3::new(0.5, 0.5, 0.5)));
        assert!(bb.contains_point(&Point3::new(0.0, 0.0, 0.0))); // edge
        assert!(!bb.contains_point(&Point3::new(1.5, 0.5, 0.5)));
    }

    #[test]
    fn test_transform_dir_identity() {
        let m = Mat4::identity();
        let v = Vec3::new(1.0, 0.0, 0.0);
        let result = transform_dir(&m, &v);
        assert_relative_eq!(result.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(result.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(result.z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_scale_factor_identity() {
        let m = Mat4::identity();
        assert_relative_eq!(scale_factor(&m), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_scale_factor_scaled() {
        let m = Mat4::new_scaling(3.0);
        assert_relative_eq!(scale_factor(&m), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn test_transform_point_translation() {
        let m = Mat4::new_translation(&Vec3::new(1.0, 2.0, 3.0));
        let p = Point3::new(0.0, 0.0, 0.0);
        let result = transform_point(&m, &p);
        assert_relative_eq!(result.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(result.y, 2.0, epsilon = 1e-12);
        assert_relative_eq!(result.z, 3.0, epsilon = 1e-12);
    }
}

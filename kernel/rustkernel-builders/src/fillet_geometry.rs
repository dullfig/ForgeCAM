//! Fillet surface determination for the rolling-ball algorithm.
//!
//! Given two adjacent surface types, determines the correct SurfaceDef for
//! the fillet face. This module handles the geometry-only computation;
//! topology is handled by euler_fillet.rs.

use rustkernel_math::{Point3, Vec3};
use rustkernel_geom::{
    CylinderSurface, Plane, SurfaceDef, TorusSurface,
};

/// Determine the fillet surface for an edge between two faces.
///
/// Given the two adjacent face surfaces, the fillet center (on the rolling ball
/// axis at one endpoint), the edge direction, and the fillet radius, returns
/// the appropriate SurfaceDef.
///
/// Surface pair → fillet surface:
/// - Plane + Plane → CylinderSurface (axis = edge_dir, origin = center)
/// - Plane + Cylinder (perpendicular) → TorusSurface
/// - Other combinations → CylinderSurface approximation (future: NURBS)
pub fn determine_fillet_surface(
    surf_a: &SurfaceDef,
    surf_b: &SurfaceDef,
    center: &Point3,
    edge_dir: &Vec3,
    radius: f64,
) -> SurfaceDef {
    match (surf_a, surf_b) {
        // Plane-Plane → Cylinder (existing behavior)
        (SurfaceDef::Plane(_), SurfaceDef::Plane(_)) => {
            SurfaceDef::Cylinder(CylinderSurface {
                origin: *center,
                axis: *edge_dir,
                radius,
            })
        }

        // Plane-Cylinder (cap-to-side edge) → Torus
        (SurfaceDef::Plane(plane), SurfaceDef::Cylinder(cyl))
        | (SurfaceDef::Cylinder(cyl), SurfaceDef::Plane(plane)) => {
            make_plane_cylinder_torus(plane, cyl, radius, center, edge_dir)
        }

        // Fallback: use cylinder approximation (works well for short polygon edges)
        _ => {
            SurfaceDef::Cylinder(CylinderSurface {
                origin: *center,
                axis: *edge_dir,
                radius,
            })
        }
    }
}

/// Construct a TorusSurface for a Plane-Cylinder fillet.
///
/// The fillet torus has:
/// - center: on the cylinder axis, offset from the plane by -R along plane normal
/// - axis: cylinder axis
/// - major_radius: cylinder_radius - R (for convex edge)
/// - minor_radius: R
fn make_plane_cylinder_torus(
    plane: &Plane,
    cyl: &CylinderSurface,
    radius: f64,
    _center: &Point3,
    _edge_dir: &Vec3,
) -> SurfaceDef {
    let r_cyl = cyl.radius.abs();

    // Check if the plane normal is (anti-)parallel to the cylinder axis
    let dot = plane.normal.normalize().dot(&cyl.axis.normalize());
    let perpendicular = dot.abs() > 0.9; // close to perpendicular (cap edge)

    if !perpendicular || r_cyl <= radius {
        // Oblique cut or degenerate — fall back to cylinder approximation
        return SurfaceDef::Cylinder(CylinderSurface {
            origin: *_center,
            axis: *_edge_dir,
            radius,
        });
    }

    // Find where the cylinder axis meets the offset plane.
    // Offset plane: origin shifted by -R along normal (toward interior).
    let offset_origin = plane.origin - plane.normal.normalize() * radius;
    let n = plane.normal.normalize();
    let a = cyl.axis.normalize();

    // Parameter t where axis line hits offset plane:
    // (cyl.origin + t * a - offset_origin) · n = 0
    let denom = a.dot(&n);
    if denom.abs() < 1e-12 {
        // Axis parallel to plane — shouldn't happen for a cap edge
        return SurfaceDef::Cylinder(CylinderSurface {
            origin: *_center,
            axis: *_edge_dir,
            radius,
        });
    }
    let t = (offset_origin - cyl.origin).dot(&n) / denom;
    let torus_center = cyl.origin + t * a;

    let major_r = r_cyl - radius;

    SurfaceDef::Torus(TorusSurface {
        center: torus_center,
        axis: cyl.axis,
        major_radius: major_r,
        minor_radius: radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plane_plane_gives_cylinder() {
        let plane_a = SurfaceDef::Plane(Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        });
        let plane_b = SurfaceDef::Plane(Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 1.0, 0.0),
        });
        let center = Point3::new(1.0, 0.0, 0.0);
        let edge_dir = Vec3::new(1.0, 0.0, 0.0);
        let result = determine_fillet_surface(&plane_a, &plane_b, &center, &edge_dir, 0.5);
        assert!(matches!(result, SurfaceDef::Cylinder(_)));
    }

    #[test]
    fn test_plane_cylinder_gives_torus() {
        let plane = SurfaceDef::Plane(Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        });
        let cyl = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        });
        let center = Point3::new(2.5, 0.0, 4.5);
        let edge_dir = Vec3::new(0.0, 1.0, 0.0);
        let result = determine_fillet_surface(&plane, &cyl, &center, &edge_dir, 0.5);
        match &result {
            SurfaceDef::Torus(t) => {
                assert!((t.major_radius - 2.5).abs() < 1e-10, "major_r = {}", t.major_radius);
                assert!((t.minor_radius - 0.5).abs() < 1e-10, "minor_r = {}", t.minor_radius);
                assert!((t.center.z - 4.5).abs() < 1e-10, "center.z = {}", t.center.z);
            }
            other => panic!("Expected Torus, got {:?}", other),
        }
    }

    #[test]
    fn test_cylinder_plane_order_independent() {
        let plane = SurfaceDef::Plane(Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        });
        let cyl = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        });
        let center = Point3::new(2.5, 0.0, 4.5);
        let edge_dir = Vec3::new(0.0, 1.0, 0.0);

        let r1 = determine_fillet_surface(&plane, &cyl, &center, &edge_dir, 0.5);
        let r2 = determine_fillet_surface(&cyl, &plane, &center, &edge_dir, 0.5);
        assert!(matches!(r1, SurfaceDef::Torus(_)));
        assert!(matches!(r2, SurfaceDef::Torus(_)));
    }
}

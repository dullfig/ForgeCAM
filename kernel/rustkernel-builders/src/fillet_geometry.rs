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
/// - Cylinder + Cylinder (parallel axes) → TorusSurface
/// - Other combinations → CylinderSurface approximation (future: NURBS)
///
/// For concave edges, the fillet surface uses a negative radius (flipped normal),
/// following the project convention of negative radius = inward-facing.
pub fn determine_fillet_surface(
    surf_a: &SurfaceDef,
    surf_b: &SurfaceDef,
    center: &Point3,
    edge_dir: &Vec3,
    radius: f64,
    is_concave: bool,
) -> SurfaceDef {
    // For concave fillets, the surface faces inward (negative radius convention).
    let signed_radius = if is_concave { -radius } else { radius };

    match (surf_a, surf_b) {
        // Plane-Plane → Cylinder
        (SurfaceDef::Plane(_), SurfaceDef::Plane(_)) => {
            SurfaceDef::Cylinder(CylinderSurface {
                origin: *center,
                axis: *edge_dir,
                radius: signed_radius,
            })
        }

        // Plane-Cylinder (cap-to-side edge) → Torus
        (SurfaceDef::Plane(plane), SurfaceDef::Cylinder(cyl))
        | (SurfaceDef::Cylinder(cyl), SurfaceDef::Plane(plane)) => {
            make_plane_cylinder_torus(plane, cyl, radius, center, edge_dir, is_concave)
        }

        // Cylinder-Cylinder (fillet between two cylindrical faces)
        (SurfaceDef::Cylinder(cyl_a), SurfaceDef::Cylinder(cyl_b)) => {
            make_cylinder_cylinder_fillet(cyl_a, cyl_b, center, edge_dir, radius, is_concave)
        }

        // Fallback: cylinder approximation (works well for short polygon edges).
        // The vertex positions from arc computation are correct; the surface
        // is only used for tessellation UV mapping and downstream SSI.
        _ => {
            tracing::debug!(
                "fillet surface: using cylinder approximation for {:?} + {:?}",
                std::mem::discriminant(surf_a),
                std::mem::discriminant(surf_b)
            );
            SurfaceDef::Cylinder(CylinderSurface {
                origin: *center,
                axis: *edge_dir,
                radius: signed_radius,
            })
        }
    }
}

/// Construct a TorusSurface for a Plane-Cylinder fillet.
///
/// The fillet torus has:
/// - center: on the cylinder axis, offset from the plane by R along the axis
/// - axis: cylinder axis
/// - major_radius: cylinder_radius ± R (depending on convexity)
/// - minor_radius: R
fn make_plane_cylinder_torus(
    plane: &Plane,
    cyl: &CylinderSurface,
    radius: f64,
    _center: &Point3,
    _edge_dir: &Vec3,
    is_concave: bool,
) -> SurfaceDef {
    let r_cyl = cyl.radius.abs();

    // Check if the plane normal is (anti-)parallel to the cylinder axis
    let dot = plane.normal.normalize().dot(&cyl.axis.normalize());
    let perpendicular = dot.abs() > 0.9; // close to perpendicular (cap edge)

    if !perpendicular || r_cyl <= radius {
        // Oblique cut or degenerate — fall back to cylinder approximation
        let signed_r = if is_concave { -radius } else { radius };
        return SurfaceDef::Cylinder(CylinderSurface {
            origin: *_center,
            axis: *_edge_dir,
            radius: signed_r,
        });
    }

    // Find where the cylinder axis meets the offset plane.
    let offset_origin = plane.origin - plane.normal.normalize() * radius;
    let n = plane.normal.normalize();
    let a = cyl.axis.normalize();

    let denom = a.dot(&n);
    if denom.abs() < 1e-12 {
        let signed_r = if is_concave { -radius } else { radius };
        return SurfaceDef::Cylinder(CylinderSurface {
            origin: *_center,
            axis: *_edge_dir,
            radius: signed_r,
        });
    }
    let t = (offset_origin - cyl.origin).dot(&n) / denom;
    let torus_center = cyl.origin + t * a;

    // For convex: fillet torus sits inside the cylinder → major = R_cyl - R
    // For concave: fillet torus sits outside → major = R_cyl + R
    let major_r = if is_concave { r_cyl + radius } else { r_cyl - radius };

    let minor_r = if is_concave { -radius } else { radius };

    SurfaceDef::Torus(TorusSurface {
        center: torus_center,
        axis: cyl.axis,
        major_radius: major_r,
        minor_radius: minor_r,
    })
}

/// Construct a fillet surface for a Cylinder-Cylinder edge.
///
/// Common case: fillet between an existing fillet surface and another face
/// that was also filleted. If the cylinder axes are parallel (e.g., two
/// fillets on parallel edges of a box), the result is a torus.
/// Otherwise, falls back to cylinder approximation.
fn make_cylinder_cylinder_fillet(
    cyl_a: &CylinderSurface,
    cyl_b: &CylinderSurface,
    center: &Point3,
    edge_dir: &Vec3,
    radius: f64,
    is_concave: bool,
) -> SurfaceDef {
    let signed_radius = if is_concave { -radius } else { radius };

    let axis_a = cyl_a.axis.normalize();
    let axis_b = cyl_b.axis.normalize();

    // Check if axes are parallel (or anti-parallel).
    let dot = axis_a.dot(&axis_b).abs();
    if dot > 0.99 {
        // Parallel cylinders — the fillet is a torus section.
        // The torus axis is the shared cylinder axis direction.
        // Major radius = distance from the axis midpoint to the fillet center.
        let axis = if axis_a.dot(&axis_b) > 0.0 { axis_a } else { -axis_a };

        // Project center onto the line between the two cylinder axes
        // to find the torus center and major radius.
        let to_center = *center - cyl_a.origin;
        let axial = to_center.dot(&axis);
        let radial = to_center - axial * axis;
        let major_r = radial.norm();

        if major_r > 1e-12 {
            let torus_center = cyl_a.origin + axial * axis;
            return SurfaceDef::Torus(TorusSurface {
                center: torus_center,
                axis,
                major_radius: major_r,
                minor_radius: signed_radius,
            });
        }
    }

    // Non-parallel or degenerate — cylinder approximation.
    tracing::debug!(
        "fillet surface: cylinder-cylinder with non-parallel axes, using cylinder approx"
    );
    SurfaceDef::Cylinder(CylinderSurface {
        origin: *center,
        axis: *edge_dir,
        radius: signed_radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plane_plane_convex_gives_cylinder() {
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
        let result = determine_fillet_surface(&plane_a, &plane_b, &center, &edge_dir, 0.5, false);
        match &result {
            SurfaceDef::Cylinder(c) => {
                assert!(c.radius > 0.0, "convex fillet should have positive radius");
            }
            other => panic!("expected Cylinder, got {:?}", other),
        }
    }

    #[test]
    fn test_plane_plane_concave_gives_negative_radius() {
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
        let result = determine_fillet_surface(&plane_a, &plane_b, &center, &edge_dir, 0.5, true);
        match &result {
            SurfaceDef::Cylinder(c) => {
                assert!(c.radius < 0.0, "concave fillet should have negative radius, got {}", c.radius);
            }
            other => panic!("expected Cylinder, got {:?}", other),
        }
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
        let result = determine_fillet_surface(&plane, &cyl, &center, &edge_dir, 0.5, false);
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

        let r1 = determine_fillet_surface(&plane, &cyl, &center, &edge_dir, 0.5, false);
        let r2 = determine_fillet_surface(&cyl, &plane, &center, &edge_dir, 0.5, false);
        assert!(matches!(r1, SurfaceDef::Torus(_)));
        assert!(matches!(r2, SurfaceDef::Torus(_)));
    }

    #[test]
    fn test_cylinder_cylinder_parallel_gives_torus() {
        let cyl_a = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        });
        let cyl_b = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::new(1.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        });
        let center = Point3::new(0.5, 0.0, 1.0);
        let edge_dir = Vec3::new(0.0, 1.0, 0.0);
        let result = determine_fillet_surface(&cyl_a, &cyl_b, &center, &edge_dir, 0.25, false);
        assert!(matches!(result, SurfaceDef::Torus(_)), "parallel cylinders should give torus");
    }

    #[test]
    fn test_concave_plane_cylinder_torus() {
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
        let result = determine_fillet_surface(&plane, &cyl, &center, &edge_dir, 0.5, true);
        match &result {
            SurfaceDef::Torus(t) => {
                // Concave: major = R_cyl + R = 3.5, minor = -R = -0.5
                assert!((t.major_radius - 3.5).abs() < 1e-10, "major_r = {}", t.major_radius);
                assert!((t.minor_radius - (-0.5)).abs() < 1e-10, "minor_r = {}", t.minor_radius);
            }
            other => panic!("Expected Torus, got {:?}", other),
        }
    }
}

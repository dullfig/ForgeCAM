// All intersection algorithms: SSI, ray-surface, tri-tri, mesh-mesh, mesh-plane
// See API.md Section 4

use serde::{Deserialize, Serialize};

use crate::mesh::{tessellate_surface_grid, TriMesh};
use crate::surfaces::{
    ConeSurface, CylinderSurface, FlippableNurbs, Plane, SphereSurface, Surface, SurfaceDef,
    SurfaceKind, TorusSurface,
};
use crate::types::*;
use crate::utils::orthonormal_basis;

// ===========================================================================
// Result Types
// ===========================================================================

#[derive(Debug, Clone)]
pub enum SurfaceSurfaceResult {
    Curves(Vec<IntersectionCurve>),
    Coincident,
    Empty,
}

#[derive(Debug, Clone)]
pub enum IntersectionCurve {
    Line(IntersectionLine),
    Circle(IntersectionCircle),
    Ellipse(IntersectionEllipse),
    Polyline(IntersectionPolyline),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntersectionLine {
    pub origin: Point3,
    pub direction: Vec3, // unit direction
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntersectionCircle {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub ref_dir: Vec3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntersectionEllipse {
    pub center: Point3,
    pub axis: Vec3,
    pub semi_major: f64,
    pub semi_minor: f64,
    pub major_dir: Vec3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntersectionPolyline {
    pub points: Vec<Point3>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RayHit {
    pub t: f64,       // parameter along ray
    pub point: Point3, // hit point
    pub normal: Vec3,  // surface normal at hit
}

#[derive(Debug, Clone)]
pub struct TriTriResult {
    pub intersects: bool,
    pub segment: Option<(Point3, Point3)>,
}

impl TriTriResult {
    fn no_intersection() -> Self {
        TriTriResult {
            intersects: false,
            segment: None,
        }
    }
}

// ===========================================================================
// Internal helpers
// ===========================================================================

/// Snap near-zero values to exactly zero.
#[inline]
fn snap(v: f64) -> f64 {
    if v.abs() < 1e-12 {
        0.0
    } else {
        v
    }
}

/// Index of the component with largest absolute value.
fn best_axis_idx(v: &Vec3) -> usize {
    let mut best = 0;
    let mut best_val = v[0].abs();
    if v[1].abs() > best_val {
        best = 1;
        best_val = v[1].abs();
    }
    if v[2].abs() > best_val {
        best = 2;
    }
    best
}

/// Find crossing points of a triangle with a plane, given signed distances.
fn triangle_plane_crossings(
    verts: [&Point3; 3],
    dists: [f64; 3],
) -> Vec<Point3> {
    let mut pts = Vec::with_capacity(2);
    for i in 0..3 {
        let j = (i + 1) % 3;
        if dists[i] == 0.0 {
            pts.push(*verts[i]);
        } else if dists[i] * dists[j] < 0.0 {
            let t = dists[i] / (dists[i] - dists[j]);
            pts.push(Point3::from(
                verts[i].coords * (1.0 - t) + verts[j].coords * t,
            ));
        }
    }
    pts
}

// ===========================================================================
// Analytical Surface-Surface Intersection
// ===========================================================================

/// Plane-Plane intersection.
/// Returns Line, Coincident, or Empty.
pub fn intersect_plane_plane(a: &Plane, b: &Plane) -> SurfaceSurfaceResult {
    let dir = a.normal.cross(&b.normal);
    let dir_len = dir.norm();

    if dir_len < 1e-10 {
        // Parallel — check coincidence
        let d = (b.origin - a.origin).dot(&a.normal);
        return if d.abs() < 1e-10 {
            SurfaceSurfaceResult::Coincident
        } else {
            SurfaceSurfaceResult::Empty
        };
    }

    let dir = dir / dir_len;

    // Find a point on the intersection line.
    // Set the coordinate with largest |dir| to 0, solve the 2x2 system.
    let d1 = a.normal.dot(&a.origin.coords);
    let d2 = b.normal.dot(&b.origin.coords);

    let axis = best_axis_idx(&dir);
    let (j, k) = match axis {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };

    let det = a.normal[j] * b.normal[k] - a.normal[k] * b.normal[j];
    let mut origin = Point3::origin();
    origin[j] = (d1 * b.normal[k] - d2 * a.normal[k]) / det;
    origin[k] = (a.normal[j] * d2 - b.normal[j] * d1) / det;

    SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Line(IntersectionLine {
        origin,
        direction: dir,
    })])
}

/// Plane-Sphere intersection.
/// Returns Circle or Empty.
pub fn intersect_plane_sphere(p: &Plane, s: &SphereSurface) -> SurfaceSurfaceResult {
    let d = (s.center - p.origin).dot(&p.normal);
    let r = s.radius.abs();

    if d.abs() > r {
        return SurfaceSurfaceResult::Empty;
    }

    let circle_center = s.center - p.normal * d;
    let circle_radius = (r * r - d * d).sqrt();
    let (ref_dir, _) = orthonormal_basis(&p.normal);

    SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Circle(IntersectionCircle {
        center: circle_center,
        axis: p.normal,
        radius: circle_radius,
        ref_dir,
    })])
}

/// Plane-Cylinder intersection.
/// Returns Circle (perpendicular), Ellipse (oblique), or Empty (parallel).
pub fn intersect_plane_cylinder(p: &Plane, c: &CylinderSurface) -> SurfaceSurfaceResult {
    let dot = p.normal.dot(&c.axis);
    let abs_dot = dot.abs();

    if abs_dot < 1e-10 {
        // Plane parallel to cylinder axis — would produce lines, not circles/ellipses
        return SurfaceSurfaceResult::Empty;
    }

    // Intersection of plane with cylinder axis
    let t = (p.origin - c.origin).dot(&p.normal) / dot;
    let center = Point3::from(c.origin.coords + c.axis * t);
    let r = c.radius.abs();

    if abs_dot > 1.0 - 1e-10 {
        // Perpendicular → circle
        let (ref_dir, _) = orthonormal_basis(&p.normal);
        SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Circle(IntersectionCircle {
            center,
            axis: p.normal,
            radius: r,
            ref_dir,
        })])
    } else {
        // Oblique → ellipse
        let semi_minor = r;
        let semi_major = r / abs_dot;

        // Major direction: projection of cylinder axis onto cutting plane
        let axis_proj = c.axis - dot * p.normal;
        let major_dir = axis_proj.normalize();

        SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Ellipse(IntersectionEllipse {
            center,
            axis: p.normal,
            semi_major,
            semi_minor,
            major_dir,
        })])
    }
}

/// Plane-Cone intersection.
/// Returns Circle (perpendicular) or Empty.
pub fn intersect_plane_cone(p: &Plane, c: &ConeSurface) -> SurfaceSurfaceResult {
    let dot = p.normal.dot(&c.axis);
    let abs_dot = dot.abs();

    if abs_dot < 1e-10 {
        return SurfaceSurfaceResult::Empty;
    }

    // Intersection of plane with cone axis
    let t = (p.origin - c.apex).dot(&p.normal) / dot;

    if t < 0.0 {
        // Behind apex
        return SurfaceSurfaceResult::Empty;
    }

    let center = Point3::from(c.apex.coords + c.axis * t);
    let alpha = c.half_angle.abs();

    if abs_dot > 1.0 - 1e-10 {
        // Perpendicular to axis → circle
        let circle_radius = t * alpha.tan();
        let (ref_dir, _) = orthonormal_basis(&p.normal);
        SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Circle(IntersectionCircle {
            center,
            axis: p.normal,
            radius: circle_radius,
            ref_dir,
        })])
    } else {
        // Oblique — check if it's an ellipse (not parabola/hyperbola)
        let cos_alpha = alpha.cos();
        let sin_alpha = alpha.sin();
        let (u_basis, v_basis) = orthonormal_basis(&c.axis);
        let q = sin_alpha * u_basis.dot(&p.normal);
        let r = sin_alpha * v_basis.dot(&p.normal);
        let p_coeff = cos_alpha * dot.abs();
        let amplitude = (q * q + r * r).sqrt();

        if p_coeff <= amplitude {
            // Parabola or hyperbola — not a closed curve
            return SurfaceSurfaceResult::Empty;
        }

        // Ellipse case: sample parametrically and return as polyline
        // (exact ellipse parameters for oblique cone cuts are complex)
        let n_samples = 64;
        let mut points = Vec::with_capacity(n_samples + 1);
        let d_val = (p.origin - c.apex).dot(&p.normal);

        for i in 0..=n_samples {
            let u = std::f64::consts::TAU * i as f64 / n_samples as f64;
            let (sin_u, cos_u) = u.sin_cos();
            let denom = p_coeff + q * cos_u + r * sin_u;
            if denom.abs() < 1e-15 {
                continue;
            }
            let v = d_val / denom;
            if v < 0.0 {
                continue;
            }
            let radial = cos_u * u_basis + sin_u * v_basis;
            let pt = c.apex + c.axis * v * cos_alpha + radial * v * sin_alpha;
            points.push(pt);
        }

        if points.len() < 3 {
            return SurfaceSurfaceResult::Empty;
        }

        SurfaceSurfaceResult::Curves(vec![IntersectionCurve::Polyline(IntersectionPolyline {
            points,
        })])
    }
}

/// Plane-Torus intersection.
/// For axis-perpendicular planes: returns pair of Circles or Empty.
/// For oblique planes: returns Empty (handled by mesh-based fallback).
pub fn intersect_plane_torus(p: &Plane, t: &TorusSurface) -> SurfaceSurfaceResult {
    let dot = p.normal.dot(&t.axis);
    let abs_dot = dot.abs();

    if abs_dot < 1.0 - 1e-10 {
        // Oblique — not handled analytically
        return SurfaceSurfaceResult::Empty;
    }

    // Plane is (approximately) perpendicular to torus axis
    let h = (p.origin - t.center).dot(&t.axis);
    let r_minor = t.minor_radius.abs();

    if h.abs() > r_minor {
        return SurfaceSurfaceResult::Empty;
    }

    let cross_r = (r_minor * r_minor - h * h).sqrt();
    let inner_radius = t.major_radius - cross_r;
    let outer_radius = t.major_radius + cross_r;

    let circle_center = Point3::from(t.center.coords + t.axis * h);
    let (ref_dir, _) = orthonormal_basis(&p.normal);

    let mut curves = Vec::new();

    if outer_radius > 0.0 {
        curves.push(IntersectionCurve::Circle(IntersectionCircle {
            center: circle_center,
            axis: p.normal,
            radius: outer_radius,
            ref_dir,
        }));
    }

    if inner_radius > 0.0 {
        curves.push(IntersectionCurve::Circle(IntersectionCircle {
            center: circle_center,
            axis: p.normal,
            radius: inner_radius,
            ref_dir,
        }));
    }

    if curves.is_empty() {
        SurfaceSurfaceResult::Empty
    } else {
        SurfaceSurfaceResult::Curves(curves)
    }
}

// ===========================================================================
// Ray-Surface Intersection
// ===========================================================================

/// Ray-Plane intersection.
pub fn ray_plane(ray_origin: &Point3, ray_dir: &Vec3, plane: &Plane) -> Option<RayHit> {
    let denom = ray_dir.dot(&plane.normal);
    if denom.abs() < 1e-15 {
        return None; // parallel
    }

    let t = (plane.origin - ray_origin).dot(&plane.normal) / denom;
    if t < 1e-10 {
        return None;
    }

    Some(RayHit {
        t,
        point: Point3::from(ray_origin.coords + ray_dir * t),
        normal: plane.normal,
    })
}

/// Ray-Sphere intersection. Returns nearest positive hit.
pub fn ray_sphere(
    ray_origin: &Point3,
    ray_dir: &Vec3,
    sphere: &SphereSurface,
) -> Option<RayHit> {
    let oc = ray_origin - sphere.center;
    let r = sphere.radius.abs();

    let a = ray_dir.dot(ray_dir);
    let b = 2.0 * oc.dot(ray_dir);
    let c = oc.dot(&oc) - r * r;

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 || a.abs() < 1e-15 {
        return None;
    }

    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);

    let t = if t1 > 1e-10 {
        t1
    } else if t2 > 1e-10 {
        t2
    } else {
        return None;
    };

    let point = Point3::from(ray_origin.coords + ray_dir * t);
    let n = (point - sphere.center).normalize();
    let normal = if sphere.radius >= 0.0 { n } else { -n };

    Some(RayHit { t, point, normal })
}

/// Ray-Cylinder intersection (infinite cylinder).
pub fn ray_cylinder(
    ray_origin: &Point3,
    ray_dir: &Vec3,
    cyl: &CylinderSurface,
) -> Option<RayHit> {
    let oc = ray_origin - cyl.origin;
    let r = cyl.radius.abs();

    // Project perpendicular to axis
    let d = ray_dir - cyl.axis * ray_dir.dot(&cyl.axis);
    let oc_perp = oc - cyl.axis * oc.dot(&cyl.axis);

    let a = d.dot(&d);
    let b = 2.0 * d.dot(&oc_perp);
    let c = oc_perp.dot(&oc_perp) - r * r;

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 || a.abs() < 1e-15 {
        return None;
    }

    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);

    let t = if t1 > 1e-10 {
        t1
    } else if t2 > 1e-10 {
        t2
    } else {
        return None;
    };

    let point = Point3::from(ray_origin.coords + ray_dir * t);

    // Normal: radial direction from axis
    let dp = point - cyl.origin;
    let radial = (dp - cyl.axis * dp.dot(&cyl.axis)).normalize();
    let normal = if cyl.radius >= 0.0 { radial } else { -radial };

    Some(RayHit { t, point, normal })
}

/// Ray-Cone intersection (infinite cone from apex).
pub fn ray_cone(
    ray_origin: &Point3,
    ray_dir: &Vec3,
    cone: &ConeSurface,
) -> Option<RayHit> {
    let co = ray_origin - cone.apex;
    let alpha = cone.half_angle.abs();
    let cos2 = alpha.cos().powi(2);

    let d_dot_a = ray_dir.dot(&cone.axis);
    let co_dot_a = co.dot(&cone.axis);

    // Quadratic: (co + t*d)·a)² = cos²α * |co + t*d|²
    let a = d_dot_a * d_dot_a - cos2 * ray_dir.dot(ray_dir);
    let b = 2.0 * (d_dot_a * co_dot_a - cos2 * co.dot(ray_dir));
    let c = co_dot_a * co_dot_a - cos2 * co.dot(&co);

    let solve_hits = |a: f64, b: f64, c: f64| -> Option<RayHit> {
        if a.abs() < 1e-15 {
            // Linear
            if b.abs() < 1e-15 {
                return None;
            }
            let t = -c / b;
            if t < 1e-10 {
                return None;
            }
            let point = Point3::from(ray_origin.coords + ray_dir * t);
            let v = (point - cone.apex).dot(&cone.axis);
            if v < 0.0 {
                return None;
            }
            return Some(RayHit {
                t,
                point,
                normal: cone_normal(cone, &point),
            });
        }

        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }

        let sqrt_disc = disc.sqrt();
        for t in [(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)] {
            if t < 1e-10 {
                continue;
            }
            let point = Point3::from(ray_origin.coords + ray_dir * t);
            let v = (point - cone.apex).dot(&cone.axis);
            if v < 0.0 {
                continue; // wrong half of double cone
            }
            return Some(RayHit {
                t,
                point,
                normal: cone_normal(cone, &point),
            });
        }
        None
    };

    solve_hits(a, b, c)
}

fn cone_normal(cone: &ConeSurface, point: &Point3) -> Vec3 {
    let d = point - cone.apex;
    let v = d.dot(&cone.axis);
    let lateral = d - v * cone.axis;
    let alpha = cone.half_angle.abs();
    let ln = lateral.norm();
    if ln < 1e-15 {
        return cone.axis;
    }
    let radial = lateral / ln;
    let outward = alpha.cos() * radial - alpha.sin() * cone.axis;
    let n = outward.normalize();
    if cone.half_angle >= 0.0 {
        n
    } else {
        -n
    }
}

/// Ray vs any SurfaceDef — dispatches to appropriate function.
pub fn ray_surface(
    ray_origin: &Point3,
    ray_dir: &Vec3,
    surface: &SurfaceDef,
) -> Option<RayHit> {
    match surface {
        SurfaceDef::Plane(s) => ray_plane(ray_origin, ray_dir, s),
        SurfaceDef::Sphere(s) => ray_sphere(ray_origin, ray_dir, s),
        SurfaceDef::Cylinder(s) => ray_cylinder(ray_origin, ray_dir, s),
        SurfaceDef::Cone(s) => ray_cone(ray_origin, ray_dir, s),
        SurfaceDef::Torus(_) | SurfaceDef::Nurbs(_) => None, // not yet supported
    }
}

// ===========================================================================
// Triangle-Triangle Intersection
// ===========================================================================

/// Möller's triangle-triangle intersection algorithm.
/// Returns intersection segment endpoints, or None.
pub fn tri_tri_intersect(
    a0: &Point3,
    a1: &Point3,
    a2: &Point3,
    b0: &Point3,
    b1: &Point3,
    b2: &Point3,
) -> TriTriResult {
    // Plane of triangle A
    let na = (a1 - a0).cross(&(a2 - a0));
    if na.norm_squared() < 1e-24 {
        return TriTriResult::no_intersection();
    }
    let da = na.dot(&a0.coords);

    // Signed distances of B's vertices to plane A
    let db = [
        snap(na.dot(&b0.coords) - da),
        snap(na.dot(&b1.coords) - da),
        snap(na.dot(&b2.coords) - da),
    ];
    if (db[0] > 0.0 && db[1] > 0.0 && db[2] > 0.0)
        || (db[0] < 0.0 && db[1] < 0.0 && db[2] < 0.0)
    {
        return TriTriResult::no_intersection();
    }

    // Plane of triangle B
    let nb = (b1 - b0).cross(&(b2 - b0));
    if nb.norm_squared() < 1e-24 {
        return TriTriResult::no_intersection();
    }
    let db_plane = nb.dot(&b0.coords);

    // Signed distances of A's vertices to plane B
    let da_arr = [
        snap(nb.dot(&a0.coords) - db_plane),
        snap(nb.dot(&a1.coords) - db_plane),
        snap(nb.dot(&a2.coords) - db_plane),
    ];
    if (da_arr[0] > 0.0 && da_arr[1] > 0.0 && da_arr[2] > 0.0)
        || (da_arr[0] < 0.0 && da_arr[1] < 0.0 && da_arr[2] < 0.0)
    {
        return TriTriResult::no_intersection();
    }

    // Intersection line direction
    let dir = na.cross(&nb);
    let dir_len = dir.norm();
    if dir_len < 1e-15 {
        // Coplanar — not handled (return no intersection)
        return TriTriResult::no_intersection();
    }
    let dir_unit = dir / dir_len;

    // Crossing points: where each triangle crosses the other's plane
    let a_cross = triangle_plane_crossings([a0, a1, a2], da_arr);
    let b_cross = triangle_plane_crossings([b0, b1, b2], db);

    if a_cross.len() < 2 || b_cross.len() < 2 {
        return if !a_cross.is_empty() && !b_cross.is_empty() {
            TriTriResult {
                intersects: true,
                segment: None,
            }
        } else {
            TriTriResult::no_intersection()
        };
    }

    // Project crossings onto intersection line
    let origin = a_cross[0];
    let ta0 = 0.0_f64; // (a_cross[0] - origin).dot(&dir_unit) = 0
    let ta1 = (a_cross[1] - origin).dot(&dir_unit);
    let tb0 = (b_cross[0] - origin).dot(&dir_unit);
    let tb1 = (b_cross[1] - origin).dot(&dir_unit);

    let (ta_min, ta_max) = if ta0 < ta1 { (ta0, ta1) } else { (ta1, ta0) };
    let (tb_min, tb_max) = if tb0 < tb1 { (tb0, tb1) } else { (tb1, tb0) };

    // Check interval overlap
    let s0 = ta_min.max(tb_min);
    let s1 = ta_max.min(tb_max);

    if s0 > s1 + 1e-12 {
        return TriTriResult::no_intersection();
    }

    let p0 = Point3::from(origin.coords + dir_unit * s0);
    let p1 = Point3::from(origin.coords + dir_unit * s1);

    TriTriResult {
        intersects: true,
        segment: Some((p0, p1)),
    }
}

// ===========================================================================
// Mesh-Level Intersection
// ===========================================================================

/// Triangle mesh vs plane intersection.
/// Returns raw segments from edge-plane crossings.
pub fn mesh_plane_intersect(mesh: &TriMesh, plane: &Plane) -> Vec<(Point3, Point3)> {
    // Precompute signed distances for all vertices
    let dists: Vec<f64> = mesh
        .positions
        .iter()
        .map(|p| (p - plane.origin).dot(&plane.normal))
        .collect();

    let mut segments = Vec::new();

    for tri in &mesh.indices {
        let i0 = tri[0] as usize;
        let i1 = tri[1] as usize;
        let i2 = tri[2] as usize;

        let verts = [&mesh.positions[i0], &mesh.positions[i1], &mesh.positions[i2]];
        let d = [dists[i0], dists[i1], dists[i2]];

        let crossings = triangle_plane_crossings(verts, d);
        if crossings.len() >= 2 {
            segments.push((crossings[0], crossings[1]));
        }
    }

    segments
}

/// Triangle mesh vs triangle mesh intersection.
/// O(n*m) brute force — calls tri_tri_intersect per pair.
pub fn mesh_mesh_intersect(a: &TriMesh, b: &TriMesh) -> Vec<(Point3, Point3)> {
    let mut segments = Vec::new();

    for tri_a in &a.indices {
        let a0 = &a.positions[tri_a[0] as usize];
        let a1 = &a.positions[tri_a[1] as usize];
        let a2 = &a.positions[tri_a[2] as usize];

        for tri_b in &b.indices {
            let b0 = &b.positions[tri_b[0] as usize];
            let b1 = &b.positions[tri_b[1] as usize];
            let b2 = &b.positions[tri_b[2] as usize];

            let result = tri_tri_intersect(a0, a1, a2, b0, b1, b2);
            if let Some((p0, p1)) = result.segment {
                segments.push((p0, p1));
            }
        }
    }

    segments
}

/// Chain raw intersection segments into continuous polylines.
/// Greedy nearest-endpoint chaining within tolerance.
pub fn chain_segments(segments: Vec<(Point3, Point3)>, tolerance: f64) -> Vec<Vec<Point3>> {
    let tol_sq = tolerance * tolerance;
    let mut remaining = segments;
    let mut chains: Vec<Vec<Point3>> = Vec::new();

    while !remaining.is_empty() {
        let (p0, p1) = remaining.swap_remove(0);
        let mut chain = vec![p0, p1];

        // Extend forward
        loop {
            let end = *chain.last().unwrap();
            let mut found = false;

            for i in 0..remaining.len() {
                let d_start = (remaining[i].0 - end).norm_squared();
                let d_end = (remaining[i].1 - end).norm_squared();

                if d_start < tol_sq {
                    let seg = remaining.swap_remove(i);
                    chain.push(seg.1);
                    found = true;
                    break;
                } else if d_end < tol_sq {
                    let seg = remaining.swap_remove(i);
                    chain.push(seg.0);
                    found = true;
                    break;
                }
            }

            if !found {
                break;
            }
        }

        // Extend backward (reverse chain, extend forward, reverse back)
        chain.reverse();
        loop {
            let end = *chain.last().unwrap();
            let mut found = false;

            for i in 0..remaining.len() {
                let d_start = (remaining[i].0 - end).norm_squared();
                let d_end = (remaining[i].1 - end).norm_squared();

                if d_start < tol_sq {
                    let seg = remaining.swap_remove(i);
                    chain.push(seg.1);
                    found = true;
                    break;
                } else if d_end < tol_sq {
                    let seg = remaining.swap_remove(i);
                    chain.push(seg.0);
                    found = true;
                    break;
                }
            }

            if !found {
                break;
            }
        }
        chain.reverse();

        chains.push(chain);
    }

    chains
}

// ===========================================================================
// NURBS / General SSI
// ===========================================================================

/// Plane vs NURBS surface intersection.
/// Tessellate NURBS → mesh_plane_intersect → chain → refine.
pub fn intersect_plane_nurbs(
    p: &Plane,
    ns: &FlippableNurbs,
    tess_divs: usize,
    chain_tol: f64,
    refine_iters: usize,
    refine_tol: f64,
) -> SurfaceSurfaceResult {
    let mesh = tessellate_surface_grid(ns, tess_divs, tess_divs);
    let raw = mesh_plane_intersect(&mesh, p);

    if raw.is_empty() {
        return SurfaceSurfaceResult::Empty;
    }

    let chains = chain_segments(raw, chain_tol);
    let plane_surf: &dyn Surface = p;
    let nurbs_surf: &dyn Surface = ns;

    let curves: Vec<IntersectionCurve> = chains
        .into_iter()
        .map(|pts| {
            let refined = refine_intersection_points(plane_surf, nurbs_surf, &pts, refine_iters, refine_tol);
            IntersectionCurve::Polyline(IntersectionPolyline { points: refined })
        })
        .collect();

    if curves.is_empty() {
        SurfaceSurfaceResult::Empty
    } else {
        SurfaceSurfaceResult::Curves(curves)
    }
}

/// General surface-surface intersection (at least one NURBS, or any pair).
/// Tessellate both → mesh_mesh_intersect → chain → refine.
pub fn intersect_surfaces(
    a: &SurfaceDef,
    b: &SurfaceDef,
    tess_divs: usize,
    chain_tol: f64,
    refine_iters: usize,
    refine_tol: f64,
) -> SurfaceSurfaceResult {
    let mesh_a = tessellate_surface_grid(a, tess_divs, tess_divs);
    let mesh_b = tessellate_surface_grid(b, tess_divs, tess_divs);
    let raw = mesh_mesh_intersect(&mesh_a, &mesh_b);

    if raw.is_empty() {
        return SurfaceSurfaceResult::Empty;
    }

    let chains = chain_segments(raw, chain_tol);
    let surf_a: &dyn Surface = a;
    let surf_b: &dyn Surface = b;

    let curves: Vec<IntersectionCurve> = chains
        .into_iter()
        .map(|pts| {
            let refined = refine_intersection_points(surf_a, surf_b, &pts, refine_iters, refine_tol);
            IntersectionCurve::Polyline(IntersectionPolyline { points: refined })
        })
        .collect();

    if curves.is_empty() {
        SurfaceSurfaceResult::Empty
    } else {
        SurfaceSurfaceResult::Curves(curves)
    }
}

/// Refine approximate intersection points to lie on both surfaces.
/// Alternating projection: project onto A → B → midpoint, repeat.
pub fn refine_intersection_points(
    a: &dyn Surface,
    b: &dyn Surface,
    points: &[Point3],
    max_iter: usize,
    tol: f64,
) -> Vec<Point3> {
    points
        .iter()
        .map(|p| {
            let mut pt = *p;
            for _ in 0..max_iter {
                // Project onto surface A
                let (ua, va) = a.inverse_uv(&pt);
                let pa = a.eval(ua, va);

                // Project onto surface B
                let (ub, vb) = b.inverse_uv(&pa);
                let pb = b.eval(ub, vb);

                // Midpoint
                pt = Point3::from((pa.coords + pb.coords) * 0.5);

                if (pa - pb).norm() < tol {
                    break;
                }
            }
            pt
        })
        .collect()
}

// ===========================================================================
// SSI Pipeline
// ===========================================================================

/// A solver that can handle a specific pair of surface types.
pub trait SsiSolver {
    /// Returns true if this solver can handle (kind_a, kind_b) in either order.
    fn can_solve(&self, kind_a: SurfaceKind, kind_b: SurfaceKind) -> bool;

    /// Solve using concrete surface definitions.
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult;
}

/// Registry of solvers tried in order. Returns first non-Empty result.
pub struct SsiPipeline {
    solvers: Vec<Box<dyn SsiSolver>>,
}

impl SsiPipeline {
    pub fn new() -> Self {
        SsiPipeline {
            solvers: Vec::new(),
        }
    }

    pub fn register(&mut self, solver: Box<dyn SsiSolver>) {
        self.solvers.push(solver);
    }

    pub fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        let ka = a.kind();
        let kb = b.kind();

        for solver in &self.solvers {
            if solver.can_solve(ka, kb) {
                let result = solver.solve(a, b);
                if !matches!(result, SurfaceSurfaceResult::Empty) {
                    return result;
                }
            }
        }

        SurfaceSurfaceResult::Empty
    }
}

impl Default for SsiPipeline {
    fn default() -> Self {
        default_pipeline()
    }
}

// ---------------------------------------------------------------------------
// Concrete solvers
// ---------------------------------------------------------------------------

struct PlanePlaneSolver;
impl SsiSolver for PlanePlaneSolver {
    fn can_solve(&self, a: SurfaceKind, b: SurfaceKind) -> bool {
        a == SurfaceKind::Plane && b == SurfaceKind::Plane
    }
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        if let (SurfaceDef::Plane(pa), SurfaceDef::Plane(pb)) = (a, b) {
            intersect_plane_plane(pa, pb)
        } else {
            SurfaceSurfaceResult::Empty
        }
    }
}

struct PlaneSphereSolver;
impl SsiSolver for PlaneSphereSolver {
    fn can_solve(&self, a: SurfaceKind, b: SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Plane)
        )
    }
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        match (a, b) {
            (SurfaceDef::Plane(p), SurfaceDef::Sphere(s)) => intersect_plane_sphere(p, s),
            (SurfaceDef::Sphere(s), SurfaceDef::Plane(p)) => intersect_plane_sphere(p, s),
            _ => SurfaceSurfaceResult::Empty,
        }
    }
}

struct PlaneCylinderSolver;
impl SsiSolver for PlaneCylinderSolver {
    fn can_solve(&self, a: SurfaceKind, b: SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane, SurfaceKind::Cylinder)
                | (SurfaceKind::Cylinder, SurfaceKind::Plane)
        )
    }
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        match (a, b) {
            (SurfaceDef::Plane(p), SurfaceDef::Cylinder(c)) => intersect_plane_cylinder(p, c),
            (SurfaceDef::Cylinder(c), SurfaceDef::Plane(p)) => intersect_plane_cylinder(p, c),
            _ => SurfaceSurfaceResult::Empty,
        }
    }
}

struct PlaneConeSolver;
impl SsiSolver for PlaneConeSolver {
    fn can_solve(&self, a: SurfaceKind, b: SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Plane)
        )
    }
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        match (a, b) {
            (SurfaceDef::Plane(p), SurfaceDef::Cone(c)) => intersect_plane_cone(p, c),
            (SurfaceDef::Cone(c), SurfaceDef::Plane(p)) => intersect_plane_cone(p, c),
            _ => SurfaceSurfaceResult::Empty,
        }
    }
}

struct PlaneTorusSolver;
impl SsiSolver for PlaneTorusSolver {
    fn can_solve(&self, a: SurfaceKind, b: SurfaceKind) -> bool {
        matches!(
            (a, b),
            (SurfaceKind::Plane, SurfaceKind::Torus) | (SurfaceKind::Torus, SurfaceKind::Plane)
        )
    }
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        match (a, b) {
            (SurfaceDef::Plane(p), SurfaceDef::Torus(t)) => intersect_plane_torus(p, t),
            (SurfaceDef::Torus(t), SurfaceDef::Plane(p)) => intersect_plane_torus(p, t),
            _ => SurfaceSurfaceResult::Empty,
        }
    }
}

/// Fallback: tessellates both surfaces and finds intersection via mesh operations.
struct NurbsFallbackSolver {
    tess_divs: usize,
    chain_tol: f64,
    refine_iters: usize,
    refine_tol: f64,
}

impl Default for NurbsFallbackSolver {
    fn default() -> Self {
        NurbsFallbackSolver {
            tess_divs: 16,
            chain_tol: 0.05,
            refine_iters: 5,
            refine_tol: 1e-8,
        }
    }
}

impl SsiSolver for NurbsFallbackSolver {
    fn can_solve(&self, _a: SurfaceKind, _b: SurfaceKind) -> bool {
        true // can handle any pair
    }
    fn solve(&self, a: &SurfaceDef, b: &SurfaceDef) -> SurfaceSurfaceResult {
        intersect_surfaces(a, b, self.tess_divs, self.chain_tol, self.refine_iters, self.refine_tol)
    }
}

/// Pre-built pipeline with all analytical + NURBS solvers registered.
pub fn default_pipeline() -> SsiPipeline {
    let mut p = SsiPipeline::new();
    p.register(Box::new(PlanePlaneSolver));
    p.register(Box::new(PlaneSphereSolver));
    p.register(Box::new(PlaneCylinderSolver));
    p.register(Box::new(PlaneConeSolver));
    p.register(Box::new(PlaneTorusSolver));
    p.register(Box::new(NurbsFallbackSolver::default()));
    p
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_plane_plane_perpendicular() {
        let a = Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let b = Plane {
            origin: Point3::origin(),
            normal: Vec3::new(1.0, 0.0, 0.0),
        };
        match intersect_plane_plane(&a, &b) {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                if let IntersectionCurve::Line(line) = &curves[0] {
                    // Direction should be along y-axis
                    assert_relative_eq!(line.direction.y.abs(), 1.0, epsilon = 1e-10);
                } else {
                    panic!("Expected line");
                }
            }
            _ => panic!("Expected Curves"),
        }
    }

    #[test]
    fn test_plane_plane_parallel() {
        let a = Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let b = Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        assert!(matches!(intersect_plane_plane(&a, &b), SurfaceSurfaceResult::Empty));
    }

    #[test]
    fn test_plane_plane_coincident() {
        let a = Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let b = Plane {
            origin: Point3::new(5.0, 3.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        assert!(matches!(
            intersect_plane_plane(&a, &b),
            SurfaceSurfaceResult::Coincident
        ));
    }

    #[test]
    fn test_plane_sphere_intersection() {
        let p = Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let s = SphereSurface {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 5.0,
        };
        match intersect_plane_sphere(&p, &s) {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                if let IntersectionCurve::Circle(c) = &curves[0] {
                    assert_relative_eq!(c.radius, 5.0, epsilon = 1e-10);
                    assert_relative_eq!(c.center.z, 0.0, epsilon = 1e-10);
                } else {
                    panic!("Expected circle");
                }
            }
            _ => panic!("Expected Curves"),
        }
    }

    #[test]
    fn test_plane_sphere_miss() {
        let p = Plane {
            origin: Point3::new(0.0, 0.0, 10.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let s = SphereSurface {
            center: Point3::origin(),
            radius: 5.0,
        };
        assert!(matches!(intersect_plane_sphere(&p, &s), SurfaceSurfaceResult::Empty));
    }

    #[test]
    fn test_plane_cylinder_perpendicular() {
        let p = Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let c = CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        };
        match intersect_plane_cylinder(&p, &c) {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                if let IntersectionCurve::Circle(circ) = &curves[0] {
                    assert_relative_eq!(circ.radius, 3.0, epsilon = 1e-10);
                    assert_relative_eq!(circ.center.z, 5.0, epsilon = 1e-10);
                } else {
                    panic!("Expected circle");
                }
            }
            _ => panic!("Expected Curves"),
        }
    }

    #[test]
    fn test_plane_cylinder_oblique() {
        let n = Vec3::new(0.0, 1.0, 1.0).normalize();
        let p = Plane {
            origin: Point3::origin(),
            normal: n,
        };
        let c = CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        match intersect_plane_cylinder(&p, &c) {
            SurfaceSurfaceResult::Curves(curves) => {
                assert_eq!(curves.len(), 1);
                if let IntersectionCurve::Ellipse(e) = &curves[0] {
                    assert_relative_eq!(e.semi_minor, 2.0, epsilon = 1e-10);
                    assert!(e.semi_major > 2.0);
                } else {
                    panic!("Expected ellipse");
                }
            }
            _ => panic!("Expected Curves"),
        }
    }

    #[test]
    fn test_ray_plane_hit() {
        let plane = Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let origin = Point3::origin();
        let dir = Vec3::new(0.0, 0.0, 1.0);
        let hit = ray_plane(&origin, &dir, &plane).unwrap();
        assert_relative_eq!(hit.t, 5.0, epsilon = 1e-10);
        assert_relative_eq!(hit.point.z, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_ray_plane_miss() {
        let plane = Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let origin = Point3::origin();
        let dir = Vec3::new(1.0, 0.0, 0.0); // parallel to plane
        assert!(ray_plane(&origin, &dir, &plane).is_none());
    }

    #[test]
    fn test_ray_sphere_hit() {
        let sphere = SphereSurface {
            center: Point3::new(0.0, 0.0, 10.0),
            radius: 3.0,
        };
        let origin = Point3::origin();
        let dir = Vec3::new(0.0, 0.0, 1.0);
        let hit = ray_sphere(&origin, &dir, &sphere).unwrap();
        assert_relative_eq!(hit.t, 7.0, epsilon = 1e-10);
        assert_relative_eq!(hit.point.z, 7.0, epsilon = 1e-10);
    }

    #[test]
    fn test_ray_sphere_miss() {
        let sphere = SphereSurface {
            center: Point3::new(10.0, 0.0, 0.0),
            radius: 1.0,
        };
        let origin = Point3::origin();
        let dir = Vec3::new(0.0, 0.0, 1.0);
        assert!(ray_sphere(&origin, &dir, &sphere).is_none());
    }

    #[test]
    fn test_ray_cylinder_hit() {
        let cyl = CylinderSurface {
            origin: Point3::new(5.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let origin = Point3::origin();
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let hit = ray_cylinder(&origin, &dir, &cyl).unwrap();
        assert_relative_eq!(hit.t, 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_tri_tri_intersect() {
        // Two triangles that cross each other
        let a0 = Point3::new(-1.0, -1.0, 0.0);
        let a1 = Point3::new(1.0, -1.0, 0.0);
        let a2 = Point3::new(0.0, 1.0, 0.0);

        let b0 = Point3::new(0.0, 0.0, -1.0);
        let b1 = Point3::new(0.0, 0.0, 1.0);
        let b2 = Point3::new(0.0, 2.0, 0.0);

        let result = tri_tri_intersect(&a0, &a1, &a2, &b0, &b1, &b2);
        assert!(result.intersects);
        assert!(result.segment.is_some());
    }

    #[test]
    fn test_tri_tri_no_intersect() {
        // Two separated triangles
        let a0 = Point3::new(0.0, 0.0, 0.0);
        let a1 = Point3::new(1.0, 0.0, 0.0);
        let a2 = Point3::new(0.0, 1.0, 0.0);

        let b0 = Point3::new(0.0, 0.0, 5.0);
        let b1 = Point3::new(1.0, 0.0, 5.0);
        let b2 = Point3::new(0.0, 1.0, 5.0);

        let result = tri_tri_intersect(&a0, &a1, &a2, &b0, &b1, &b2);
        assert!(!result.intersects);
    }

    #[test]
    fn test_chain_segments() {
        let segments = vec![
            (
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ),
            (
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
            ),
            (
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ),
        ];
        let chains = chain_segments(segments, 0.1);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].len(), 4);
    }

    #[test]
    fn test_chain_segments_two_chains() {
        let segments = vec![
            (
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ),
            (
                Point3::new(10.0, 0.0, 0.0),
                Point3::new(11.0, 0.0, 0.0),
            ),
        ];
        let chains = chain_segments(segments, 0.1);
        assert_eq!(chains.len(), 2);
    }

    #[test]
    fn test_ssi_pipeline_plane_plane() {
        let pipeline = default_pipeline();
        let a = SurfaceDef::Plane(Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        });
        let b = SurfaceDef::Plane(Plane {
            origin: Point3::origin(),
            normal: Vec3::new(1.0, 0.0, 0.0),
        });
        let result = pipeline.solve(&a, &b);
        assert!(matches!(result, SurfaceSurfaceResult::Curves(_)));
    }
}

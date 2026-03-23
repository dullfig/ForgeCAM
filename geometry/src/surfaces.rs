// Surface definitions, Surface trait, and operations
// See API.md Section 3

use std::f64::consts::TAU;

use serde::{Deserialize, Serialize};

use curvo::prelude::Transformable;

use crate::mesh::TriMesh;
use crate::types::*;
use crate::utils::{orthonormal_basis, scale_factor, transform_dir, transform_point};

// ---------------------------------------------------------------------------
// Re-export curvo NURBS surface type
// ---------------------------------------------------------------------------

pub use curvo::prelude::NurbsSurface3D;

// ---------------------------------------------------------------------------
// Surface Structs
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Plane {
    pub origin: Point3,
    pub normal: Vec3, // unit normal
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CylinderSurface {
    pub origin: Point3, // point on axis
    pub axis: Vec3,     // unit direction along axis
    pub radius: f64,    // negative = inward-facing normal
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SphereSurface {
    pub center: Point3,
    pub radius: f64, // negative = inward-facing normal
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConeSurface {
    pub apex: Point3,
    pub axis: Vec3,     // unit direction from apex toward base
    pub half_angle: f64, // radians; negative = inward-facing normal
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TorusSurface {
    pub center: Point3,
    pub axis: Vec3,        // unit axis of revolution
    pub major_radius: f64,
    pub minor_radius: f64, // negative = inward-facing normal
}

/// NURBS surface wrapper with orientation flag.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlippableNurbs {
    pub surface: NurbsSurface3D<f64>,
    pub flipped: bool,
}

// ---------------------------------------------------------------------------
// SurfaceDef Enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceDef {
    Plane(Plane),
    Cylinder(CylinderSurface),
    Sphere(SphereSurface),
    Cone(ConeSurface),
    Torus(TorusSurface),
    Nurbs(FlippableNurbs),
}

// ---------------------------------------------------------------------------
// SurfaceKind discriminant (no geometric data)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SurfaceKind {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Nurbs,
    Unknown,
}

impl SurfaceDef {
    pub fn kind(&self) -> SurfaceKind {
        match self {
            SurfaceDef::Plane(_) => SurfaceKind::Plane,
            SurfaceDef::Cylinder(_) => SurfaceKind::Cylinder,
            SurfaceDef::Sphere(_) => SurfaceKind::Sphere,
            SurfaceDef::Cone(_) => SurfaceKind::Cone,
            SurfaceDef::Torus(_) => SurfaceKind::Torus,
            SurfaceDef::Nurbs(_) => SurfaceKind::Nurbs,
        }
    }
}

// ---------------------------------------------------------------------------
// Surface Trait
// ---------------------------------------------------------------------------

pub trait Surface {
    /// Evaluate position at parameters (u, v).
    fn eval(&self, u: f64, v: f64) -> Point3;

    /// Unit outward normal at (u, v).
    fn normal(&self, u: f64, v: f64) -> Vec3;

    /// Parameter domain: ((u_min, u_max), (v_min, v_max)).
    fn domain(&self) -> ((f64, f64), (f64, f64));

    /// Inverse parametric mapping: given a 3D point near the surface,
    /// find the closest (u, v) parameters.
    fn inverse_uv(&self, point: &Point3) -> (f64, f64);
}

// ---------------------------------------------------------------------------
// Surface impl — Plane
// ---------------------------------------------------------------------------

impl Plane {
    /// Get the local (u, v) basis vectors for the plane.
    #[inline]
    fn basis(&self) -> (Vec3, Vec3) {
        orthonormal_basis(&self.normal)
    }
}

impl Surface for Plane {
    fn eval(&self, u: f64, v: f64) -> Point3 {
        let (bu, bv) = self.basis();
        self.origin + bu * u + bv * v
    }

    fn normal(&self, _u: f64, _v: f64) -> Vec3 {
        self.normal
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (
            (f64::NEG_INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, f64::INFINITY),
        )
    }

    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        let (bu, bv) = self.basis();
        let d = point - self.origin;
        (d.dot(&bu), d.dot(&bv))
    }
}

// ---------------------------------------------------------------------------
// Surface impl — Cylinder
// ---------------------------------------------------------------------------

impl CylinderSurface {
    #[inline]
    fn ref_frame(&self) -> (Vec3, Vec3) {
        orthonormal_basis(&self.axis)
    }
}

impl Surface for CylinderSurface {
    /// u = angle [0, 2π], v = axial distance.
    fn eval(&self, u: f64, v: f64) -> Point3 {
        let (u_dir, v_dir) = self.ref_frame();
        let (sin_u, cos_u) = u.sin_cos();
        let r = self.radius.abs();
        self.origin + self.axis * v + r * (cos_u * u_dir + sin_u * v_dir)
    }

    fn normal(&self, u: f64, _v: f64) -> Vec3 {
        let (u_dir, v_dir) = self.ref_frame();
        let (sin_u, cos_u) = u.sin_cos();
        let n = cos_u * u_dir + sin_u * v_dir;
        if self.radius >= 0.0 { n } else { -n }
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), (f64::NEG_INFINITY, f64::INFINITY))
    }

    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        let d = point - self.origin;
        let v = d.dot(&self.axis);
        let lateral = d - v * self.axis;
        let (u_dir, v_dir) = self.ref_frame();
        let mut u = lateral.dot(&v_dir).atan2(lateral.dot(&u_dir));
        if u < 0.0 {
            u += TAU;
        }
        (u, v)
    }
}

// ---------------------------------------------------------------------------
// Surface impl — Sphere
// ---------------------------------------------------------------------------

impl Surface for SphereSurface {
    /// u = azimuth [0, 2π], v = elevation [-π/2, π/2].
    /// Uses z-axis as polar axis.
    fn eval(&self, u: f64, v: f64) -> Point3 {
        let r = self.radius.abs();
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_v, cos_v) = v.sin_cos();
        self.center + Vec3::new(r * cos_v * cos_u, r * cos_v * sin_u, r * sin_v)
    }

    fn normal(&self, u: f64, v: f64) -> Vec3 {
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_v, cos_v) = v.sin_cos();
        let n = Vec3::new(cos_v * cos_u, cos_v * sin_u, sin_v);
        if self.radius >= 0.0 { n } else { -n }
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2))
    }

    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        let d = point - self.center;
        let r = d.norm();
        if r < 1e-15 {
            return (0.0, 0.0);
        }
        let v = (d.z / r).clamp(-1.0, 1.0).asin();
        let mut u = d.y.atan2(d.x);
        if u < 0.0 {
            u += TAU;
        }
        (u, v)
    }
}

// ---------------------------------------------------------------------------
// Surface impl — Cone
// ---------------------------------------------------------------------------

impl ConeSurface {
    #[inline]
    fn ref_frame(&self) -> (Vec3, Vec3) {
        orthonormal_basis(&self.axis)
    }
}

impl Surface for ConeSurface {
    /// u = angle [0, 2π], v = distance along axis from apex.
    fn eval(&self, u: f64, v: f64) -> Point3 {
        let (u_dir, v_dir) = self.ref_frame();
        let (sin_u, cos_u) = u.sin_cos();
        let alpha = self.half_angle.abs();
        let r = v * alpha.tan();
        self.apex + self.axis * v + r * (cos_u * u_dir + sin_u * v_dir)
    }

    fn normal(&self, u: f64, _v: f64) -> Vec3 {
        let (u_dir, v_dir) = self.ref_frame();
        let (sin_u, cos_u) = u.sin_cos();
        let alpha = self.half_angle.abs();
        let radial = cos_u * u_dir + sin_u * v_dir;
        let outward = alpha.cos() * radial - alpha.sin() * self.axis;
        if self.half_angle >= 0.0 {
            outward.normalize()
        } else {
            -outward.normalize()
        }
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), (0.0, f64::INFINITY))
    }

    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        let d = point - self.apex;
        let v = d.dot(&self.axis);
        let lateral = d - v * self.axis;
        let (u_dir, v_dir) = self.ref_frame();
        let mut u = lateral.dot(&v_dir).atan2(lateral.dot(&u_dir));
        if u < 0.0 {
            u += TAU;
        }
        (u, v)
    }
}

// ---------------------------------------------------------------------------
// Surface impl — Torus
// ---------------------------------------------------------------------------

impl TorusSurface {
    #[inline]
    fn ref_frame(&self) -> (Vec3, Vec3) {
        orthonormal_basis(&self.axis)
    }
}

impl Surface for TorusSurface {
    /// u = major angle [0, 2π], v = minor angle [0, 2π].
    fn eval(&self, u: f64, v: f64) -> Point3 {
        let (u_dir, v_dir) = self.ref_frame();
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_v, cos_v) = v.sin_cos();
        let major_dir = cos_u * u_dir + sin_u * v_dir;
        let r_minor = self.minor_radius.abs();
        self.center
            + (self.major_radius + r_minor * cos_v) * major_dir
            + r_minor * sin_v * self.axis
    }

    fn normal(&self, u: f64, v: f64) -> Vec3 {
        let (u_dir, v_dir) = self.ref_frame();
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_v, cos_v) = v.sin_cos();
        let major_dir = cos_u * u_dir + sin_u * v_dir;
        let n = cos_v * major_dir + sin_v * self.axis;
        if self.minor_radius >= 0.0 { n } else { -n }
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), (0.0, TAU))
    }

    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        let d = point - self.center;
        // Project onto equatorial plane
        let d_proj = d - d.dot(&self.axis) * self.axis;
        let (u_dir, v_dir) = self.ref_frame();
        let mut u = d_proj.dot(&v_dir).atan2(d_proj.dot(&u_dir));
        if u < 0.0 {
            u += TAU;
        }

        // Find tube center and compute minor angle
        let proj_norm = d_proj.norm();
        let major_dir = if proj_norm > 1e-15 {
            d_proj / proj_norm
        } else {
            u_dir
        };
        let tube_center = self.center + self.major_radius * major_dir;
        let d_minor = point - tube_center;
        let mut v = d_minor.dot(&self.axis).atan2(d_minor.dot(&major_dir));
        if v < 0.0 {
            v += TAU;
        }
        (u, v)
    }
}

// ---------------------------------------------------------------------------
// Surface impl — FlippableNurbs
// ---------------------------------------------------------------------------

impl Surface for FlippableNurbs {
    fn eval(&self, u: f64, v: f64) -> Point3 {
        let p = self.surface.point_at(u, v);
        Point3::new(p.x, p.y, p.z)
    }

    fn normal(&self, u: f64, v: f64) -> Vec3 {
        let n = self.surface.normal_at(u, v);
        let nv = Vec3::new(n.x, n.y, n.z).normalize();
        if self.flipped { -nv } else { nv }
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        self.surface.knots_domain()
    }

    /// 8x8 grid search + 5-iteration Gauss-Newton with finite-difference Jacobian.
    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        let ((u_min, u_max), (v_min, v_max)) = self.domain();
        let grid = 8;
        let h = 1e-6;

        // Grid search for closest starting point
        let mut best_u = u_min;
        let mut best_v = v_min;
        let mut best_dist = f64::MAX;

        for i in 0..=grid {
            for j in 0..=grid {
                let u = u_min + (u_max - u_min) * i as f64 / grid as f64;
                let v = v_min + (v_max - v_min) * j as f64 / grid as f64;
                let p = self.eval(u, v);
                let dist = (p - point).norm_squared();
                if dist < best_dist {
                    best_dist = dist;
                    best_u = u;
                    best_v = v;
                }
            }
        }

        // Gauss-Newton refinement (5 iterations)
        let mut u = best_u;
        let mut v = best_v;
        for _ in 0..5 {
            let p = self.eval(u, v);
            let f = p - point;

            // Finite-difference Jacobian (3x2)
            let pu = (self.eval(u + h, v) - self.eval(u - h, v)) / (2.0 * h);
            let pv = (self.eval(u, v + h) - self.eval(u, v - h)) / (2.0 * h);

            // J^T * J (2x2) and J^T * f (2x1)
            let a11 = pu.dot(&pu);
            let a12 = pu.dot(&pv);
            let a22 = pv.dot(&pv);
            let b1 = -pu.dot(&f);
            let b2 = -pv.dot(&f);

            // Solve 2x2 system via Cramer's rule
            let det = a11 * a22 - a12 * a12;
            if det.abs() < 1e-30 {
                break;
            }
            let du = (b1 * a22 - b2 * a12) / det;
            let dv = (a11 * b2 - a12 * b1) / det;

            u += du;
            v += dv;

            // Clamp to domain
            u = u.clamp(u_min, u_max);
            v = v.clamp(v_min, v_max);
        }

        // Warn if final error is too large
        let final_dist = (self.eval(u, v) - point).norm();
        if final_dist > 1e-4 {
            tracing::warn!(
                "NURBS inverse_uv: final distance {:.6} > 1e-4",
                final_dist
            );
        }

        (u, v)
    }
}

// ---------------------------------------------------------------------------
// Surface impl — SurfaceDef dispatch
// ---------------------------------------------------------------------------

impl Surface for SurfaceDef {
    fn eval(&self, u: f64, v: f64) -> Point3 {
        match self {
            SurfaceDef::Plane(s) => s.eval(u, v),
            SurfaceDef::Cylinder(s) => s.eval(u, v),
            SurfaceDef::Sphere(s) => s.eval(u, v),
            SurfaceDef::Cone(s) => s.eval(u, v),
            SurfaceDef::Torus(s) => s.eval(u, v),
            SurfaceDef::Nurbs(s) => s.eval(u, v),
        }
    }

    fn normal(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceDef::Plane(s) => s.normal(u, v),
            SurfaceDef::Cylinder(s) => s.normal(u, v),
            SurfaceDef::Sphere(s) => s.normal(u, v),
            SurfaceDef::Cone(s) => s.normal(u, v),
            SurfaceDef::Torus(s) => s.normal(u, v),
            SurfaceDef::Nurbs(s) => s.normal(u, v),
        }
    }

    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        match self {
            SurfaceDef::Plane(s) => s.domain(),
            SurfaceDef::Cylinder(s) => s.domain(),
            SurfaceDef::Sphere(s) => s.domain(),
            SurfaceDef::Cone(s) => s.domain(),
            SurfaceDef::Torus(s) => s.domain(),
            SurfaceDef::Nurbs(s) => s.domain(),
        }
    }

    fn inverse_uv(&self, point: &Point3) -> (f64, f64) {
        match self {
            SurfaceDef::Plane(s) => s.inverse_uv(point),
            SurfaceDef::Cylinder(s) => s.inverse_uv(point),
            SurfaceDef::Sphere(s) => s.inverse_uv(point),
            SurfaceDef::Cone(s) => s.inverse_uv(point),
            SurfaceDef::Torus(s) => s.inverse_uv(point),
            SurfaceDef::Nurbs(s) => s.inverse_uv(point),
        }
    }
}

// ---------------------------------------------------------------------------
// Surface Operations
// ---------------------------------------------------------------------------

impl SurfaceDef {
    /// Flip the surface normal direction.
    pub fn flip_normal(&mut self) {
        match self {
            SurfaceDef::Plane(s) => s.normal = -s.normal,
            SurfaceDef::Cylinder(s) => s.radius = -s.radius,
            SurfaceDef::Sphere(s) => s.radius = -s.radius,
            SurfaceDef::Cone(s) => s.half_angle = -s.half_angle,
            SurfaceDef::Torus(s) => s.minor_radius = -s.minor_radius,
            SurfaceDef::Nurbs(s) => s.flipped = !s.flipped,
        }
    }

    /// Apply affine transform (rigid + uniform scale).
    pub fn apply_transform(&mut self, m: &Mat4) {
        let sf = scale_factor(m);
        match self {
            SurfaceDef::Plane(s) => {
                s.origin = transform_point(m, &s.origin);
                s.normal = transform_dir(m, &s.normal);
            }
            SurfaceDef::Cylinder(s) => {
                s.origin = transform_point(m, &s.origin);
                s.axis = transform_dir(m, &s.axis);
                // Preserve sign (negative = inward normal)
                s.radius = s.radius.signum() * s.radius.abs() * sf;
            }
            SurfaceDef::Sphere(s) => {
                s.center = transform_point(m, &s.center);
                s.radius = s.radius.signum() * s.radius.abs() * sf;
            }
            SurfaceDef::Cone(s) => {
                s.apex = transform_point(m, &s.apex);
                s.axis = transform_dir(m, &s.axis);
                // half_angle unchanged by rigid + uniform scale
            }
            SurfaceDef::Torus(s) => {
                s.center = transform_point(m, &s.center);
                s.axis = transform_dir(m, &s.axis);
                s.major_radius *= sf;
                s.minor_radius = s.minor_radius.signum() * s.minor_radius.abs() * sf;
            }
            SurfaceDef::Nurbs(s) => {
                s.surface = s.surface.transformed(m);
            }
        }
    }

    /// Translate surface by offset vector.
    pub fn translate(&mut self, delta: &Vec3) {
        match self {
            SurfaceDef::Plane(s) => s.origin = Point3::from(s.origin.coords + delta),
            SurfaceDef::Cylinder(s) => s.origin = Point3::from(s.origin.coords + delta),
            SurfaceDef::Sphere(s) => s.center = Point3::from(s.center.coords + delta),
            SurfaceDef::Cone(s) => s.apex = Point3::from(s.apex.coords + delta),
            SurfaceDef::Torus(s) => s.center = Point3::from(s.center.coords + delta),
            SurfaceDef::Nurbs(s) => {
                let m = Mat4::new_translation(delta);
                s.surface = s.surface.transformed(&m);
            }
        }
    }

    /// Return a new surface offset by `distance` along the surface normal.
    ///
    /// Positive distance offsets outward (in normal direction), negative inward.
    /// The sign convention (negative radius = flipped normal) is preserved.
    ///
    /// Returns `None` if the offset degenerates the surface (e.g., sphere radius → 0,
    /// cone with near-zero half-angle).  NURBS offset is not yet implemented.
    pub fn offset(&self, distance: f64) -> Option<SurfaceDef> {
        match self {
            SurfaceDef::Plane(p) => {
                let n = p.normal.normalize();
                Some(SurfaceDef::Plane(Plane {
                    origin: p.origin + n * distance,
                    normal: p.normal,
                }))
            }
            SurfaceDef::Cylinder(c) => {
                let new_radius = c.radius + distance;
                if new_radius.abs() < 1e-12 {
                    return None;
                }
                Some(SurfaceDef::Cylinder(CylinderSurface {
                    origin: c.origin,
                    axis: c.axis,
                    radius: new_radius,
                }))
            }
            SurfaceDef::Sphere(s) => {
                let new_radius = s.radius + distance;
                if new_radius.abs() < 1e-12 {
                    return None;
                }
                Some(SurfaceDef::Sphere(SphereSurface {
                    center: s.center,
                    radius: new_radius,
                }))
            }
            SurfaceDef::Cone(c) => {
                let sin_a = c.half_angle.sin();
                if sin_a.abs() < 1e-12 {
                    return None;
                }
                let axis_unit = c.axis.normalize();
                let new_apex = c.apex - axis_unit * (distance / sin_a);
                Some(SurfaceDef::Cone(ConeSurface {
                    apex: new_apex,
                    axis: c.axis,
                    half_angle: c.half_angle,
                }))
            }
            SurfaceDef::Torus(t) => {
                let new_minor = t.minor_radius + distance;
                if new_minor.abs() < 1e-12 {
                    return None;
                }
                Some(SurfaceDef::Torus(TorusSurface {
                    center: t.center,
                    axis: t.axis,
                    major_radius: t.major_radius,
                    minor_radius: new_minor,
                }))
            }
            SurfaceDef::Nurbs(ns) => nurbs_offset(ns, distance),
        }
    }

    /// Generate a triangle mesh of the surface over its domain.
    /// Returns None for analytical surfaces (tessellation is driven by face boundaries
    /// in the kernel). Returns Some(TriMesh) for NURBS.
    pub fn tessellate(&self, divs_u: usize, divs_v: usize) -> Option<TriMesh> {
        match self {
            SurfaceDef::Nurbs(s) => Some(crate::mesh::tessellate_surface_grid(s, divs_u, divs_v)),
            _ => None,
        }
    }
}

/// Compute AABB inflation amount for a curved surface.
pub fn surface_aabb_inflation(surface: &SurfaceDef) -> f64 {
    match surface {
        SurfaceDef::Plane(_) => 0.0,
        SurfaceDef::Cylinder(c) => c.radius.abs(),
        SurfaceDef::Sphere(s) => s.radius.abs(),
        SurfaceDef::Cone(c) => {
            // Conservative: use tangent of half-angle times a scale factor
            c.half_angle.abs().tan() * 10.0
        }
        SurfaceDef::Torus(t) => t.minor_radius.abs(),
        SurfaceDef::Nurbs(_) => 0.0, // NURBS AABB is computed from control points
    }
}

// ---------------------------------------------------------------------------
// NURBS surface offset (approximate via sample + loft)
// ---------------------------------------------------------------------------

/// Offset a NURBS surface by sampling S(u,v) + d*N(u,v), then refitting via loft.
fn nurbs_offset(ns: &FlippableNurbs, distance: f64) -> Option<SurfaceDef> {
    use curvo::prelude::Interpolation;
    use nalgebra::Point3 as NPoint3;

    let ((u_min, u_max), (v_min, v_max)) = ns.domain();

    // Grid density: at least 16, scale up with control point complexity
    let divs = 32_usize;

    // Sample offset points on a grid
    let mut grid: Vec<Vec<NPoint3<f64>>> = Vec::with_capacity(divs + 1);

    for j in 0..=divs {
        let v = v_min + (v_max - v_min) * j as f64 / divs as f64;
        let mut row = Vec::with_capacity(divs + 1);

        for i in 0..=divs {
            let u = u_min + (u_max - u_min) * i as f64 / divs as f64;
            let p = ns.eval(u, v);
            let n = ns.normal(u, v);
            let n_len = n.norm();

            // Validate normal — skip degenerate points by interpolating later
            let offset_n = if n_len > 1e-10 {
                n / n_len
            } else {
                // Degenerate normal (pole/seam). Use a neighbor's normal.
                // Try small perturbation in u.
                let h = (u_max - u_min) * 0.001;
                let n_alt = ns.normal(u + h, v);
                let n_alt_len = n_alt.norm();
                if n_alt_len > 1e-10 {
                    n_alt / n_alt_len
                } else {
                    // Last resort: surface z-axis
                    tracing::warn!(
                        "NURBS offset: degenerate normal at u={:.4}, v={:.4}",
                        u, v
                    );
                    Vec3::new(0.0, 0.0, 1.0)
                }
            };

            let offset_p = p + offset_n * distance;
            row.push(NPoint3::new(offset_p.x, offset_p.y, offset_p.z));
        }

        grid.push(row);
    }

    // Interpolate each v-row as a NURBS curve in the u direction
    let degree = 3.min(divs);
    let mut u_curves: Vec<crate::curves::NurbsCurve3D<f64>> = Vec::with_capacity(divs + 1);

    for row in &grid {
        match crate::curves::NurbsCurve3D::interpolate(row, degree) {
            Ok(curve) => u_curves.push(curve),
            Err(e) => {
                tracing::warn!("NURBS offset: curve interpolation failed: {}", e);
                return None;
            }
        }
    }

    // Loft the u-curves into a surface
    match NurbsSurface3D::try_loft(&u_curves, Some(degree)) {
        Ok(surface) => Some(SurfaceDef::Nurbs(FlippableNurbs {
            surface,
            flipped: ns.flipped,
        })),
        Err(e) => {
            tracing::warn!("NURBS offset: loft failed: {}", e);
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn test_plane_eval() {
        let plane = Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let p = plane.eval(1.0, 2.0);
        // Should be in the XY plane
        assert_relative_eq!(p.z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_plane_normal() {
        let plane = Plane {
            origin: Point3::origin(),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let n = plane.normal(0.0, 0.0);
        assert_relative_eq!(n.z, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_plane_inverse_uv_roundtrip() {
        let plane = Plane {
            origin: Point3::new(1.0, 2.0, 3.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let u = 5.0;
        let v = -3.0;
        let p = plane.eval(u, v);
        let (u2, v2) = plane.inverse_uv(&p);
        assert_relative_eq!(u, u2, epsilon = 1e-10);
        assert_relative_eq!(v, v2, epsilon = 1e-10);
    }

    #[test]
    fn test_cylinder_eval_and_normal() {
        let cyl = CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let p = cyl.eval(0.0, 0.0);
        // At u=0, v=0: should be at distance |radius| from axis
        let lateral = Vec2::new(p.x, p.y);
        assert_relative_eq!(lateral.norm(), 2.0, epsilon = 1e-10);

        let n = cyl.normal(0.0, 0.0);
        assert_relative_eq!(n.norm(), 1.0, epsilon = 1e-10);
        // Normal should point radially outward
        assert_relative_eq!(n.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cylinder_negative_radius_flips_normal() {
        let cyl_pos = CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let cyl_neg = CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: -1.0,
        };
        let n_pos = cyl_pos.normal(0.0, 0.0);
        let n_neg = cyl_neg.normal(0.0, 0.0);
        assert_relative_eq!(n_pos.x, -n_neg.x, epsilon = 1e-12);
        assert_relative_eq!(n_pos.y, -n_neg.y, epsilon = 1e-12);
    }

    #[test]
    fn test_cylinder_inverse_uv_roundtrip() {
        let cyl = CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        };
        let u = 1.2;
        let v = 5.0;
        let p = cyl.eval(u, v);
        let (u2, v2) = cyl.inverse_uv(&p);
        assert_relative_eq!(u, u2, epsilon = 1e-10);
        assert_relative_eq!(v, v2, epsilon = 1e-10);
    }

    #[test]
    fn test_sphere_eval_poles() {
        let sphere = SphereSurface {
            center: Point3::origin(),
            radius: 1.0,
        };
        // North pole: v = pi/2
        let north = sphere.eval(0.0, FRAC_PI_2);
        assert_relative_eq!(north.z, 1.0, epsilon = 1e-10);
        // South pole: v = -pi/2
        let south = sphere.eval(0.0, -FRAC_PI_2);
        assert_relative_eq!(south.z, -1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_sphere_inverse_uv_roundtrip() {
        let sphere = SphereSurface {
            center: Point3::origin(),
            radius: 2.0,
        };
        let u = 1.0;
        let v = 0.3;
        let p = sphere.eval(u, v);
        let (u2, v2) = sphere.inverse_uv(&p);
        assert_relative_eq!(u, u2, epsilon = 1e-10);
        assert_relative_eq!(v, v2, epsilon = 1e-10);
    }

    #[test]
    fn test_cone_eval() {
        let cone = ConeSurface {
            apex: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: PI / 4.0, // 45 degrees
        };
        // At v=1, u=0: r = tan(45°) = 1, so point should be at distance 1 from axis, z=1
        let p = cone.eval(0.0, 1.0);
        assert_relative_eq!(p.z, 1.0, epsilon = 1e-10);
        let lateral = (p.x * p.x + p.y * p.y).sqrt();
        assert_relative_eq!(lateral, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_torus_eval() {
        let torus = TorusSurface {
            center: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // At u=0, v=0: outermost point
        let p = torus.eval(0.0, 0.0);
        let dist_from_axis = (p.x * p.x + p.y * p.y).sqrt();
        assert_relative_eq!(dist_from_axis, 6.0, epsilon = 1e-10);
        assert_relative_eq!(p.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_torus_inverse_uv_roundtrip() {
        let torus = TorusSurface {
            center: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let u = 0.7;
        let v = 2.1;
        let p = torus.eval(u, v);
        let (u2, v2) = torus.inverse_uv(&p);
        assert_relative_eq!(u, u2, epsilon = 1e-8);
        assert_relative_eq!(v, v2, epsilon = 1e-8);
    }

    #[test]
    fn test_flip_normal() {
        let mut sd = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        });
        let n1 = sd.normal(0.0, 0.0);
        sd.flip_normal();
        let n2 = sd.normal(0.0, 0.0);
        assert_relative_eq!(n1.dot(&n2), -1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_surface_aabb_inflation() {
        let plane = SurfaceDef::Plane(Plane {
            origin: Point3::origin(),
            normal: Vec3::z(),
        });
        assert_relative_eq!(surface_aabb_inflation(&plane), 0.0);

        let cyl = SurfaceDef::Cylinder(CylinderSurface {
            origin: Point3::origin(),
            axis: Vec3::z(),
            radius: -3.0,
        });
        assert_relative_eq!(surface_aabb_inflation(&cyl), 3.0);
    }
}

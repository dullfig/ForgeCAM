// Curve definitions, Curve trait, and operations
// See API.md Section 2

use std::f64::consts::TAU;

use serde::{Deserialize, Serialize};

use curvo::prelude::{TrimRange, Transformable};

use crate::types::*;
use crate::utils::{scale_factor, transform_dir, transform_point, Aabb};

// ---------------------------------------------------------------------------
// Re-export curvo NURBS curve type
// ---------------------------------------------------------------------------

pub use curvo::prelude::NurbsCurve3D;

// ---------------------------------------------------------------------------
// Curve Structs
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineSegment {
    pub start: Point3,
    pub end: Point3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircleCurve {
    pub center: Point3,
    pub axis: Vec3,    // unit normal to circle plane
    pub radius: f64,
    pub ref_dir: Vec3, // unit vector in circle plane, defines angle=0
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircularArcCurve {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub ref_dir: Vec3,
    pub start_angle: f64, // radians
    pub end_angle: f64,   // radians
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EllipseCurve {
    pub center: Point3,
    pub axis: Vec3,      // unit normal to ellipse plane
    pub semi_major: f64,
    pub semi_minor: f64,
    pub major_dir: Vec3, // unit vector along semi-major axis
}

// ---------------------------------------------------------------------------
// CurveEnd — identifies which end of a curve
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CurveEnd {
    Start,
    End,
}

// ---------------------------------------------------------------------------
// CurveDef Enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CurveDef {
    LineSegment(LineSegment),
    Circle(CircleCurve),
    CircularArc(CircularArcCurve),
    Ellipse(EllipseCurve),
    Nurbs(NurbsCurve3D<f64>),
}

// ---------------------------------------------------------------------------
// CurveKind discriminant (no geometric data)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CurveKind {
    LineSegment,
    Circle,
    CircularArc,
    Ellipse,
    Nurbs,
    Unknown,
}

impl CurveDef {
    pub fn kind(&self) -> CurveKind {
        match self {
            CurveDef::LineSegment(_) => CurveKind::LineSegment,
            CurveDef::Circle(_) => CurveKind::Circle,
            CurveDef::CircularArc(_) => CurveKind::CircularArc,
            CurveDef::Ellipse(_) => CurveKind::Ellipse,
            CurveDef::Nurbs(_) => CurveKind::Nurbs,
        }
    }
}

// ---------------------------------------------------------------------------
// Curve Trait
// ---------------------------------------------------------------------------

pub trait Curve {
    /// Evaluate position at parameter t.
    fn eval(&self, t: f64) -> Point3;

    /// Unit tangent vector at parameter t.
    fn tangent(&self, t: f64) -> Vec3;

    /// Parameter domain: (t_min, t_max).
    fn domain(&self) -> (f64, f64);
}

// ---------------------------------------------------------------------------
// Curve impls — LineSegment
// ---------------------------------------------------------------------------

impl Curve for LineSegment {
    /// t in [0, 1]: linear interpolation start→end.
    fn eval(&self, t: f64) -> Point3 {
        let d = self.end - self.start;
        self.start + d * t
    }

    fn tangent(&self, _t: f64) -> Vec3 {
        let d = self.end - self.start;
        let n = d.norm();
        if n > 1e-15 {
            d / n
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        }
    }

    fn domain(&self) -> (f64, f64) {
        (0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Curve impls — CircleCurve
// ---------------------------------------------------------------------------

impl CircleCurve {
    /// Perpendicular direction in circle plane: axis × ref_dir.
    #[inline]
    fn perp(&self) -> Vec3 {
        self.axis.cross(&self.ref_dir)
    }
}

impl Curve for CircleCurve {
    /// t in [0, 1] maps to angle [0, 2π].
    fn eval(&self, t: f64) -> Point3 {
        let theta = t * TAU;
        let (sin_t, cos_t) = theta.sin_cos();
        self.center + self.radius * (cos_t * self.ref_dir + sin_t * self.perp())
    }

    fn tangent(&self, t: f64) -> Vec3 {
        let theta = t * TAU;
        let (sin_t, cos_t) = theta.sin_cos();
        let raw = -sin_t * self.ref_dir + cos_t * self.perp();
        raw.normalize()
    }

    fn domain(&self) -> (f64, f64) {
        (0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Curve impls — CircularArcCurve
// ---------------------------------------------------------------------------

impl CircularArcCurve {
    #[inline]
    fn perp(&self) -> Vec3 {
        self.axis.cross(&self.ref_dir)
    }
}

impl Curve for CircularArcCurve {
    /// t in [0, 1] maps to [start_angle, end_angle].
    fn eval(&self, t: f64) -> Point3 {
        let theta = self.start_angle + t * (self.end_angle - self.start_angle);
        let (sin_t, cos_t) = theta.sin_cos();
        self.center + self.radius * (cos_t * self.ref_dir + sin_t * self.perp())
    }

    fn tangent(&self, t: f64) -> Vec3 {
        let theta = self.start_angle + t * (self.end_angle - self.start_angle);
        let (sin_t, cos_t) = theta.sin_cos();
        let raw = -sin_t * self.ref_dir + cos_t * self.perp();
        raw.normalize()
    }

    fn domain(&self) -> (f64, f64) {
        (0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Curve impls — EllipseCurve
// ---------------------------------------------------------------------------

impl EllipseCurve {
    /// Perpendicular direction: axis × major_dir (semi-minor direction).
    #[inline]
    fn minor_dir(&self) -> Vec3 {
        self.axis.cross(&self.major_dir)
    }
}

impl Curve for EllipseCurve {
    /// t in [0, 1] maps to angle [0, 2π].
    fn eval(&self, t: f64) -> Point3 {
        let theta = t * TAU;
        let (sin_t, cos_t) = theta.sin_cos();
        self.center + self.semi_major * cos_t * self.major_dir
            + self.semi_minor * sin_t * self.minor_dir()
    }

    fn tangent(&self, t: f64) -> Vec3 {
        let theta = t * TAU;
        let (sin_t, cos_t) = theta.sin_cos();
        let raw =
            -self.semi_major * sin_t * self.major_dir + self.semi_minor * cos_t * self.minor_dir();
        raw.normalize()
    }

    fn domain(&self) -> (f64, f64) {
        (0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Curve impls — NurbsCurve3D (delegates to curvo)
// ---------------------------------------------------------------------------

impl Curve for NurbsCurve3D<f64> {
    fn eval(&self, t: f64) -> Point3 {
        let p = self.point_at(t);
        Point3::new(p.x, p.y, p.z)
    }

    fn tangent(&self, t: f64) -> Vec3 {
        let d = self.tangent_at(t);
        let v = Vec3::new(d.x, d.y, d.z);
        let n = v.norm();
        if n > 1e-15 {
            v / n
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        }
    }

    fn domain(&self) -> (f64, f64) {
        self.knots_domain()
    }
}

// ---------------------------------------------------------------------------
// Curve impls — CurveDef dispatch
// ---------------------------------------------------------------------------

impl Curve for CurveDef {
    fn eval(&self, t: f64) -> Point3 {
        match self {
            CurveDef::LineSegment(c) => c.eval(t),
            CurveDef::Circle(c) => c.eval(t),
            CurveDef::CircularArc(c) => c.eval(t),
            CurveDef::Ellipse(c) => c.eval(t),
            CurveDef::Nurbs(c) => c.eval(t),
        }
    }

    fn tangent(&self, t: f64) -> Vec3 {
        match self {
            CurveDef::LineSegment(c) => c.tangent(t),
            CurveDef::Circle(c) => c.tangent(t),
            CurveDef::CircularArc(c) => c.tangent(t),
            CurveDef::Ellipse(c) => c.tangent(t),
            CurveDef::Nurbs(c) => c.tangent(t),
        }
    }

    fn domain(&self) -> (f64, f64) {
        match self {
            CurveDef::LineSegment(c) => c.domain(),
            CurveDef::Circle(c) => c.domain(),
            CurveDef::CircularArc(c) => c.domain(),
            CurveDef::Ellipse(c) => c.domain(),
            CurveDef::Nurbs(c) => c.domain(),
        }
    }
}

// ---------------------------------------------------------------------------
// Curve Operations (methods on CurveDef)
// ---------------------------------------------------------------------------

impl CurveDef {
    /// Uniform parameter sampling — returns n+1 points from domain start to end.
    pub fn sample(&self, n: usize) -> Vec<Point3> {
        let (t0, t1) = self.domain();
        (0..=n)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                self.eval(t)
            })
            .collect()
    }

    /// Apply affine transform (rigid + uniform scale).
    pub fn apply_transform(&mut self, m: &Mat4) {
        match self {
            CurveDef::LineSegment(c) => {
                c.start = transform_point(m, &c.start);
                c.end = transform_point(m, &c.end);
            }
            CurveDef::Circle(c) => {
                c.center = transform_point(m, &c.center);
                c.axis = transform_dir(m, &c.axis);
                c.ref_dir = transform_dir(m, &c.ref_dir);
                c.radius *= scale_factor(m);
            }
            CurveDef::CircularArc(c) => {
                c.center = transform_point(m, &c.center);
                c.axis = transform_dir(m, &c.axis);
                c.ref_dir = transform_dir(m, &c.ref_dir);
                c.radius *= scale_factor(m);
            }
            CurveDef::Ellipse(c) => {
                c.center = transform_point(m, &c.center);
                c.axis = transform_dir(m, &c.axis);
                c.major_dir = transform_dir(m, &c.major_dir);
                c.semi_major *= scale_factor(m);
                c.semi_minor *= scale_factor(m);
            }
            CurveDef::Nurbs(c) => {
                *c = c.transformed(m);
            }
        }
    }

    /// Translate by offset vector.
    pub fn translate(&mut self, delta: &Vec3) {
        match self {
            CurveDef::LineSegment(c) => {
                c.start = Point3::from(c.start.coords + delta);
                c.end = Point3::from(c.end.coords + delta);
            }
            CurveDef::Circle(c) => {
                c.center = Point3::from(c.center.coords + delta);
            }
            CurveDef::CircularArc(c) => {
                c.center = Point3::from(c.center.coords + delta);
            }
            CurveDef::Ellipse(c) => {
                c.center = Point3::from(c.center.coords + delta);
            }
            CurveDef::Nurbs(c) => {
                let m = Mat4::new_translation(delta);
                *c = c.transformed(&m);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Split / Trim / Extend / Reverse
    // -----------------------------------------------------------------------

    /// Split curve at parameter t. Returns (left, right).
    pub fn split(&self, t: f64) -> (CurveDef, CurveDef) {
        match self {
            CurveDef::LineSegment(c) => {
                let mid = c.eval(t);
                (
                    CurveDef::LineSegment(LineSegment { start: c.start, end: mid }),
                    CurveDef::LineSegment(LineSegment { start: mid, end: c.end }),
                )
            }
            CurveDef::Circle(c) => {
                let theta = t * TAU;
                (
                    CurveDef::CircularArc(CircularArcCurve {
                        center: c.center, axis: c.axis, radius: c.radius,
                        ref_dir: c.ref_dir, start_angle: 0.0, end_angle: theta,
                    }),
                    CurveDef::CircularArc(CircularArcCurve {
                        center: c.center, axis: c.axis, radius: c.radius,
                        ref_dir: c.ref_dir, start_angle: theta, end_angle: TAU,
                    }),
                )
            }
            CurveDef::CircularArc(c) => {
                let theta = c.start_angle + t * (c.end_angle - c.start_angle);
                (
                    CurveDef::CircularArc(CircularArcCurve {
                        center: c.center, axis: c.axis, radius: c.radius,
                        ref_dir: c.ref_dir, start_angle: c.start_angle, end_angle: theta,
                    }),
                    CurveDef::CircularArc(CircularArcCurve {
                        center: c.center, axis: c.axis, radius: c.radius,
                        ref_dir: c.ref_dir, start_angle: theta, end_angle: c.end_angle,
                    }),
                )
            }
            CurveDef::Ellipse(_) => {
                // Ellipse split produces two sampled NURBS curves (no EllipticalArc type)
                let n = 32;
                let (t0, t1) = self.domain();
                let t_split = t0 + t * (t1 - t0);
                let left = self.sample_range_to_nurbs(t0, t_split, n);
                let right = self.sample_range_to_nurbs(t_split, t1, n);
                (left, right)
            }
            CurveDef::Nurbs(c) => {
                let (d0, d1) = c.knots_domain();
                let t_split = d0 + t * (d1 - d0);
                let left = c.try_trim_range((d0, t_split))
                    .ok()
                    .and_then(|v| v.into_iter().next())
                    .map(CurveDef::Nurbs);
                let right = c.try_trim_range((t_split, d1))
                    .ok()
                    .and_then(|v| v.into_iter().next())
                    .map(CurveDef::Nurbs);
                match (left, right) {
                    (Some(l), Some(r)) => (l, r),
                    _ => {
                        let n = 32;
                        (
                            self.sample_range_to_nurbs(d0, t_split, n),
                            self.sample_range_to_nurbs(t_split, d1, n),
                        )
                    }
                }
            }
        }
    }

    /// Trim curve to parameter range [t0, t1].
    pub fn trim(&self, t0: f64, t1: f64) -> CurveDef {
        match self {
            CurveDef::LineSegment(c) => {
                CurveDef::LineSegment(LineSegment {
                    start: c.eval(t0),
                    end: c.eval(t1),
                })
            }
            CurveDef::Circle(c) => {
                let a0 = t0 * TAU;
                let a1 = t1 * TAU;
                CurveDef::CircularArc(CircularArcCurve {
                    center: c.center, axis: c.axis, radius: c.radius,
                    ref_dir: c.ref_dir, start_angle: a0, end_angle: a1,
                })
            }
            CurveDef::CircularArc(c) => {
                let span = c.end_angle - c.start_angle;
                CurveDef::CircularArc(CircularArcCurve {
                    center: c.center, axis: c.axis, radius: c.radius,
                    ref_dir: c.ref_dir,
                    start_angle: c.start_angle + t0 * span,
                    end_angle: c.start_angle + t1 * span,
                })
            }
            CurveDef::Ellipse(_) => {
                let (d0, d1) = self.domain();
                let ta = d0 + t0 * (d1 - d0);
                let tb = d0 + t1 * (d1 - d0);
                self.sample_range_to_nurbs(ta, tb, 32)
            }
            CurveDef::Nurbs(c) => {
                let (d0, d1) = c.knots_domain();
                let ta = d0 + t0 * (d1 - d0);
                let tb = d0 + t1 * (d1 - d0);
                match c.try_trim_range((ta, tb)) {
                    Ok(v) if !v.is_empty() => CurveDef::Nurbs(v.into_iter().next().unwrap()),
                    _ => self.sample_range_to_nurbs(ta, tb, 32),
                }
            }
        }
    }

    /// Extend curve beyond its domain by distance d at start or end.
    pub fn extend(&self, d: f64, end: CurveEnd) -> CurveDef {
        match self {
            CurveDef::LineSegment(c) => {
                let dir = (c.end - c.start).normalize();
                match end {
                    CurveEnd::Start => CurveDef::LineSegment(LineSegment {
                        start: Point3::from(c.start.coords - dir * d),
                        end: c.end,
                    }),
                    CurveEnd::End => CurveDef::LineSegment(LineSegment {
                        start: c.start,
                        end: Point3::from(c.end.coords + dir * d),
                    }),
                }
            }
            CurveDef::CircularArc(c) => {
                let d_angle = d / c.radius.abs();
                match end {
                    CurveEnd::Start => CurveDef::CircularArc(CircularArcCurve {
                        center: c.center, axis: c.axis, radius: c.radius,
                        ref_dir: c.ref_dir,
                        start_angle: c.start_angle - d_angle,
                        end_angle: c.end_angle,
                    }),
                    CurveEnd::End => CurveDef::CircularArc(CircularArcCurve {
                        center: c.center, axis: c.axis, radius: c.radius,
                        ref_dir: c.ref_dir,
                        start_angle: c.start_angle,
                        end_angle: c.end_angle + d_angle,
                    }),
                }
            }
            _ => {
                // Generic: tangent-linear extension
                let (t0, t1) = self.domain();
                match end {
                    CurveEnd::Start => {
                        let p = self.eval(t0);
                        let tang = self.tangent(t0);
                        let new_start = Point3::from(p.coords - tang * d);
                        // Prepend a line segment
                        // For simplicity, return just the original with extended endpoints
                        // by sampling to NURBS
                        let mut pts = vec![new_start];
                        pts.extend(self.sample(32));
                        self.fit_nurbs_from_points(&pts)
                    }
                    CurveEnd::End => {
                        let p = self.eval(t1);
                        let tang = self.tangent(t1);
                        let new_end = Point3::from(p.coords + tang * d);
                        let mut pts = self.sample(32);
                        pts.push(new_end);
                        self.fit_nurbs_from_points(&pts)
                    }
                }
            }
        }
    }

    /// Reverse curve direction (swap start and end).
    pub fn reverse(&mut self) {
        match self {
            CurveDef::LineSegment(c) => std::mem::swap(&mut c.start, &mut c.end),
            CurveDef::Circle(c) => c.axis = -c.axis,
            CurveDef::CircularArc(c) => {
                std::mem::swap(&mut c.start_angle, &mut c.end_angle);
            }
            CurveDef::Ellipse(c) => c.axis = -c.axis,
            CurveDef::Nurbs(c) => {
                // Sample from the NURBS directly, reverse, refit
                let n = 64;
                let (t0, t1) = c.knots_domain();
                let pts: Vec<Point3> = (0..=n)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * i as f64 / n as f64;
                        let p = c.point_at(t);
                        Point3::new(p.x, p.y, p.z)
                    })
                    .rev()
                    .collect();
                if let Ok(nurbs) = try_interpolate_nurbs(&pts) {
                    *c = nurbs;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Offset
    // -----------------------------------------------------------------------

    /// Offset curve by distance d in the plane perpendicular to `axis`.
    /// Positive d offsets in the direction of tangent x axis (right side
    /// looking along curve direction with axis up — outward for CCW contours).
    pub fn offset(&self, d: f64, axis: &Vec3) -> CurveDef {
        match self {
            CurveDef::LineSegment(c) => {
                let tang = (c.end - c.start).normalize();
                let offset_dir = tang.cross(axis).normalize();
                CurveDef::LineSegment(LineSegment {
                    start: Point3::from(c.start.coords + offset_dir * d),
                    end: Point3::from(c.end.coords + offset_dir * d),
                })
            }
            CurveDef::Circle(c) => {
                let sign = if c.axis.dot(axis) >= 0.0 { 1.0 } else { -1.0 };
                let new_radius = c.radius + sign * d;
                if new_radius.abs() < 1e-15 {
                    // Degenerate — offset to a point. Return zero-length line.
                    CurveDef::LineSegment(LineSegment {
                        start: c.center,
                        end: c.center,
                    })
                } else {
                    CurveDef::Circle(CircleCurve {
                        center: c.center,
                        axis: c.axis,
                        radius: new_radius,
                        ref_dir: c.ref_dir,
                    })
                }
            }
            CurveDef::CircularArc(c) => {
                let sign = if c.axis.dot(axis) >= 0.0 { 1.0 } else { -1.0 };
                let new_radius = c.radius + sign * d;
                if new_radius.abs() < 1e-15 {
                    CurveDef::LineSegment(LineSegment {
                        start: c.center,
                        end: c.center,
                    })
                } else {
                    CurveDef::CircularArc(CircularArcCurve {
                        center: c.center,
                        axis: c.axis,
                        radius: new_radius,
                        ref_dir: c.ref_dir,
                        start_angle: c.start_angle,
                        end_angle: c.end_angle,
                    })
                }
            }
            _ => {
                // General: sample, offset each point along local offset direction
                let n = 64;
                let (t0, t1) = self.domain();
                let mut offset_pts = Vec::with_capacity(n + 1);
                for i in 0..=n {
                    let t = t0 + (t1 - t0) * i as f64 / n as f64;
                    let p = self.eval(t);
                    let tang = self.tangent(t);
                    let offset_dir = tang.cross(axis).normalize();
                    offset_pts.push(Point3::from(p.coords + offset_dir * d));
                }
                self.fit_nurbs_from_points(&offset_pts)
            }
        }
    }

    // -----------------------------------------------------------------------
    // Length / Bounding Box
    // -----------------------------------------------------------------------

    /// Total arc length of the curve.
    pub fn length(&self) -> f64 {
        match self {
            CurveDef::LineSegment(c) => (c.end - c.start).norm(),
            CurveDef::Circle(c) => c.radius.abs() * TAU,
            CurveDef::CircularArc(c) => {
                c.radius.abs() * (c.end_angle - c.start_angle).abs()
            }
            _ => {
                // Chord-length approximation with 256 samples
                let n = 256;
                let (t0, t1) = self.domain();
                let dt = (t1 - t0) / n as f64;
                let mut total = 0.0;
                let mut prev = self.eval(t0);
                for i in 1..=n {
                    let curr = self.eval(t0 + dt * i as f64);
                    total += (curr - prev).norm();
                    prev = curr;
                }
                total
            }
        }
    }

    /// Axis-aligned bounding box.
    pub fn bounding_box(&self) -> Aabb {
        Aabb::from_points(&self.sample(64))
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Sample a parameter sub-range and fit a NURBS through the points.
    fn sample_range_to_nurbs(&self, t0: f64, t1: f64, n: usize) -> CurveDef {
        let pts: Vec<Point3> = (0..=n)
            .map(|i| self.eval(t0 + (t1 - t0) * i as f64 / n as f64))
            .collect();
        self.fit_nurbs_from_points(&pts)
    }

    /// Fit a NURBS curve through a set of points. Falls back to a
    /// line segment if interpolation fails or there are too few points.
    fn fit_nurbs_from_points(&self, pts: &[Point3]) -> CurveDef {
        if pts.len() < 2 {
            return CurveDef::LineSegment(LineSegment {
                start: pts.first().copied().unwrap_or(Point3::origin()),
                end: pts.last().copied().unwrap_or(Point3::origin()),
            });
        }
        match try_interpolate_nurbs(pts) {
            Ok(nurbs) => CurveDef::Nurbs(nurbs),
            Err(_) => {
                // Fallback: return line from first to last point
                CurveDef::LineSegment(LineSegment {
                    start: pts[0],
                    end: *pts.last().unwrap(),
                })
            }
        }
    }
}

/// Try to interpolate a NURBS curve through a set of 3D points.
fn try_interpolate_nurbs(pts: &[Point3]) -> Result<NurbsCurve3D<f64>, ()> {
    use curvo::prelude::Interpolation;
    use nalgebra::Point3 as NPoint3;
    let points: Vec<NPoint3<f64>> = pts
        .iter()
        .map(|p| NPoint3::new(p.x, p.y, p.z))
        .collect();
    let degree = 3.min(points.len() - 1);
    NurbsCurve3D::interpolate(&points, degree)
        .map_err(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_line_segment_eval() {
        let seg = LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(2.0, 0.0, 0.0),
        };
        let mid = seg.eval(0.5);
        assert_relative_eq!(mid.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(mid.y, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_line_segment_tangent() {
        let seg = LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(3.0, 4.0, 0.0),
        };
        let t = seg.tangent(0.0);
        assert_relative_eq!(t.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(t.x, 0.6, epsilon = 1e-12);
        assert_relative_eq!(t.y, 0.8, epsilon = 1e-12);
    }

    #[test]
    fn test_circle_eval_start() {
        let circle = CircleCurve {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            ref_dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let p = circle.eval(0.0);
        assert_relative_eq!(p.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(p.y, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_circle_eval_quarter() {
        let circle = CircleCurve {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            ref_dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let p = circle.eval(0.25); // 90 degrees
        assert_relative_eq!(p.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(p.y, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_circular_arc_eval() {
        let arc = CircularArcCurve {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            ref_dir: Vec3::new(1.0, 0.0, 0.0),
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        // t=0 should be at (2, 0, 0)
        let p0 = arc.eval(0.0);
        assert_relative_eq!(p0.x, 2.0, epsilon = 1e-12);
        assert_relative_eq!(p0.y, 0.0, epsilon = 1e-12);
        // t=1 should be at (0, 2, 0)
        let p1 = arc.eval(1.0);
        assert_relative_eq!(p1.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(p1.y, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn test_ellipse_eval() {
        let ell = EllipseCurve {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            semi_major: 3.0,
            semi_minor: 1.0,
            major_dir: Vec3::new(1.0, 0.0, 0.0),
        };
        // t=0: (3, 0, 0)
        let p0 = ell.eval(0.0);
        assert_relative_eq!(p0.x, 3.0, epsilon = 1e-12);
        assert_relative_eq!(p0.y, 0.0, epsilon = 1e-12);
        // t=0.25: (0, 1, 0)
        let p1 = ell.eval(0.25);
        assert_relative_eq!(p1.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(p1.y, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_sample_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(10.0, 0.0, 0.0),
        });
        let pts = seg.sample(10);
        assert_eq!(pts.len(), 11);
        assert_relative_eq!(pts[0].x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(pts[5].x, 5.0, epsilon = 1e-12);
        assert_relative_eq!(pts[10].x, 10.0, epsilon = 1e-12);
    }

    #[test]
    fn test_translate_line() {
        let mut seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(1.0, 0.0, 0.0),
        });
        seg.translate(&Vec3::new(5.0, 10.0, 0.0));
        let p = seg.eval(0.0);
        assert_relative_eq!(p.x, 5.0, epsilon = 1e-12);
        assert_relative_eq!(p.y, 10.0, epsilon = 1e-12);
    }

    #[test]
    fn test_apply_transform_identity() {
        let mut seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(1.0, 2.0, 3.0),
            end: Point3::new(4.0, 5.0, 6.0),
        });
        seg.apply_transform(&Mat4::identity());
        let p = seg.eval(0.0);
        assert_relative_eq!(p.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(p.y, 2.0, epsilon = 1e-12);
        assert_relative_eq!(p.z, 3.0, epsilon = 1e-12);
    }

    // -- Split / Trim / Reverse / Offset / Length / BBox --

    #[test]
    fn test_split_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(10.0, 0.0, 0.0),
        });
        let (left, right) = seg.split(0.3);
        assert_relative_eq!(left.eval(1.0).x, 3.0, epsilon = 1e-10);
        assert_relative_eq!(right.eval(0.0).x, 3.0, epsilon = 1e-10);
        assert_relative_eq!(right.eval(1.0).x, 10.0, epsilon = 1e-10);
    }

    #[test]
    fn test_split_arc() {
        let arc = CurveDef::CircularArc(CircularArcCurve {
            center: Point3::origin(),
            axis: Vec3::z(),
            radius: 1.0,
            ref_dir: Vec3::x(),
            start_angle: 0.0,
            end_angle: std::f64::consts::PI,
        });
        let (left, right) = arc.split(0.5);
        // Split at 90 degrees
        let mid = left.eval(1.0);
        assert_relative_eq!(mid.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(mid.y, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_trim_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(10.0, 0.0, 0.0),
        });
        let trimmed = seg.trim(0.2, 0.8);
        assert_relative_eq!(trimmed.eval(0.0).x, 2.0, epsilon = 1e-10);
        assert_relative_eq!(trimmed.eval(1.0).x, 8.0, epsilon = 1e-10);
    }

    #[test]
    fn test_reverse_line() {
        let mut seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(10.0, 0.0, 0.0),
        });
        seg.reverse();
        assert_relative_eq!(seg.eval(0.0).x, 10.0, epsilon = 1e-10);
        assert_relative_eq!(seg.eval(1.0).x, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_reverse_arc() {
        let mut arc = CurveDef::CircularArc(CircularArcCurve {
            center: Point3::origin(),
            axis: Vec3::z(),
            radius: 1.0,
            ref_dir: Vec3::x(),
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        });
        let end_before = arc.eval(1.0);
        arc.reverse();
        let start_after = arc.eval(0.0);
        assert_relative_eq!(start_after.x, end_before.x, epsilon = 1e-10);
        assert_relative_eq!(start_after.y, end_before.y, epsilon = 1e-10);
    }

    #[test]
    fn test_offset_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(10.0, 0.0, 0.0),
        });
        let axis = Vec3::z();
        let offset = seg.offset(2.0, &axis);
        // Line along +X, axis = +Z: tangent x axis = (1,0,0)x(0,0,1) = (0,-1,0)? No...
        // (1,0,0) cross (0,0,1) = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0,-1,0)
        // So offset by +2 moves in -Y direction
        assert_relative_eq!(offset.eval(0.0).y, -2.0, epsilon = 1e-10);
        assert_relative_eq!(offset.eval(0.0).x, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_offset_circle() {
        let circle = CurveDef::Circle(CircleCurve {
            center: Point3::origin(),
            axis: Vec3::z(),
            radius: 5.0,
            ref_dir: Vec3::x(),
        });
        let bigger = circle.offset(2.0, &Vec3::z());
        // For CCW circle with axis=Z, offset by +2 should increase radius
        assert_relative_eq!(bigger.eval(0.0).x, 7.0, epsilon = 1e-10);
    }

    #[test]
    fn test_length_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(3.0, 4.0, 0.0),
        });
        assert_relative_eq!(seg.length(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn test_length_circle() {
        let circle = CurveDef::Circle(CircleCurve {
            center: Point3::origin(),
            axis: Vec3::z(),
            radius: 1.0,
            ref_dir: Vec3::x(),
        });
        assert_relative_eq!(circle.length(), TAU, epsilon = 1e-12);
    }

    #[test]
    fn test_length_arc() {
        let arc = CurveDef::CircularArc(CircularArcCurve {
            center: Point3::origin(),
            axis: Vec3::z(),
            radius: 2.0,
            ref_dir: Vec3::x(),
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        });
        assert_relative_eq!(arc.length(), std::f64::consts::PI, epsilon = 1e-12);
    }

    #[test]
    fn test_bounding_box_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(-1.0, -2.0, 0.0),
            end: Point3::new(3.0, 4.0, 5.0),
        });
        let bb = seg.bounding_box();
        assert!(bb.contains_point(&Point3::new(0.0, 0.0, 1.0)));
        assert!(!bb.contains_point(&Point3::new(10.0, 0.0, 0.0)));
    }

    #[test]
    fn test_extend_line() {
        let seg = CurveDef::LineSegment(LineSegment {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(10.0, 0.0, 0.0),
        });
        let extended = seg.extend(5.0, CurveEnd::End);
        assert_relative_eq!(extended.eval(1.0).x, 15.0, epsilon = 1e-10);
    }
}

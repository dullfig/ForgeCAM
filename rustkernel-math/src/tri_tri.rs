use crate::{Point3, Vec3};

/// Result of a triangle-triangle intersection test.
#[derive(Debug, Clone)]
pub struct TriTriResult {
    /// Whether the triangles intersect.
    pub intersects: bool,
    /// If intersecting, the segment endpoints on the intersection line.
    pub segment: Option<(Point3, Point3)>,
}

/// Möller triangle-triangle intersection.
///
/// Given two triangles (a0,a1,a2) and (b0,b1,b2), determines whether they
/// intersect. If they do, returns the intersection segment endpoints.
///
/// Algorithm:
/// 1. Compute plane of each triangle.
/// 2. Test vertices of each triangle against the other's plane.
/// 3. If all vertices on one side → no intersection.
/// 4. Compute line of plane-plane intersection.
/// 5. Project triangle edges onto the line, compute overlap interval.
/// 6. Return segment endpoints.
pub fn tri_tri_intersect(
    a0: &Point3, a1: &Point3, a2: &Point3,
    b0: &Point3, b1: &Point3, b2: &Point3,
) -> TriTriResult {
    let no_hit = TriTriResult { intersects: false, segment: None };

    // Plane of triangle B
    let nb = (b1 - b0).cross(&(b2 - b0));
    let nb_len = nb.norm();
    if nb_len < 1e-15 {
        return no_hit; // degenerate triangle B
    }
    let db = nb.dot(&b0.coords);

    // Signed distances of A vertices to plane of B
    let da0 = nb.dot(&a0.coords) - db;
    let da1 = nb.dot(&a1.coords) - db;
    let da2 = nb.dot(&a2.coords) - db;

    // Snap near-zero to zero
    let da0 = if da0.abs() < 1e-12 { 0.0 } else { da0 };
    let da1 = if da1.abs() < 1e-12 { 0.0 } else { da1 };
    let da2 = if da2.abs() < 1e-12 { 0.0 } else { da2 };

    // All on same side → no intersection
    if da0 > 0.0 && da1 > 0.0 && da2 > 0.0 { return no_hit; }
    if da0 < 0.0 && da1 < 0.0 && da2 < 0.0 { return no_hit; }

    // Plane of triangle A
    let na = (a1 - a0).cross(&(a2 - a0));
    let na_len = na.norm();
    if na_len < 1e-15 {
        return no_hit; // degenerate triangle A
    }
    let da = na.dot(&a0.coords);

    // Signed distances of B vertices to plane of A
    let db0 = na.dot(&b0.coords) - da;
    let db1 = na.dot(&b1.coords) - da;
    let db2 = na.dot(&b2.coords) - da;

    let db0 = if db0.abs() < 1e-12 { 0.0 } else { db0 };
    let db1 = if db1.abs() < 1e-12 { 0.0 } else { db1 };
    let db2 = if db2.abs() < 1e-12 { 0.0 } else { db2 };

    if db0 > 0.0 && db1 > 0.0 && db2 > 0.0 { return no_hit; }
    if db0 < 0.0 && db1 < 0.0 && db2 < 0.0 { return no_hit; }

    // Direction of the intersection line
    let line_dir = na.cross(&nb);
    let line_dir_len = line_dir.norm();
    if line_dir_len < 1e-15 {
        // Coplanar triangles — not handled (return no intersection)
        return no_hit;
    }
    let d = line_dir / line_dir_len;

    // Project vertices onto the intersection line direction to get scalar parameters.
    // We use the component with largest magnitude for numerical stability.
    let ax = d.x.abs();
    let ay = d.y.abs();
    let az = d.z.abs();
    let project = |p: &Point3| -> f64 {
        if ax >= ay && ax >= az { p.x * d.x.signum() }
        else if ay >= az { p.y * d.y.signum() }
        else { p.z * d.z.signum() }
    };

    // Compute intervals on the intersection line for each triangle.
    let interval_a = compute_interval(a0, a1, a2, da0, da1, da2, &project);
    let interval_b = compute_interval(b0, b1, b2, db0, db1, db2, &project);

    let (ia, ib) = match (interval_a, interval_b) {
        (Some(a), Some(b)) => (a, b),
        _ => return no_hit,
    };

    // Overlap of the two intervals
    let start = ia.0.max(ib.0);
    let end = ia.1.min(ib.1);

    if end - start < -1e-12 {
        return no_hit;
    }

    // Reconstruct 3D points on the intersection line.
    // Find a point on the line: solve for the intersection of the two planes.
    let line_origin = find_line_point(&na, da, &nb, db, &d);

    let p0 = line_origin + start * d;
    let p1 = line_origin + end * d;

    if (p1 - p0).norm() < 1e-14 {
        // Point contact
        return TriTriResult {
            intersects: true,
            segment: Some((p0, p0)),
        };
    }

    TriTriResult {
        intersects: true,
        segment: Some((p0, p1)),
    }
}

/// Compute the interval [t_min, t_max] where a triangle's edges cross the
/// opposite triangle's plane, projected onto the intersection line direction.
fn compute_interval(
    v0: &Point3, v1: &Point3, v2: &Point3,
    d0: f64, d1: f64, d2: f64,
    project: &dyn Fn(&Point3) -> f64,
) -> Option<(f64, f64)> {
    let mut ts = Vec::with_capacity(2);

    // For each edge, if the endpoints straddle the plane, compute the crossing point.
    maybe_crossing(v0, v1, d0, d1, project, &mut ts);
    maybe_crossing(v1, v2, d1, d2, project, &mut ts);
    maybe_crossing(v2, v0, d2, d0, project, &mut ts);

    // Also include any vertex exactly on the plane.
    if d0 == 0.0 { ts.push(project(v0)); }
    if d1 == 0.0 { ts.push(project(v1)); }
    if d2 == 0.0 { ts.push(project(v2)); }

    if ts.is_empty() {
        return None;
    }

    let mut t_min = ts[0];
    let mut t_max = ts[0];
    for &t in &ts[1..] {
        if t < t_min { t_min = t; }
        if t > t_max { t_max = t; }
    }
    Some((t_min, t_max))
}

/// If edge (va, vb) straddles the plane (d_a and d_b have opposite signs),
/// compute the crossing point and push its projection onto the intersection line.
fn maybe_crossing(
    va: &Point3, vb: &Point3,
    d_a: f64, d_b: f64,
    project: &dyn Fn(&Point3) -> f64,
    ts: &mut Vec<f64>,
) {
    if (d_a > 0.0 && d_b < 0.0) || (d_a < 0.0 && d_b > 0.0) {
        let t = d_a / (d_a - d_b);
        let crossing = Point3::from(va.coords * (1.0 - t) + vb.coords * t);
        ts.push(project(&crossing));
    }
}

/// Find a point on the intersection line of two planes.
/// Planes: na · x = da, nb · x = db.
/// Line direction: d = na × nb.
/// Solve for a point on the line by finding the point closest to origin.
fn find_line_point(na: &Vec3, da: f64, nb: &Vec3, db: f64, d: &Vec3) -> Point3 {
    // Use the formula: p = ((da * nb - db * na) × d) / |d|²
    let d_sq = d.dot(d);
    if d_sq < 1e-30 {
        return Point3::origin();
    }
    let p = (da * nb - db * na).cross(d) / d_sq;
    Point3::from(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perpendicular_crossing() {
        // Two triangles crossing each other at right angles
        let a0 = Point3::new(-1.0, 0.0, 0.0);
        let a1 = Point3::new(1.0, 0.0, 0.0);
        let a2 = Point3::new(0.0, 0.0, 1.0);

        let b0 = Point3::new(0.0, -1.0, 0.25);
        let b1 = Point3::new(0.0, 1.0, 0.25);
        let b2 = Point3::new(0.0, 0.0, 0.75);

        let result = tri_tri_intersect(&a0, &a1, &a2, &b0, &b1, &b2);
        assert!(result.intersects);
        assert!(result.segment.is_some());
        let (p0, p1) = result.segment.unwrap();
        // Intersection should be along x=0, in the z range [0.25, 0.75]
        assert!(p0.x.abs() < 1e-10);
        assert!(p1.x.abs() < 1e-10);
    }

    #[test]
    fn test_disjoint_triangles() {
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
    fn test_coplanar_triangles() {
        // Coplanar — we return no intersection (handled elsewhere)
        let a0 = Point3::new(0.0, 0.0, 0.0);
        let a1 = Point3::new(1.0, 0.0, 0.0);
        let a2 = Point3::new(0.0, 1.0, 0.0);

        let b0 = Point3::new(0.5, 0.5, 0.0);
        let b1 = Point3::new(1.5, 0.5, 0.0);
        let b2 = Point3::new(0.5, 1.5, 0.0);

        let result = tri_tri_intersect(&a0, &a1, &a2, &b0, &b1, &b2);
        assert!(!result.intersects);
    }

    #[test]
    fn test_edge_touching() {
        // Triangle A in XY plane, triangle B touches along shared edge
        let a0 = Point3::new(0.0, 0.0, 0.0);
        let a1 = Point3::new(1.0, 0.0, 0.0);
        let a2 = Point3::new(0.0, 1.0, 0.0);

        let b0 = Point3::new(0.0, 0.0, 0.0);
        let b1 = Point3::new(1.0, 0.0, 0.0);
        let b2 = Point3::new(0.5, 0.0, 1.0);

        let result = tri_tri_intersect(&a0, &a1, &a2, &b0, &b1, &b2);
        // Edge touching should register as intersecting
        assert!(result.intersects);
    }

    #[test]
    fn test_separated_by_plane() {
        // All vertices of B on one side of A's plane
        let a0 = Point3::new(0.0, 0.0, 0.0);
        let a1 = Point3::new(1.0, 0.0, 0.0);
        let a2 = Point3::new(0.0, 1.0, 0.0);

        let b0 = Point3::new(0.0, 0.0, 1.0);
        let b1 = Point3::new(1.0, 0.0, 2.0);
        let b2 = Point3::new(0.0, 1.0, 3.0);

        let result = tri_tri_intersect(&a0, &a1, &a2, &b0, &b1, &b2);
        assert!(!result.intersects);
    }
}

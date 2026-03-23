use rustkernel_math::polygon2d::{PointClassification, Polygon2D};
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::face_util::{face_boundary_points, polygon_centroid_and_normal, PlaneFrame};
use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::{FaceIdx, SolidIdx};

/// Classification of a point relative to a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointVsSolid {
    Inside,
    Outside,
    OnBoundary,
}

/// Classification of a face relative to another solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FacePosition {
    /// Entirely outside the other solid.
    Outside,
    /// Entirely inside the other solid.
    Inside,
    /// Coplanar with a face of the other solid, normals in same direction.
    OnSame,
    /// Coplanar with a face of the other solid, normals in opposite direction.
    OnOpposite,
}

/// Classify where a point lies relative to a solid using ray casting.
pub fn classify_point_vs_solid(
    point: &Point3,
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    solid: SolidIdx,
) -> PointVsSolid {
    let shell_idx = topo.solids.get(solid).outer_shell();
    let faces = &topo.shells.get(shell_idx).faces;

    // Try ray directions; if we hit a degenerate case, perturb.
    let directions = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.57735, 0.57735, 0.57735), // (1,1,1)/sqrt(3)
    ];

    for ray_dir in &directions {
        match ray_cast_classify(point, ray_dir, topo, geom, faces) {
            Some(result) => return result,
            None => continue, // degenerate — try another direction
        }
    }

    // All directions degenerate — very unlikely, default to outside.
    PointVsSolid::Outside
}

/// Compute ray-surface intersection. Returns (t, plane_origin_for_frame, plane_normal_for_frame).
/// Returns None if the ray misses or the surface type is not supported.
fn ray_surface_intersection(
    point: &Point3,
    ray_dir: &Vec3,
    kind: &SurfaceKind,
) -> Option<(f64, Point3, Vec3)> {
    match kind {
        SurfaceKind::Plane { origin, normal } => {
            let denom = ray_dir.dot(normal);
            if denom.abs() < 1e-12 {
                return None;
            }
            let t = (origin - point).dot(normal) / denom;
            Some((t, *origin, *normal))
        }
        SurfaceKind::Sphere { center, radius } => {
            // Quadratic: |P + t*d - C|^2 = r^2
            let r = radius.abs();
            let oc = point - center;
            let a = ray_dir.dot(ray_dir);
            let b = 2.0 * oc.dot(ray_dir);
            let c = oc.dot(&oc) - r * r;
            let disc = b * b - 4.0 * a * c;
            if disc < 0.0 {
                return None;
            }
            let sqrt_disc = disc.sqrt();
            let t1 = (-b - sqrt_disc) / (2.0 * a);
            let t2 = (-b + sqrt_disc) / (2.0 * a);
            // Pick the smallest positive t.
            let t = if t1 > 1e-10 { t1 } else if t2 > 1e-10 { t2 } else { return None };
            let hit = point + t * ray_dir;
            let normal = (hit - center).normalize();
            Some((t, hit, normal))
        }
        SurfaceKind::Cylinder { origin: co, axis, radius } => {
            // Quadratic in the component perpendicular to axis.
            let r = radius.abs();
            let a_dir = axis.normalize();
            let oc = point - co;
            let d_perp = ray_dir - ray_dir.dot(&a_dir) * a_dir;
            let oc_perp = oc - oc.dot(&a_dir) * a_dir;
            let a = d_perp.dot(&d_perp);
            let b = 2.0 * oc_perp.dot(&d_perp);
            let c = oc_perp.dot(&oc_perp) - r * r;
            let disc = b * b - 4.0 * a * c;
            if disc < 0.0 {
                return None;
            }
            let sqrt_disc = disc.sqrt();
            let t1 = (-b - sqrt_disc) / (2.0 * a);
            let t2 = (-b + sqrt_disc) / (2.0 * a);
            let t = if t1 > 1e-10 { t1 } else if t2 > 1e-10 { t2 } else { return None };
            let hit = point + t * ray_dir;
            let to_hit = hit - co;
            let normal = (to_hit - to_hit.dot(&a_dir) * a_dir).normalize();
            Some((t, hit, normal))
        }
        SurfaceKind::Cone { apex, axis, half_angle } => {
            // Quadratic for ray-cone intersection.
            let a_dir = axis.normalize();
            let ha = half_angle.abs();
            let cos2 = ha.cos() * ha.cos();
            let sin2 = ha.sin() * ha.sin();
            let oc = point - apex;

            let d_dot_a = ray_dir.dot(&a_dir);
            let oc_dot_a = oc.dot(&a_dir);

            let a = d_dot_a * d_dot_a - cos2 * ray_dir.dot(ray_dir) / (cos2 + sin2)
                + sin2 * d_dot_a * d_dot_a / (cos2 + sin2);
            // Simplified: use standard cone intersection formula.
            let a_coeff = d_dot_a * d_dot_a * cos2 - ray_dir.dot(ray_dir) * sin2
                + d_dot_a * d_dot_a * sin2;
            // Actually, let's use the correct formula:
            // For cone: (P·a)^2 * cos^2(α) = |P|^2 * sin^2(α) where P = hit - apex
            // Substituting P = oc + t*d:
            // (oc·a + t*d·a)^2 * cos^2 - (|oc + t*d|^2) * sin^2 = 0
            // Expanding:
            let qa = d_dot_a * d_dot_a * cos2 - (ray_dir.dot(ray_dir)) * sin2;
            let qb = 2.0 * (d_dot_a * oc_dot_a * cos2 - ray_dir.dot(&oc) * sin2);
            let qc = oc_dot_a * oc_dot_a * cos2 - oc.dot(&oc) * sin2;

            let disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0 || qa.abs() < 1e-20 {
                return None;
            }
            let sqrt_disc = disc.sqrt();
            let t1 = (-qb - sqrt_disc) / (2.0 * qa);
            let t2 = (-qb + sqrt_disc) / (2.0 * qa);
            // Only hits on the positive half of the cone (same side as axis direction).
            let valid = |t: f64| -> bool {
                if t < 1e-10 { return false; }
                let p = oc + t * ray_dir;
                p.dot(&a_dir) > 0.0
            };
            let t = if valid(t1) { t1 } else if valid(t2) { t2 } else { return None };
            let hit = point + t * ray_dir;
            let p = hit - apex;
            let along = p.dot(&a_dir);
            let radial = p - along * a_dir;
            let normal = if radial.norm() > 1e-12 {
                (ha.cos() * radial.normalize() - ha.sin() * a_dir).normalize()
            } else {
                a_dir
            };
            Some((t, hit, normal))
        }
        SurfaceKind::Torus { .. } => {
            // Quartic — defer (return None, skip face).
            None
        }
        SurfaceKind::Nurbs | SurfaceKind::Unknown => None,
    }
}

/// Perform a single ray cast and count crossings.
/// Returns None if the ray hits a degenerate case (edge/vertex).
fn ray_cast_classify(
    point: &Point3,
    ray_dir: &Vec3,
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    faces: &[FaceIdx],
) -> Option<PointVsSolid> {
    let mut crossings = 0i32;

    for &face_idx in faces {
        let surface_id = topo.faces.get(face_idx).surface_id;
        let kind = geom.surface_kind(surface_id);

        // Compute ray-surface intersection parameter t.
        let (t, plane_origin, plane_normal) = match ray_surface_intersection(point, ray_dir, &kind) {
            Some(result) => result,
            None => {
                // Fallback for NURBS (and torus): use best-fit plane from boundary polygon.
                if matches!(kind, SurfaceKind::Nurbs | SurfaceKind::Torus { .. }) {
                    let boundary = face_boundary_points(topo, geom, face_idx);
                    if boundary.len() < 3 { continue; }
                    let (origin, normal) = polygon_centroid_and_normal(&boundary);
                    let denom = ray_dir.dot(&normal);
                    if denom.abs() < 1e-12 { continue; }
                    let t = (origin - point).dot(&normal) / denom;
                    (t, origin, normal)
                } else {
                    continue;
                }
            }
        };

        if t < 1e-10 {
            continue;
        }

        // Compute hit point.
        let hit = point + t * ray_dir;

        // Check if hit point is inside the face's polygon.
        let frame = PlaneFrame::from_normal(plane_origin, plane_normal);
        let boundary_3d = face_boundary_points(topo, geom, face_idx);
        let poly = Polygon2D {
            vertices: boundary_3d.iter().map(|p| frame.project_to_2d(p)).collect(),
        };
        let hit_2d = frame.project_to_2d(&hit);

        match poly.classify_point(hit_2d, 1e-8) {
            PointClassification::Inside => {
                crossings += 1;
            }
            PointClassification::OnBoundary => {
                // Hit an edge — this direction is degenerate.
                return None;
            }
            PointClassification::Outside => {
                // Missed the face.
            }
        }
    }

    if crossings % 2 == 1 {
        Some(PointVsSolid::Inside)
    } else {
        Some(PointVsSolid::Outside)
    }
}

/// Classify a face of one solid relative to another solid.
///
/// If the face has intersection segments, the caller should split the face
/// first and classify each sub-face independently.
///
/// For unsplit faces: pick the face centroid and ray-cast against the other solid.
/// For coplanar faces: compare normals.
pub fn classify_face(
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    face: FaceIdx,
    other_solid: SolidIdx,
) -> FacePosition {
    let surface_id = topo.faces.get(face).surface_id;
    let kind = geom.surface_kind(surface_id);
    let (face_origin, face_normal) = match &kind {
        SurfaceKind::Plane { origin, normal } => (*origin, *normal),
        _ => {
            // For curved faces, use centroid and surface normal at centroid.
            let boundary = face_boundary_points(topo, geom, face);
            let centroid = compute_centroid(&boundary);
            let (u, v) = geom.surface_inverse_uv(surface_id, &centroid);
            let normal = geom.surface_normal(surface_id, u, v);
            (centroid, normal)
        }
    };

    // Check for coplanarity with any face of the other solid (only for planar faces).
    let is_planar = matches!(kind, SurfaceKind::Plane { .. });
    let other_shell = topo.solids.get(other_solid).outer_shell();
    let other_faces = &topo.shells.get(other_shell).faces;

    if is_planar {
    for &other_face_idx in other_faces {
        let other_sid = topo.faces.get(other_face_idx).surface_id;
        let other_kind = geom.surface_kind(other_sid);
        let (other_origin, other_normal) = match other_kind {
            SurfaceKind::Plane { origin, normal } => (origin, normal),
            _ => continue,
        };

        // Check if normals are parallel.
        let cross = face_normal.cross(&other_normal);
        if cross.norm() < 1e-8 {
            // Parallel normals — check if same plane (signed distance).
            let signed_dist = face_normal.dot(&(other_origin - face_origin));
            if signed_dist.abs() < 1e-8 {
                // Same plane. Check if the faces actually overlap by testing if
                // the centroid of one face lies inside the other face's polygon.
                let frame = PlaneFrame::from_normal(face_origin, face_normal);
                let my_boundary = face_boundary_points(topo, geom, face);
                let other_boundary = face_boundary_points(topo, geom, other_face_idx);

                let my_centroid = compute_centroid(&my_boundary);
                let other_poly = Polygon2D {
                    vertices: other_boundary.iter().map(|p| frame.project_to_2d(p)).collect(),
                };
                let my_centroid_2d = frame.project_to_2d(&my_centroid);

                let overlaps = matches!(
                    other_poly.classify_point(my_centroid_2d, 1e-8),
                    PointClassification::Inside | PointClassification::OnBoundary
                );

                if overlaps {
                    if face_normal.dot(&other_normal) > 0.0 {
                        return FacePosition::OnSame;
                    } else {
                        return FacePosition::OnOpposite;
                    }
                }
            }
        }
    }
    } // end if is_planar

    // Not coplanar — classify by centroid.
    let boundary = face_boundary_points(topo, geom, face);
    let centroid = compute_centroid(&boundary);

    // Offset centroid slightly inward (along negative face normal) to avoid
    // landing exactly on the face surface of the other solid.
    let test_point = centroid - 1e-6 * face_normal;

    match classify_point_vs_solid(&test_point, topo, geom, other_solid) {
        PointVsSolid::Inside => FacePosition::Inside,
        PointVsSolid::Outside => FacePosition::Outside,
        PointVsSolid::OnBoundary => {
            // Nudge didn't help — try the raw centroid.
            match classify_point_vs_solid(&centroid, topo, geom, other_solid) {
                PointVsSolid::Inside => FacePosition::Inside,
                _ => FacePosition::Outside,
            }
        }
    }
}

/// Compute the centroid of a polygon (average of vertices).
fn compute_centroid(points: &[Point3]) -> Point3 {
    let n = points.len() as f64;
    let sum = points.iter().fold(Vec3::zeros(), |acc, p| acc + p.coords);
    Point3::from(sum / n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_geom::AnalyticalGeomStore;
    use rustkernel_topology::store::TopoStore;

    #[test]
    fn test_point_inside_box() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let result = classify_point_vs_solid(&Point3::new(0.0, 0.0, 0.0), &topo, &geom, solid);
        assert_eq!(result, PointVsSolid::Inside);
    }

    #[test]
    fn test_point_outside_box() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let result = classify_point_vs_solid(&Point3::new(10.0, 0.0, 0.0), &topo, &geom, solid);
        assert_eq!(result, PointVsSolid::Outside);
    }

    #[test]
    fn test_point_inside_offset_box() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::new(5.0, 5.0, 5.0), 2.0, 2.0, 2.0);
        let result = classify_point_vs_solid(&Point3::new(5.0, 5.0, 5.0), &topo, &geom, solid);
        assert_eq!(result, PointVsSolid::Inside);
    }

    #[test]
    fn test_classify_face_outside() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(5.0, 0.0, 0.0), 2.0, 2.0, 2.0);

        // All faces of A should be outside B (no overlap).
        let shell_a = topo.solids.get(a).outer_shell();
        for &face in &topo.shells.get(shell_a).faces.clone() {
            let pos = classify_face(&topo, &geom, face, b);
            assert_eq!(pos, FacePosition::Outside, "Face should be outside");
        }
    }

    #[test]
    fn test_classify_face_inside() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let big = make_box_into(&mut topo, &mut geom, Point3::origin(), 10.0, 10.0, 10.0);
        let small = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // All faces of small should be inside big.
        let shell_small = topo.solids.get(small).outer_shell();
        for &face in &topo.shells.get(shell_small).faces.clone() {
            let pos = classify_face(&topo, &geom, face, big);
            assert_eq!(pos, FacePosition::Inside, "Small box face should be inside big box");
        }
    }
}

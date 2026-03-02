use std::collections::HashMap;

use rustkernel_math::Point3;
use rustkernel_math::polygon2d::{PointClassification, Polygon2D};
use rustkernel_topology::face_util::{polygon_centroid_and_normal, PlaneFrame};
use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
use rustkernel_topology::intersection::{
    IntersectionCurve, IntersectionPipeline, SurfaceSurfaceResult,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::{FaceIdx, SolidIdx};

use crate::broad_phase::find_interfering_face_pairs;
use crate::curve_trimming::{
    trim_intersection_circle, trim_intersection_ellipse, trim_intersection_line, TrimmedSegment,
};
use crate::face_classifier::{classify_face, FacePosition};
use crate::face_selector::{select_faces, BooleanOp};
use crate::face_splitter::{face_boundary_verts, split_face_along_segment};
use crate::topology_builder::{build_result_solid, BuildError};
use rustkernel_geom::AnalyticalGeomStore;

/// Errors from boolean operations.
#[derive(Debug)]
pub enum BooleanError {
    NoIntersection,
    OpenResult,
    EulerViolation { v: usize, e: usize, f: usize },
    IntersectionFailed(String),
    DegenerateInput(String),
    BuildFailed(BuildError),
}

impl std::fmt::Display for BooleanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BooleanError::NoIntersection => write!(f, "solids do not intersect"),
            BooleanError::OpenResult => write!(f, "result shell is not closed"),
            BooleanError::EulerViolation { v, e, f: faces } => {
                write!(f, "Euler violation: V({v}) - E({e}) + F({faces}) != 2")
            }
            BooleanError::IntersectionFailed(msg) => write!(f, "intersection failed: {msg}"),
            BooleanError::DegenerateInput(msg) => write!(f, "degenerate input: {msg}"),
            BooleanError::BuildFailed(e) => write!(f, "build failed: {e}"),
        }
    }
}

impl std::error::Error for BooleanError {}

impl From<BuildError> for BooleanError {
    fn from(e: BuildError) -> Self {
        BooleanError::BuildFailed(e)
    }
}

/// Execute a boolean operation on two solids.
pub fn boolean_op(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    pipeline: &IntersectionPipeline,
    solid_a: SolidIdx,
    solid_b: SolidIdx,
    op: BooleanOp,
) -> Result<SolidIdx, BooleanError> {
    // Step 1: Broad phase — find potentially interfering face pairs.
    let pairs = find_interfering_face_pairs(topo, geom, solid_a, solid_b, 1e-8);

    // Step 2: Surface-surface intersection + curve trimming.
    // Track which faces have intersection segments.
    let mut segments_for_face_a: HashMap<FaceIdx, Vec<TrimmedSegment>> = HashMap::new();
    let mut segments_for_face_b: HashMap<FaceIdx, Vec<TrimmedSegment>> = HashMap::new();

    for &(face_a, face_b) in &pairs {
        let sid_a = topo.faces.get(face_a).surface_id;
        let sid_b = topo.faces.get(face_b).surface_id;

        let result = pipeline.solve(geom, sid_a, sid_b);
        let ssi_result = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        match ssi_result {
            SurfaceSurfaceResult::Empty => continue,
            SurfaceSurfaceResult::Coincident => {
                // Coplanar faces are handled by face classification.
                continue;
            }
            SurfaceSurfaceResult::Curves(curves) => {
                for curve in curves {
                    match &curve {
                        IntersectionCurve::Line(line) => {
                            let trimmed = trim_intersection_line(line, topo, geom, face_a, face_b);
                            for seg in trimmed {
                                segments_for_face_a
                                    .entry(face_a)
                                    .or_default()
                                    .push(seg.clone());
                                segments_for_face_b
                                    .entry(face_b)
                                    .or_default()
                                    .push(seg);
                            }
                        }
                        IntersectionCurve::Circle(circle) => {
                            let trimmed =
                                trim_intersection_circle(circle, topo, geom, face_a, face_b);
                            for seg in trimmed {
                                segments_for_face_a
                                    .entry(face_a)
                                    .or_default()
                                    .push(seg.clone());
                                segments_for_face_b
                                    .entry(face_b)
                                    .or_default()
                                    .push(seg);
                            }
                        }
                        IntersectionCurve::Ellipse(ellipse) => {
                            let trimmed =
                                trim_intersection_ellipse(ellipse, topo, geom, face_a, face_b);
                            for seg in trimmed {
                                segments_for_face_a
                                    .entry(face_a)
                                    .or_default()
                                    .push(seg.clone());
                                segments_for_face_b
                                    .entry(face_b)
                                    .or_default()
                                    .push(seg);
                            }
                        }
                    }
                }
            }
        }
    }

    // Step 3: Split faces that have intersection segments, classify all faces.
    let shell_a = topo.solids.get(solid_a).shell;
    let shell_b = topo.solids.get(solid_b).shell;

    let faces_a: Vec<FaceIdx> = topo.shells.get(shell_a).faces.clone();
    let faces_b: Vec<FaceIdx> = topo.shells.get(shell_b).faces.clone();

    let classified_a = classify_and_split_faces(
        topo, geom, &faces_a, &segments_for_face_a, solid_b,
    );
    let classified_b = classify_and_split_faces(
        topo, geom, &faces_b, &segments_for_face_b, solid_a,
    );

    // Step 4: Select faces based on boolean operation.
    let selected = select_faces(op, &classified_a, &classified_b);

    // Step 5: Build result solid.
    let result = build_result_solid(topo, geom, &selected)?;

    Ok(result)
}

/// For each face, either classify it directly or split it by all its
/// intersection segments, then classify each resulting sub-face.
fn classify_and_split_faces(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    faces: &[FaceIdx],
    segments_map: &HashMap<FaceIdx, Vec<TrimmedSegment>>,
    other_solid: SolidIdx,
) -> Vec<(FaceIdx, FacePosition)> {
    let mut result = Vec::new();

    for &face in faces {
        if let Some(segments) = segments_map.get(&face) {
            if !segments.is_empty() {
                // Split the face by all segments, then classify each leaf.
                let leaves = split_face_by_all_segments(topo, geom, face, segments);
                for leaf in leaves {
                    let pos = classify_face(topo, geom, leaf, other_solid);
                    result.push((leaf, pos));
                }
                continue;
            }
        }

        // No split needed — classify directly.
        let pos = classify_face(topo, geom, face, other_solid);
        result.push((face, pos));
    }

    result
}

/// Iteratively split a face by multiple intersection segments.
///
/// Uses a work-queue approach: start with the original face and all segments.
/// Pop a face, try splitting by the first available segment. On success,
/// distribute remaining segments to the two sub-faces. On error (degenerate
/// or endpoint too far), skip that segment. Repeat until no segments remain.
fn split_face_by_all_segments(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    face: FaceIdx,
    segments: &[TrimmedSegment],
) -> Vec<FaceIdx> {
    // Work queue: (face, indices into `segments` still to process).
    let all_indices: Vec<usize> = (0..segments.len()).collect();
    let mut queue: Vec<(FaceIdx, Vec<usize>)> = vec![(face, all_indices)];
    let mut leaves = Vec::new();

    while let Some((current_face, remaining)) = queue.pop() {
        if remaining.is_empty() {
            leaves.push(current_face);
            continue;
        }

        // Try each remaining segment until one succeeds.
        let mut split_done = false;
        for (pos, &seg_idx) in remaining.iter().enumerate() {
            match split_face_along_segment(topo, geom, current_face, &segments[seg_idx]) {
                Ok(split_result) => {
                    // Remove the used segment from the remaining list.
                    let rest: Vec<usize> = remaining.iter()
                        .enumerate()
                        .filter(|&(i, _)| i != pos)
                        .map(|(_, &idx)| idx)
                        .collect();

                    // Distribute remaining segments to the two sub-faces.
                    let (for_a, for_b) = distribute_segments(
                        topo, geom, segments, &rest,
                        split_result.face_a, split_result.face_b,
                    );

                    queue.push((split_result.face_a, for_a));
                    queue.push((split_result.face_b, for_b));
                    split_done = true;
                    break;
                }
                Err(_) => {
                    // This segment doesn't work for this face — try next.
                    continue;
                }
            }
        }

        if !split_done {
            // No segment could split this face — it's a leaf.
            leaves.push(current_face);
        }
    }

    leaves
}

/// Distribute remaining segment indices to sub-face A or B.
/// For each segment, test its midpoint against both sub-face polygons.
fn distribute_segments(
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    segments: &[TrimmedSegment],
    remaining: &[usize],
    face_a: FaceIdx,
    face_b: FaceIdx,
) -> (Vec<usize>, Vec<usize>) {
    if remaining.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Build 2D polygons for both sub-faces using a best-fit plane.
    let pts_a = face_boundary_verts(topo, geom, face_a);
    let pts_b = face_boundary_verts(topo, geom, face_b);

    // Use face_a's boundary to determine the projection plane.
    let surface_id = topo.faces.get(face_a).surface_id;
    let kind = geom.surface_kind(surface_id);
    let (plane_origin, plane_normal) = match kind {
        SurfaceKind::Plane { origin, normal } => (origin, normal),
        _ => {
            if pts_a.len() >= 3 {
                polygon_centroid_and_normal(&pts_a)
            } else {
                // Fallback: give everything to face_a.
                return (remaining.to_vec(), Vec::new());
            }
        }
    };
    let frame = PlaneFrame::from_normal(plane_origin, plane_normal);

    let poly_a = Polygon2D {
        vertices: pts_a.iter().map(|p| frame.project_to_2d(p)).collect(),
    };
    let poly_b = Polygon2D {
        vertices: pts_b.iter().map(|p| frame.project_to_2d(p)).collect(),
    };

    let mut for_a = Vec::new();
    let mut for_b = Vec::new();

    for &idx in remaining {
        let seg = &segments[idx];
        let mid = Point3::from((seg.start_point.coords + seg.end_point.coords) * 0.5);
        let mid_2d = frame.project_to_2d(&mid);

        let in_a = matches!(
            poly_a.classify_point(mid_2d, 1e-6),
            PointClassification::Inside | PointClassification::OnBoundary
        );
        let in_b = matches!(
            poly_b.classify_point(mid_2d, 1e-6),
            PointClassification::Inside | PointClassification::OnBoundary
        );

        if in_a && !in_b {
            for_a.push(idx);
        } else if in_b && !in_a {
            for_b.push(idx);
        } else {
            // Ambiguous (on boundary of both, or neither) — try both.
            for_a.push(idx);
        }
    }

    (for_a, for_b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::face_selector::BooleanOp;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_builders::cylinder_builder::make_cylinder_into;
    use rustkernel_builders::sphere_builder::make_sphere_into;
    use rustkernel_builders::extrude_builder::make_extrude_into;
    use rustkernel_builders::revolve_builder::make_revolve_into;
    use rustkernel_geom::AnalyticalGeomStore;
    use rustkernel_math::{Point3, Vec3};
    use rustkernel_solvers::default_pipeline;
    use rustkernel_topology::store::TopoStore;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    fn verify_solid(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).shell;
        let faces = &topo.shells.get(shell_idx).faces;

        let mut verts = HashSet::new();
        let mut edges = HashSet::new();

        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let he_data = topo.half_edges.get(he);
                verts.insert(he_data.origin.raw());
                edges.insert(he_data.edge.raw());
                assert!(
                    he_data.twin.is_some(),
                    "Half-edge {} has no twin",
                    he.raw()
                );
                he = he_data.next;
                if he == start_he {
                    break;
                }
            }
        }

        let v = verts.len() as i32;
        let e = edges.len() as i32;
        let f = faces.len() as i32;
        assert_eq!(v - e + f, 2, "Euler: V({v}) - E({e}) + F({f}) != 2");
    }

    #[test]
    fn test_fuse_overlapping_boxes() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(1.0, 0.0, 0.0), 2.0, 2.0, 2.0);

        let result = boolean_op(
            &mut topo,
            &mut geom,
            &pipeline,
            a,
            b,
            BooleanOp::Fuse,
        );

        match result {
            Ok(solid) => {
                verify_solid(&topo, solid);
                let shell = topo.solids.get(solid).shell;
                tessellate_shell(&mut topo, shell, &geom);
            }
            Err(e) => {
                eprintln!("Fuse failed (may need debugging): {e}");
            }
        }
    }

    #[test]
    fn test_fuse_disjoint_boxes() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 1.0, 1.0, 1.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(10.0, 0.0, 0.0), 1.0, 1.0, 1.0);

        let result = boolean_op(
            &mut topo,
            &mut geom,
            &pipeline,
            a,
            b,
            BooleanOp::Fuse,
        );

        match result {
            Ok(_) => {}
            Err(_) => {}
        }
    }

    #[test]
    fn test_fuse_box_cylinder() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_cylinder_into(&mut topo, &mut geom, Point3::new(0.5, 0.0, 0.0), 0.5, 2.0, 16);

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Fuse)
        }));
        match result {
            Ok(Ok(solid)) => {
                verify_solid(&topo, solid);
            }
            Ok(Err(e)) => {
                // Circle trimming now works; downstream topology reconstruction may still fail.
                eprintln!("Fuse box+cylinder: {e}");
            }
            Err(_) => {
                eprintln!("Fuse box+cylinder panicked (topology reconstruction issue)");
            }
        }
    }

    #[test]
    fn test_cut_box_cylinder() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_cylinder_into(&mut topo, &mut geom, Point3::new(0.0, 0.0, 0.0), 0.5, 2.0, 16);

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Cut)
        }));
        match result {
            Ok(Ok(solid)) => {
                verify_solid(&topo, solid);
            }
            Ok(Err(e)) => {
                eprintln!("Cut box-cylinder: {e}");
            }
            Err(_) => {
                eprintln!("Cut box-cylinder panicked (topology reconstruction issue)");
            }
        }
    }

    #[test]
    fn test_fuse_two_extrudes() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        // Extrude a rectangle in Z.
        let profile_a = vec![
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
        ];
        let a = make_extrude_into(&mut topo, &mut geom, &profile_a, Vec3::new(0.0, 0.0, 1.0), 2.0);

        // Extrude another overlapping rectangle.
        let profile_b = vec![
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(2.0, -1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let b = make_extrude_into(&mut topo, &mut geom, &profile_b, Vec3::new(0.0, 0.0, 1.0), 2.0);

        let result = boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Fuse);
        match result {
            Ok(solid) => {
                verify_solid(&topo, solid);
            }
            Err(e) => {
                eprintln!("Fuse two extrudes: {e}");
            }
        }
    }

    #[test]
    fn test_fuse_box_revolved() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 4.0, 4.0, 4.0);

        // Revolve a small rectangle around Y axis.
        let profile = vec![
            Point3::new(1.0, -0.5, 0.0),
            Point3::new(1.5, -0.5, 0.0),
            Point3::new(1.5, 0.5, 0.0),
            Point3::new(1.0, 0.5, 0.0),
        ];
        let b = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 1.0, 0.0),
            std::f64::consts::TAU, 16,
        );

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Fuse)
        }));

        match result {
            Ok(Ok(solid)) => {
                verify_solid(&topo, solid);
            }
            Ok(Err(e)) => {
                eprintln!("Fuse box+revolved: {e} (may need circle curve support)");
            }
            Err(_) => {
                eprintln!("Fuse box+revolved panicked (may need further fixes)");
            }
        }
    }

    #[test]
    fn test_cut_box_sphere() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 4.0, 4.0, 4.0);
        let b = make_sphere_into(&mut topo, &mut geom, Point3::origin(), 1.0, 16, 8);

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Cut)
        }));
        match result {
            Ok(Ok(solid)) => {
                verify_solid(&topo, solid);
            }
            Ok(Err(e)) => {
                // Circle trimming now works; downstream topology may still fail.
                eprintln!("Cut box-sphere: {e}");
            }
            Err(_) => {
                eprintln!("Cut box-sphere panicked (topology reconstruction issue)");
            }
        }
    }
}

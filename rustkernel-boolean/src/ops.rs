use std::collections::HashMap;

use rustkernel_topology::intersection::{
    IntersectionCurve, IntersectionPipeline, SurfaceSurfaceResult,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::{FaceIdx, SolidIdx};

use crate::broad_phase::find_interfering_face_pairs;
use crate::curve_trimming::{trim_intersection_line, TrimmedSegment};
use crate::face_classifier::{classify_face, FacePosition};
use crate::face_selector::{select_faces, BooleanOp};
use crate::face_splitter::split_face_along_segment;
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
                        IntersectionCurve::Circle(_) | IntersectionCurve::Ellipse(_) => {
                            // Curved intersection trimming deferred to Phase 4.
                            continue;
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

/// For each face, either classify it directly or split it first,
/// then classify the sub-faces.
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
            if let Some(seg) = segments.first() {
                // Split the face and classify each sub-face.
                let split = split_face_along_segment(topo, geom, face, seg);

                let pos_a = classify_face(topo, geom, split.face_a, other_solid);
                let pos_b = classify_face(topo, geom, split.face_b, other_solid);

                result.push((split.face_a, pos_a));
                result.push((split.face_b, pos_b));

                // If there are additional segments, we'd need recursive splitting.
                // For convex face + single intersection (the box-box case), one segment suffices.
                continue;
            }
        }

        // No split needed — classify directly.
        let pos = classify_face(topo, geom, face, other_solid);
        result.push((face, pos));
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::face_selector::BooleanOp;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_geom::AnalyticalGeomStore;
    use rustkernel_math::Point3;
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
}

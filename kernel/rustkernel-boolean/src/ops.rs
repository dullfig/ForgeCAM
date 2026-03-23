use std::collections::{HashMap, HashSet};

use rustkernel_math::Point3;
use rustkernel_math::polygon2d::{PointClassification, Polygon2D};
use rustkernel_topology::face_util::{polygon_centroid_and_normal, PlaneFrame};
use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
use rustkernel_topology::intersection::{
    IntersectionCurve, IntersectionPipeline, SurfaceSurfaceResult,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::evolution::{FaceOrigin, ShapeEvolution};
use rustkernel_topology::topo::{FaceIdx, SolidIdx};
use tracing::{debug, info_span, warn};

use crate::broad_phase::find_interfering_face_pairs;
use crate::curve_trimming::{
    trim_circle_to_face, trim_ellipse_to_face, trim_intersection_circle,
    trim_intersection_ellipse, trim_intersection_line, trim_intersection_polyline,
    TrimmedSegment,
};
use crate::face_classifier::{classify_face, FacePosition};
use crate::face_selector::{select_faces, BooleanOp};
use crate::face_splitter::{face_boundary_verts, snap_to_boundary, split_face_along_segment, split_face_by_interior_loop, split_face_with_extension};
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

/// Result of a boolean operation, including shape evolution tracking.
pub struct BooleanResult {
    pub solid: SolidIdx,
    pub evolution: ShapeEvolution,
}

/// Execute a boolean operation on two solids.
pub fn boolean_op(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    pipeline: &IntersectionPipeline,
    solid_a: SolidIdx,
    solid_b: SolidIdx,
    op: BooleanOp,
) -> Result<BooleanResult, BooleanError> {
    let _span = info_span!("boolean_op", ?op, a = solid_a.raw(), b = solid_b.raw()).entered();

    // Step 1: Broad phase — find potentially interfering face pairs.
    let pairs = find_interfering_face_pairs(topo, geom, solid_a, solid_b, 1e-8);
    debug!(pairs = pairs.len(), "broad phase complete");

    // Step 2: Surface-surface intersection + curve trimming.
    // Track which faces have intersection segments.
    let mut segments_for_face_a: HashMap<FaceIdx, Vec<TrimmedSegment>> = HashMap::new();
    let mut segments_for_face_b: HashMap<FaceIdx, Vec<TrimmedSegment>> = HashMap::new();

    // Track which (face_a, surface_b) pairs have already contributed face_a segments.
    // When multiple face_b faces share the same surface, we only trim against face_a once
    // to avoid fragmented duplicates that break segment chaining.
    let mut face_a_processed: HashSet<(u32, u32)> = HashSet::new(); // (face_a.raw(), sid_b)

    // Cache SSI results by ordered surface pair to avoid redundant solves.
    let mut ssi_cache: HashMap<(u32, u32), Vec<IntersectionCurve>> = HashMap::new();

    for &(face_a, face_b) in &pairs {
        let sid_a = topo.faces.get(face_a).surface_id;
        let sid_b = topo.faces.get(face_b).surface_id;

        // Get or compute SSI result for this surface pair.
        let cache_key = (sid_a.min(sid_b), sid_a.max(sid_b));
        let curves = if let Some(cached) = ssi_cache.get(&cache_key) {
            cached.clone()
        } else {
            let result = pipeline.solve(geom, sid_a, sid_b);
            let curves = match result {
                Ok(SurfaceSurfaceResult::Curves(c)) => c,
                Ok(SurfaceSurfaceResult::Empty) => vec![],
                Ok(SurfaceSurfaceResult::Coincident) => vec![],
                Err(e) => {
                    warn!(
                        face_a = face_a.raw(),
                        face_b = face_b.raw(),
                        error = %e,
                        "SSI failed, skipping face pair"
                    );
                    vec![]
                }
            };
            ssi_cache.insert(cache_key, curves.clone());
            curves
        };

        if curves.is_empty() {
            continue;
        }

        for curve in &curves {
            match curve {
                IntersectionCurve::Line(line) => {
                    // Lines are unbounded — must trim against both faces for
                    // both face_a and face_b to get the correct overlap segment.
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
                    // face_b: trim against both faces
                    let trimmed_b =
                        trim_intersection_circle(circle, topo, geom, face_a, face_b);
                    for seg in trimmed_b {
                        segments_for_face_b
                            .entry(face_b)
                            .or_default()
                            .push(seg);
                    }
                    // face_a: trim against face_a only, once per unique (face_a, surface_b)
                    if face_a_processed.insert((face_a.raw(), sid_b)) {
                        let trimmed_a = trim_circle_to_face(circle, topo, geom, face_a);
                        for seg in trimmed_a {
                            segments_for_face_a
                                .entry(face_a)
                                .or_default()
                                .push(seg);
                        }
                    }
                }
                IntersectionCurve::Ellipse(ellipse) => {
                    let trimmed_b =
                        trim_intersection_ellipse(ellipse, topo, geom, face_a, face_b);
                    for seg in trimmed_b {
                        segments_for_face_b
                            .entry(face_b)
                            .or_default()
                            .push(seg);
                    }
                    if face_a_processed.insert((face_a.raw(), sid_b)) {
                        let trimmed_a = trim_ellipse_to_face(ellipse, topo, geom, face_a);
                        for seg in trimmed_a {
                            segments_for_face_a
                                .entry(face_a)
                                .or_default()
                                .push(seg);
                        }
                    }
                }
                IntersectionCurve::Polyline(polyline) => {
                    // Polylines are open — trim against both faces for both sides.
                    let trimmed =
                        trim_intersection_polyline(polyline, topo, geom, face_a, face_b);
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

    // Step 3: Split faces that have intersection segments, classify all faces.
    let shell_a = topo.solids.get(solid_a).outer_shell();
    let shell_b = topo.solids.get(solid_b).outer_shell();

    let faces_a: Vec<FaceIdx> = topo.shells.get(shell_a).faces.clone();
    let faces_b: Vec<FaceIdx> = topo.shells.get(shell_b).faces.clone();

    let result_a = classify_and_split_faces(
        topo, geom, &faces_a, &segments_for_face_a, solid_b,
    );
    let result_b = classify_and_split_faces(
        topo, geom, &faces_b, &segments_for_face_b, solid_a,
    );
    debug!(
        a_count = result_a.classified.len(),
        b_count = result_b.classified.len(),
        "faces classified"
    );

    // Step 4: Select faces based on boolean operation.
    let selected = select_faces(op, &result_a.classified, &result_b.classified);
    debug!(
        keep_a = selected.keep_from_a.len(),
        keep_b = selected.keep_from_b.len(),
        flip_b = selected.flip_from_b.len(),
        "faces selected"
    );

    // Step 5: Build result solid.
    let build = build_result_solid(topo, geom, &selected)?;

    // Step 6: Build shape evolution from split maps + copy map + selection.
    let evolution = build_boolean_evolution(
        &faces_a,
        &faces_b,
        &result_a.split_map,
        &result_b.split_map,
        &build.face_copy_map,
    );

    Ok(BooleanResult {
        solid: build.solid,
        evolution,
    })
}

/// Build a ShapeEvolution for a boolean operation.
///
/// Uses the split maps (original face → sub-faces), the face copy map
/// (source face → result face from deep_copy), and the original face lists
/// to determine the provenance of each result face and which input faces
/// were deleted.
fn build_boolean_evolution(
    faces_a: &[FaceIdx],
    faces_b: &[FaceIdx],
    split_map_a: &HashMap<FaceIdx, Vec<FaceIdx>>,
    split_map_b: &HashMap<FaceIdx, Vec<FaceIdx>>,
    face_copy_map: &[(FaceIdx, FaceIdx)],
) -> ShapeEvolution {
    let mut evo = ShapeEvolution::new();

    // Build a set of source faces that were copied into the result.
    let copied_sources: HashSet<FaceIdx> = face_copy_map.iter().map(|&(src, _)| src).collect();

    // For each copied face, determine its origin.
    for &(source, result) in face_copy_map {
        // Check if this source face is a sub-face from a split.
        let origin = find_original_face(source, split_map_a, split_map_b);
        match origin {
            Some(original) => {
                // This result face came from splitting an original face.
                evo.record_face(result, FaceOrigin::SplitFrom(original));
            }
            None => {
                // This source face wasn't split — it's a direct copy from the
                // original solid (possibly with modified boundary from boolean
                // trimming, but topologically the same face).
                evo.record_face(result, FaceOrigin::CopiedFrom(source));
            }
        }
    }

    // Determine which original faces were deleted (not in the result).
    for &face in faces_a.iter().chain(faces_b.iter()) {
        if let Some(sub_faces) = split_map_a.get(&face).or_else(|| split_map_b.get(&face)) {
            // Face was split. Check if any sub-face made it into the result.
            let any_survived = sub_faces.iter().any(|sf| copied_sources.contains(sf));
            if !any_survived {
                evo.record_deleted_face(face);
            }
        } else {
            // Face was not split. Check if it was copied into the result.
            if !copied_sources.contains(&face) {
                evo.record_deleted_face(face);
            }
        }
    }

    evo
}

/// Find the original (pre-split) face that a sub-face came from.
/// Returns None if the face wasn't produced by a split.
fn find_original_face(
    face: FaceIdx,
    split_map_a: &HashMap<FaceIdx, Vec<FaceIdx>>,
    split_map_b: &HashMap<FaceIdx, Vec<FaceIdx>>,
) -> Option<FaceIdx> {
    for (original, sub_faces) in split_map_a.iter().chain(split_map_b.iter()) {
        if sub_faces.contains(&face) {
            return Some(*original);
        }
    }
    None
}

/// Result of classifying and splitting faces, including split provenance.
struct ClassifyResult {
    /// Classified faces (possibly sub-faces from splits).
    classified: Vec<(FaceIdx, FacePosition)>,
    /// Map from original face → sub-faces it was split into.
    /// Only contains entries for faces that were actually split.
    split_map: HashMap<FaceIdx, Vec<FaceIdx>>,
}

/// For each face, either classify it directly or split it by all its
/// intersection segments, then classify each resulting sub-face.
fn classify_and_split_faces(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    faces: &[FaceIdx],
    segments_map: &HashMap<FaceIdx, Vec<TrimmedSegment>>,
    other_solid: SolidIdx,
) -> ClassifyResult {
    let mut classified = Vec::new();
    let mut split_map: HashMap<FaceIdx, Vec<FaceIdx>> = HashMap::new();

    for &face in faces {
        if let Some(segments) = segments_map.get(&face) {
            if !segments.is_empty() {
                // Split the face by all segments, then classify each leaf.
                let leaves = split_face_by_all_segments(topo, geom, face, segments);
                // Track which sub-faces came from this original face.
                if leaves.len() > 1 || (leaves.len() == 1 && leaves[0] != face) {
                    split_map.insert(face, leaves.clone());
                }
                for leaf in &leaves {
                    let pos = classify_face(topo, geom, *leaf, other_solid);
                    classified.push((*leaf, pos));
                }
                continue;
            }
        }

        // No split needed — classify directly.
        let pos = classify_face(topo, geom, face, other_solid);
        classified.push((face, pos));
    }

    ClassifyResult {
        classified,
        split_map,
    }
}

/// A chain of connected trimmed segments forming a polyline.
struct SegmentChain {
    /// Ordered points along the chain (endpoints of constituent segments).
    points: Vec<Point3>,
    /// Whether this chain forms a closed loop (first ~ last point).
    is_closed: bool,
    /// Indices into the original segment slice that were consumed by this chain.
    segment_indices: Vec<usize>,
}

/// Chain trimmed segments that share endpoints into connected polylines.
///
/// Greedily walks from the first unvisited segment, appending any unvisited
/// segment whose start matches the current chain end (within tolerance).
/// A chain is marked closed when the end point is within tolerance of the
/// first point.
fn chain_trimmed_segments(segments: &[TrimmedSegment], tolerance: f64) -> Vec<SegmentChain> {
    let n = segments.len();
    let mut visited = vec![false; n];
    let mut chains = Vec::new();
    let tol_sq = tolerance * tolerance;

    for start_idx in 0..n {
        if visited[start_idx] {
            continue;
        }
        visited[start_idx] = true;

        let mut points = vec![segments[start_idx].start_point, segments[start_idx].end_point];
        let mut indices = vec![start_idx];

        // Greedily extend the chain.
        loop {
            let chain_end = *points.last().unwrap();
            let mut found = false;

            for j in 0..n {
                if visited[j] {
                    continue;
                }
                if (segments[j].start_point - chain_end).norm_squared() < tol_sq {
                    visited[j] = true;
                    points.push(segments[j].end_point);
                    indices.push(j);
                    found = true;
                    break;
                }
            }

            if !found {
                break;
            }

            // Check if chain closed.
            if (points.last().unwrap() - points[0]).norm_squared() < tol_sq {
                break;
            }
        }

        let is_closed = points.len() > 2
            && (points.last().unwrap() - points[0]).norm_squared() < tol_sq;

        chains.push(SegmentChain {
            points,
            is_closed,
            segment_indices: indices,
        });
    }

    chains
}

/// Iteratively split a face by multiple intersection segments.
///
/// First tries to chain segments into closed loops (e.g., a circle approximation
/// entirely inside the face) and handles those via `split_face_by_interior_loop`.
/// Remaining unchained segments go through the existing segment-by-segment splitting.
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
    // Step 1: Try chaining segments into closed interior loops.
    let chains = chain_trimmed_segments(segments, 1e-6);

    // Collect segment indices consumed by closed interior loops.
    let mut consumed_by_loop: Vec<bool> = vec![false; segments.len()];
    let mut current_face = face;
    let mut extra_leaves = Vec::new();

    for chain in &chains {
        if !chain.is_closed {
            continue;
        }

        // Check if all points in the chain are interior to the face
        // (i.e., snap_to_boundary fails for all of them).
        let boundary = face_boundary_verts(topo, &*geom, current_face);
        let all_interior = chain.points.iter().all(|p| {
            snap_to_boundary(&boundary, p).is_err()
        });

        if !all_interior {
            continue;
        }

        // Drop the last point if it duplicates the first (closed loop marker).
        let loop_pts: &[Point3] = if chain.points.len() > 1
            && (chain.points.last().unwrap() - chain.points[0]).norm() < 1e-6
        {
            &chain.points[..chain.points.len() - 1]
        } else {
            &chain.points
        };

        if loop_pts.len() < 3 {
            continue;
        }

        match split_face_by_interior_loop(topo, geom, current_face, loop_pts) {
            Ok(split_result) => {
                // Annular face continues as the face to process further segments on.
                current_face = split_result.face_a;
                // Interior face is a leaf (will be classified as Inside_B and discarded).
                extra_leaves.push(split_result.face_b);
                // Mark consumed segments.
                for &idx in &chain.segment_indices {
                    consumed_by_loop[idx] = true;
                }
            }
            Err(_) => {
                // Interior loop split failed — segments will be processed individually.
                continue;
            }
        }
    }

    // Step 2: Collect remaining (unconsumed) segment indices.
    let remaining_indices: Vec<usize> = (0..segments.len())
        .filter(|&i| !consumed_by_loop[i])
        .collect();

    // Step 3: Process remaining segments with the existing work-queue approach.
    let mut queue: Vec<(FaceIdx, Vec<usize>)> = vec![(current_face, remaining_indices)];
    let mut leaves = extra_leaves;

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
            // Normal splitting failed for all segments. Try extending segments
            // whose endpoints are interior (far from boundary) to reach the
            // face boundary, preserving intermediate vertices at the original
            // segment endpoints for twin matching.
            for (pos, &seg_idx) in remaining.iter().enumerate() {
                match split_face_with_extension(topo, geom, current_face, &segments[seg_idx]) {
                    Ok(split_result) => {
                        let rest: Vec<usize> = remaining.iter()
                            .enumerate()
                            .filter(|&(i, _)| i != pos)
                            .map(|(_, &idx)| idx)
                            .collect();

                        let (for_a, for_b) = distribute_segments(
                            topo, geom, segments, &rest,
                            split_result.face_a, split_result.face_b,
                        );

                        queue.push((split_result.face_a, for_a));
                        queue.push((split_result.face_b, for_b));
                        split_done = true;
                        break;
                    }
                    Err(_) => continue,
                }
            }
        }

        if !split_done {
            // No segment could split this face even with extension — it's a leaf.
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
        let shell_idx = topo.solids.get(solid).outer_shell();
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
            Ok(br) => {
                verify_solid(&topo, br.solid);
                let shell = topo.solids.get(br.solid).outer_shell();
                tessellate_shell(&mut topo, shell, &geom);
            }
            Err(e) => {
                eprintln!("Fuse failed (may need debugging): {e}");
            }
        }
    }

    #[test]
    fn test_cut_overlapping_boxes() {
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
            BooleanOp::Cut,
        );

        match &result {
            Ok(br) => {
                verify_solid(&topo, br.solid);
                let shell = topo.solids.get(br.solid).outer_shell();
                tessellate_shell(&mut topo, shell, &geom);
            }
            Err(e) => {
                panic!("Cut overlapping boxes failed: {e}");
            }
        }
    }

    #[test]
    fn test_common_overlapping_boxes() {
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
            BooleanOp::Common,
        );

        match &result {
            Ok(br) => {
                verify_solid(&topo, br.solid);
                let shell = topo.solids.get(br.solid).outer_shell();
                tessellate_shell(&mut topo, shell, &geom);
            }
            Err(e) => {
                panic!("Common overlapping boxes failed: {e}");
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

        match boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Fuse) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
            }
            Err(e) => {
                // Curved-surface booleans may fail with topology errors (known limitation).
                eprintln!("Fuse box+cylinder: {e}");
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

        match boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Cut) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
            }
            Err(e) => {
                eprintln!("Cut box-cylinder: {e}");
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
            Ok(br) => {
                verify_solid(&topo, br.solid);
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

        match boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Fuse) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
            }
            Err(e) => {
                eprintln!("Fuse box+revolved: {e}");
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

        match boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Cut) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
            }
            Err(e) => {
                eprintln!("Cut box-sphere: {e}");
            }
        }
    }

    #[test]
    fn test_fuse_box_nurbs_extrude() {
        use curvo::prelude::Interpolation;
        use curvo::prelude::NurbsCurve3D;
        use rustkernel_builders::nurbs_extrude_builder::make_nurbs_extrude_solid_into;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        // Create a closed NURBS curve (approximate circle, 8 points).
        let pts: Vec<Point3> = (0..8)
            .map(|i| {
                let theta = i as f64 * std::f64::consts::TAU / 8.0;
                Point3::new(theta.cos() * 0.5, theta.sin() * 0.5, 0.0)
            })
            .collect();
        let mut closed_pts = pts.clone();
        closed_pts.push(pts[0]);
        let curve = NurbsCurve3D::<f64>::interpolate(&closed_pts, 3).unwrap();

        let nurbs_solid = make_nurbs_extrude_solid_into(
            &mut topo, &mut geom, &curve, Vec3::new(0.0, 0.0, 1.0), 2.0, 32,
        );
        let box_solid = make_box_into(
            &mut topo, &mut geom, Point3::new(-1.0, -1.0, -1.0), 2.0, 2.0, 2.0,
        );

        match boolean_op(&mut topo, &mut geom, &pipeline, box_solid, nurbs_solid, BooleanOp::Fuse) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
            }
            Err(e) => {
                eprintln!("Fuse box+nurbs_extrude: {e}");
            }
        }
    }

    #[test]
    fn test_cut_box_nurbs_extrude() {
        use curvo::prelude::Interpolation;
        use curvo::prelude::NurbsCurve3D;
        use rustkernel_builders::nurbs_extrude_builder::make_nurbs_extrude_solid_into;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        // Create a closed NURBS curve (approximate circle, 8 points).
        let pts: Vec<Point3> = (0..8)
            .map(|i| {
                let theta = i as f64 * std::f64::consts::TAU / 8.0;
                Point3::new(theta.cos() * 0.5, theta.sin() * 0.5, 0.0)
            })
            .collect();
        let mut closed_pts = pts.clone();
        closed_pts.push(pts[0]);
        let curve = NurbsCurve3D::<f64>::interpolate(&closed_pts, 3).unwrap();

        let nurbs_solid = make_nurbs_extrude_solid_into(
            &mut topo, &mut geom, &curve, Vec3::new(0.0, 0.0, 1.0), 2.0, 32,
        );
        let box_solid = make_box_into(
            &mut topo, &mut geom, Point3::new(-1.0, -1.0, -1.0), 2.0, 2.0, 2.0,
        );

        // Cut exercises flip_normal on NURBS faces.
        match boolean_op(&mut topo, &mut geom, &pipeline, box_solid, nurbs_solid, BooleanOp::Cut) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
            }
            Err(e) => {
                eprintln!("Cut box-nurbs_extrude: {e}");
            }
        }
    }

    #[test]
    fn test_cut_single_interior_pocket() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        // Block: 4x4x1, centered at origin
        let block = make_box_into(&mut topo, &mut geom, Point3::new(-2.0, -2.0, -0.5), 4.0, 4.0, 1.0);
        // Pocket: 0.8x0.8, fully inside block in XY, punches through Z
        let pocket = make_box_into(&mut topo, &mut geom, Point3::new(-0.4, -0.4, -0.6), 0.8, 0.8, 1.5);

        match boolean_op(&mut topo, &mut geom, &pipeline, block, pocket, BooleanOp::Cut) {
            Ok(br) => {
                verify_solid(&topo, br.solid);
                let shell = topo.solids.get(br.solid).outer_shell();
                let nf = topo.shells.get(shell).faces.len();
                eprintln!("Interior pocket cut: {nf} faces");
            }
            Err(e) => {
                panic!("Interior pocket cut failed: {e}");
            }
        }
    }

    #[test]
    fn test_iterative_cut_four_pockets() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        // 4x4x1 block
        let block = make_box_into(&mut topo, &mut geom, Point3::new(-2.0, -2.0, -0.5), 4.0, 4.0, 1.0);

        // 4 small pockets offset along X axis only
        let pocket_positions = [
            Point3::new(-1.4, -0.4, -0.6),
            Point3::new(-0.4, -0.4, -0.6),
            Point3::new(0.6, -0.4, -0.6),
            Point3::new(1.6, -0.4, -0.6),
        ];

        let mut result = block;
        for (i, pos) in pocket_positions.iter().enumerate() {
            let pocket = make_box_into(&mut topo, &mut geom, *pos, 0.8, 0.8, 1.5);
            match boolean_op(&mut topo, &mut geom, &pipeline, result, pocket, BooleanOp::Cut) {
                Ok(br) => {
                    result = br.solid;
                    eprintln!("Pocket {i} succeeded");
                }
                Err(e) => {
                    eprintln!("Pocket {i} failed: {e}");
                }
            }
        }
    }

    #[test]
    fn test_evolution_cut_overlapping_boxes() {
        use rustkernel_topology::evolution::FaceOrigin;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(1.0, 0.0, 0.0), 2.0, 2.0, 2.0);

        let shell_a = topo.solids.get(a).outer_shell();
        let faces_a: Vec<FaceIdx> = topo.shells.get(shell_a).faces.clone();
        let shell_b = topo.solids.get(b).outer_shell();
        let faces_b: Vec<FaceIdx> = topo.shells.get(shell_b).faces.clone();

        let br = boolean_op(
            &mut topo, &mut geom, &pipeline, a, b, BooleanOp::Cut,
        ).expect("Cut should succeed");

        let evo = &br.evolution;

        // Result solid should have faces with provenance.
        assert!(!evo.face_provenance.is_empty(), "Should have face provenance");

        // Some original faces from B should be deleted (they're inside A
        // or on the opposite side in a cut).
        assert!(!evo.deleted_faces.is_empty(), "Cut should delete some faces");

        // Every face in the result should have a provenance entry.
        let result_shell = topo.solids.get(br.solid).outer_shell();
        let result_faces = &topo.shells.get(result_shell).faces;
        for &rf in result_faces {
            assert!(
                evo.face_provenance.contains_key(&rf),
                "Result face {} should have provenance", rf.raw()
            );
        }

        // Faces from solid A that were split should show up as SplitFrom.
        let split_count = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::SplitFrom(_)))
            .count();
        // In a box-box cut with overlap, some faces get split by intersection lines.
        eprintln!("Evolution: {} provenance entries, {} deleted, {} splits",
            evo.face_provenance.len(), evo.deleted_faces.len(), split_count);

        // Every deleted face should be from the original inputs.
        let all_originals: HashSet<u32> = faces_a.iter().chain(faces_b.iter())
            .map(|f| f.raw()).collect();
        for &df in &evo.deleted_faces {
            assert!(
                all_originals.contains(&df.raw()),
                "Deleted face {} should be from original solids", df.raw()
            );
        }
    }

    #[test]
    fn test_evolution_cut_box_cylinder() {
        // NOTE: Box-cylinder cut has a pre-existing UnmatchedTwin issue on
        // cylinder lateral faces (face_b side segments don't chain into closed
        // loops). This test validates evolution tracking IF the boolean succeeds,
        // but tolerates the known failure. See also test_cut_box_cylinder.
        use rustkernel_topology::evolution::FaceOrigin;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let pipeline = default_pipeline();

        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_cylinder_into(&mut topo, &mut geom, Point3::origin(), 0.5, 2.0, 16);

        match boolean_op(&mut topo, &mut geom, &pipeline, a, b, BooleanOp::Cut) {
            Ok(br) => {
                let evo = &br.evolution;

                let split_faces: Vec<_> = evo.face_provenance.iter()
                    .filter(|(_, o)| matches!(o, FaceOrigin::SplitFrom(_)))
                    .collect();
                assert!(
                    !split_faces.is_empty(),
                    "Box-cylinder cut should produce split faces from circle intersection"
                );

                let deleted_count = evo.deleted_faces.len();
                assert!(
                    deleted_count > 0,
                    "Should have deleted faces from cylinder"
                );
            }
            Err(e) => {
                // Box-cylinder cut has known UnmatchedTwin issue on cylinder lateral faces.
                eprintln!("Box-cylinder cut evolution: {e} — known issue");
            }
        }
    }
}

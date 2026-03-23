use rustkernel_math::Point3;
use rustkernel_topology::arena::Idx;
use rustkernel_topology::geom_store::GeomAccess;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use crate::curve_trimming::TrimmedSegment;
use rustkernel_geom::AnalyticalGeomStore;
use rustkernel_geom::LineSegment;

/// Errors from face splitting.
#[derive(Debug)]
pub enum SplitError {
    /// A segment endpoint is too far from any face boundary edge.
    EndpointNotOnBoundary { dist: f64 },
    /// The split is degenerate (enter and exit on same edge at same position).
    DegenerateSplit,
}

impl std::fmt::Display for SplitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SplitError::EndpointNotOnBoundary { dist } => {
                write!(f, "segment endpoint not on face boundary (dist={dist})")
            }
            SplitError::DegenerateSplit => write!(f, "degenerate split (enter == exit)"),
        }
    }
}

/// Result of splitting a face along an intersection segment.
pub struct SplitResult {
    /// The two new faces created by the split.
    pub face_a: FaceIdx,
    pub face_b: FaceIdx,
    /// The vertices created at the split points.
    pub vertex_start: VertexIdx,
    pub vertex_end: VertexIdx,
}

/// Split a face along a line segment that crosses it.
/// The original face is not modified (new topology is created).
///
/// Returns the two new faces, or `SplitError` if endpoints can't be snapped
/// to the boundary or the split is degenerate.
pub fn split_face_along_segment(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    face: FaceIdx,
    segment: &TrimmedSegment,
) -> Result<SplitResult, SplitError> {
    let surface_id = topo.faces.get(face).surface_id;
    let shell_idx = topo.faces.get(face).shell;

    // Get face boundary vertices and half-edges.
    let loop_idx = topo.faces.get(face).outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    let mut he_list = Vec::new();
    let mut vert_positions = Vec::new();
    let mut he = start_he;
    loop {
        let vert_idx = topo.half_edges.get(he).origin;
        let pid = topo.vertices.get(vert_idx).point_id;
        let pos = geom.point(pid);
        vert_positions.push(pos);
        he_list.push(he);
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }

    let n = vert_positions.len();

    // Snap segment endpoints to the face boundary (relaxed tolerance).
    let (enter_edge, snapped_start) = snap_to_boundary(&vert_positions, &segment.start_point)?;
    let (exit_edge, snapped_end) = snap_to_boundary(&vert_positions, &segment.end_point)?;

    // Guard against degenerate splits.
    // Same edge: segment runs along a boundary edge — always degenerate for convex faces.
    if enter_edge == exit_edge {
        return Err(SplitError::DegenerateSplit);
    }
    // Adjacent edges with snap at the shared vertex: segment runs along an edge.
    if (snapped_start - snapped_end).norm() < 1e-10 {
        return Err(SplitError::DegenerateSplit);
    }
    // Check if both snap points coincide with the same vertex (vertex-to-itself).
    let next_enter = (enter_edge + 1) % n;
    if next_enter == exit_edge && (snapped_start - vert_positions[next_enter]).norm() < 1e-10 {
        return Err(SplitError::DegenerateSplit);
    }

    // Create vertices for the snapped split points.
    let start_pid = geom.add_point(snapped_start);
    let end_pid = geom.add_point(snapped_end);
    let start_vert = topo.vertices.alloc(Vertex { point_id: start_pid });
    let end_vert = topo.vertices.alloc(Vertex { point_id: end_pid });

    // Build two sub-polygon vertex sequences.
    // Polygon A: from enter_point, walk forward to exit_point.
    // Polygon B: from exit_point, walk forward to enter_point.
    let poly_a_verts = build_sub_polygon(
        topo, geom, &vert_positions, enter_edge, exit_edge,
        start_vert, end_vert, n,
    );
    let poly_b_verts = build_sub_polygon(
        topo, geom, &vert_positions, exit_edge, enter_edge,
        end_vert, start_vert, n,
    );

    // Reject if either sub-polygon is degenerate (fewer than 3 vertices).
    if poly_a_verts.len() < 3 || poly_b_verts.len() < 3 {
        return Err(SplitError::DegenerateSplit);
    }

    // Build face topology for each sub-polygon.
    let face_a = build_face_from_verts(topo, geom, &poly_a_verts, surface_id, shell_idx);
    let face_b = build_face_from_verts(topo, geom, &poly_b_verts, surface_id, shell_idx);

    Ok(SplitResult {
        face_a,
        face_b,
        vertex_start: start_vert,
        vertex_end: end_vert,
    })
}

/// Snap a point to the nearest face boundary edge.
/// Returns `(edge_index, snapped_point)` or `SplitError` if too far.
///
/// Tolerance: 1e-6 (relaxed from the original 1e-8 to handle floating-point
/// drift from 2D→3D coordinate reconstruction in curve trimming).
///
/// When the snap point lands at a face vertex (end of an edge), the edge index
/// is advanced to the NEXT edge so that `build_sub_polygon` starts walking from
/// the vertex AFTER the snap point, avoiding duplicate vertices.
pub(crate) fn snap_to_boundary(verts: &[Point3], point: &Point3) -> Result<(usize, Point3), SplitError> {
    let tol = 1e-6;
    let n = verts.len();
    let mut best_edge = 0;
    let mut best_dist = f64::MAX;
    let mut best_projected = *point;

    for i in 0..n {
        let j = (i + 1) % n;
        let (dist, projected) = project_point_onto_segment(point, &verts[i], &verts[j]);
        if dist < best_dist {
            best_dist = dist;
            best_edge = i;
            best_projected = projected;
        }
    }

    if best_dist > tol {
        return Err(SplitError::EndpointNotOnBoundary { dist: best_dist });
    }

    // If the snapped point is at the END of the edge (close to the next vertex),
    // advance to the next edge. This ensures build_sub_polygon doesn't create
    // a zero-length edge from the snap vertex to the coincident boundary vertex.
    let next_vert = (best_edge + 1) % n;
    if (best_projected - verts[next_vert]).norm() < tol {
        best_edge = next_vert;
        best_projected = verts[next_vert];
    }
    // If snapped point is at the START of the edge, snap to exact vertex position.
    if (best_projected - verts[best_edge]).norm() < tol {
        best_projected = verts[best_edge];
    }

    Ok((best_edge, best_projected))
}

/// Project a point onto a line segment, returning (distance, projected_point).
fn project_point_onto_segment(p: &Point3, a: &Point3, b: &Point3) -> (f64, Point3) {
    let ab = b - a;
    let ap = p - a;
    let len_sq = ab.dot(&ab);
    if len_sq < 1e-20 {
        return (ap.norm(), *a);
    }
    let t = (ap.dot(&ab) / len_sq).clamp(0.0, 1.0);
    let proj = a + t * ab;
    let dist = (p - proj).norm();
    (dist, proj)
}

/// Build a sub-polygon vertex list by walking from edge `from_edge` to edge `to_edge`.
/// The first vertex is `start_vert` (the split point on `from_edge`),
/// then the original vertices from `from_edge + 1` to `to_edge`,
/// then `end_vert` (the split point on `to_edge`).
///
/// Boundary vertices that coincide with the start or end split point are skipped
/// to avoid creating zero-length edges (self-loops).
fn build_sub_polygon(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    vert_positions: &[Point3],
    from_edge: usize,
    to_edge: usize,
    start_vert: VertexIdx,
    end_vert: VertexIdx,
    n: usize,
) -> Vec<VertexIdx> {
    let tol = 1e-8;
    let start_pos = geom.point(topo.vertices.get(start_vert).point_id);
    let end_pos = geom.point(topo.vertices.get(end_vert).point_id);

    let mut verts = Vec::new();
    verts.push(start_vert);

    // Walk from (from_edge + 1) to (to_edge) inclusive, wrapping around.
    let mut i = (from_edge + 1) % n;
    loop {
        let pos = vert_positions[i];
        // Skip boundary vertices that coincide with start or end snap points.
        let skip = (pos - start_pos).norm() < tol || (pos - end_pos).norm() < tol;

        if !skip {
            let pid = geom.add_point(pos);
            let vert = topo.vertices.alloc(Vertex { point_id: pid });
            verts.push(vert);
        }

        if i == to_edge {
            break;
        }
        i = (i + 1) % n;
    }

    verts.push(end_vert);
    verts
}

/// Build a complete face (Face, Loop, HalfEdges, Edges) from a vertex sequence.
pub(crate) fn build_face_from_verts(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    verts: &[VertexIdx],
    surface_id: u32,
    shell_idx: ShellIdx,
) -> FaceIdx {
    let n = verts.len();
    assert!(n >= 3);

    let face_idx = topo.faces.alloc(Face {
        outer_loop: Idx::from_raw(0), // placeholder
        surface_id,
        mesh_cache: None,
        shell: shell_idx,
    });
    let loop_idx = topo.loops.alloc(Loop {
        half_edge: Idx::from_raw(0), // placeholder
        face: face_idx,
    });
    topo.faces.get_mut(face_idx).outer_loop = loop_idx;

    let mut he_idxs = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let origin = verts[i];
        let dest = verts[j];

        let start_pt = geom.point(topo.vertices.get(origin).point_id);
        let end_pt = geom.point(topo.vertices.get(dest).point_id);
        let curve_id = geom.add_line_segment(LineSegment {
            start: start_pt,
            end: end_pt,
        });

        let edge_idx = topo.edges.alloc(Edge {
            half_edges: [Idx::from_raw(0), Idx::from_raw(0)],
            curve_id,
        });

        let he_idx = topo.half_edges.alloc(HalfEdge {
            origin,
            twin: None,
            next: Idx::from_raw(0),
            edge: edge_idx,
            loop_ref: loop_idx,
        });

        topo.edges.get_mut(edge_idx).half_edges[0] = he_idx;
        he_idxs.push(he_idx);
    }

    // Wire up next pointers.
    for i in 0..n {
        let next = he_idxs[(i + 1) % n];
        topo.half_edges.get_mut(he_idxs[i]).next = next;
    }
    topo.loops.get_mut(loop_idx).half_edge = he_idxs[0];

    face_idx
}

/// Find the nearest face boundary vertex to a point.
fn nearest_boundary_vertex(verts: &[Point3], point: &Point3) -> Point3 {
    *verts
        .iter()
        .min_by(|a, b| {
            let da = (*a - point).norm_squared();
            let db = (*b - point).norm_squared();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("nearest_boundary_vertex called with empty slice")
}

/// Split a face along a segment with interior endpoint(s) by bridging them
/// to the nearest face vertex.
///
/// When a segment endpoint is too far from any boundary edge to snap, it is
/// "bridged" to the nearest face vertex instead. This avoids creating new
/// vertices on shared boundary edges (which would break twin matching with
/// adjacent faces that don't have matching vertices).
///
/// The split path goes: entry → orig_start → orig_end → exit
/// where entry/exit are either snapped boundary points (if the original
/// endpoint was on the boundary) or the nearest face vertex (if interior).
pub fn split_face_with_extension(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    face: FaceIdx,
    segment: &TrimmedSegment,
) -> Result<SplitResult, SplitError> {
    let surface_id = topo.faces.get(face).surface_id;
    let shell_idx = topo.faces.get(face).shell;

    // Get face boundary.
    let loop_idx = topo.faces.get(face).outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;
    let mut vert_positions = Vec::new();
    let mut he = start_he;
    loop {
        let vid = topo.half_edges.get(he).origin;
        let pid = topo.vertices.get(vid).point_id;
        vert_positions.push(geom.point(pid));
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }
    let n = vert_positions.len();

    // For each endpoint: if it snaps to the boundary, use the snap point.
    // Otherwise, bridge to the nearest face vertex.
    let start_on_boundary = snap_to_boundary(&vert_positions, &segment.start_point).is_ok();
    let end_on_boundary = snap_to_boundary(&vert_positions, &segment.end_point).is_ok();

    if start_on_boundary && end_on_boundary {
        // Both on boundary — use normal split.
        return split_face_along_segment(topo, geom, face, segment);
    }

    let entry_point = if start_on_boundary {
        segment.start_point
    } else {
        nearest_boundary_vertex(&vert_positions, &segment.start_point)
    };
    let exit_point = if end_on_boundary {
        segment.end_point
    } else {
        nearest_boundary_vertex(&vert_positions, &segment.end_point)
    };

    // Snap entry and exit to boundary edges.
    let (enter_edge, snapped_entry) = snap_to_boundary(&vert_positions, &entry_point)?;
    let (exit_edge, snapped_exit) = snap_to_boundary(&vert_positions, &exit_point)?;

    if enter_edge == exit_edge {
        return Err(SplitError::DegenerateSplit);
    }

    // Create vertices for the split path.
    let entry_vert = topo
        .vertices
        .alloc(Vertex { point_id: geom.add_point(snapped_entry) });
    let exit_vert = topo
        .vertices
        .alloc(Vertex { point_id: geom.add_point(snapped_exit) });

    // Check if intermediate vertices are distinct from entry/exit vertices.
    let skip_ostart = (segment.start_point - snapped_entry).norm() < 1e-8;
    let skip_oend = (segment.end_point - snapped_exit).norm() < 1e-8;

    let ostart_vert = if skip_ostart {
        entry_vert
    } else {
        topo.vertices
            .alloc(Vertex { point_id: geom.add_point(segment.start_point) })
    };
    let oend_vert = if skip_oend {
        exit_vert
    } else {
        topo.vertices
            .alloc(Vertex { point_id: geom.add_point(segment.end_point) })
    };

    let tol = 1e-8;
    let entry_pos = snapped_entry;
    let exit_pos = snapped_exit;

    // Poly A: entry → boundary walk to exit → exit → (oend) → (ostart) → [wraps to entry]
    let mut poly_a = Vec::new();
    poly_a.push(entry_vert);
    let mut i = (enter_edge + 1) % n;
    loop {
        let pos = vert_positions[i];
        let skip = (pos - entry_pos).norm() < tol
            || (pos - exit_pos).norm() < tol
            || (!skip_ostart && (pos - segment.start_point).norm() < tol)
            || (!skip_oend && (pos - segment.end_point).norm() < tol);
        if !skip {
            poly_a.push(topo.vertices.alloc(Vertex { point_id: geom.add_point(pos) }));
        }
        if i == exit_edge {
            break;
        }
        i = (i + 1) % n;
    }
    poly_a.push(exit_vert);
    if !skip_oend {
        poly_a.push(oend_vert);
    }
    if !skip_ostart {
        poly_a.push(ostart_vert);
    }

    // Poly B: exit → boundary walk to entry → entry → (ostart) → (oend) → [wraps to exit]
    let mut poly_b = Vec::new();
    poly_b.push(exit_vert);
    i = (exit_edge + 1) % n;
    loop {
        let pos = vert_positions[i];
        let skip = (pos - entry_pos).norm() < tol
            || (pos - exit_pos).norm() < tol
            || (!skip_ostart && (pos - segment.start_point).norm() < tol)
            || (!skip_oend && (pos - segment.end_point).norm() < tol);
        if !skip {
            poly_b.push(topo.vertices.alloc(Vertex { point_id: geom.add_point(pos) }));
        }
        if i == enter_edge {
            break;
        }
        i = (i + 1) % n;
    }
    poly_b.push(entry_vert);
    if !skip_ostart {
        poly_b.push(ostart_vert);
    }
    if !skip_oend {
        poly_b.push(oend_vert);
    }

    if poly_a.len() < 3 || poly_b.len() < 3 {
        return Err(SplitError::DegenerateSplit);
    }

    let face_a = build_face_from_verts(topo, geom, &poly_a, surface_id, shell_idx);
    let face_b = build_face_from_verts(topo, geom, &poly_b, surface_id, shell_idx);

    Ok(SplitResult {
        face_a,
        face_b,
        vertex_start: ostart_vert,
        vertex_end: oend_vert,
    })
}

/// Split a face by a closed polyline that lies entirely inside the face.
/// Creates a "slit" from the nearest boundary vertex to the nearest loop point,
/// then builds two sub-faces: an annular region and an interior region.
///
/// The annular face polygon is:
///   boundary_walk(slit_vert → ... → slit_vert_copy) + slit_forward(slit_vert_copy → first_loop)
///   + loop_walk(first_loop → ... → last_loop) + slit_backward(last_loop → slit_vert)
///
/// The slit_vert and slit_vert_copy are geometrically identical; after
/// `merge_coincident_vertices` in topology_builder, they merge so the
/// slit forward/backward edges become proper twins.
///
/// The interior face is just the loop vertices in reverse winding.
pub fn split_face_by_interior_loop(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    face: FaceIdx,
    loop_points: &[Point3],
) -> Result<SplitResult, SplitError> {
    if loop_points.len() < 3 {
        return Err(SplitError::DegenerateSplit);
    }

    let surface_id = topo.faces.get(face).surface_id;
    let shell_idx = topo.faces.get(face).shell;

    // Get face boundary vertices and positions.
    let loop_idx = topo.faces.get(face).outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;
    let mut boundary_positions = Vec::new();
    let mut he = start_he;
    loop {
        let vid = topo.half_edges.get(he).origin;
        let pid = topo.vertices.get(vid).point_id;
        boundary_positions.push(geom.point(pid));
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }

    let n_boundary = boundary_positions.len();
    if n_boundary < 3 {
        return Err(SplitError::DegenerateSplit);
    }

    // Find the boundary vertex closest to any loop point.
    let mut best_boundary_idx = 0usize;
    let mut best_loop_idx = 0usize;
    let mut best_dist_sq = f64::MAX;

    for (bi, bpos) in boundary_positions.iter().enumerate() {
        for (li, lpos) in loop_points.iter().enumerate() {
            let d = (bpos - lpos).norm_squared();
            if d < best_dist_sq {
                best_dist_sq = d;
                best_boundary_idx = bi;
                best_loop_idx = li;
            }
        }
    }

    let n_loop = loop_points.len();

    // Create loop vertices (one set for the annular face).
    let mut loop_verts: Vec<VertexIdx> = Vec::with_capacity(n_loop);
    for lp in loop_points {
        let pid = geom.add_point(*lp);
        loop_verts.push(topo.vertices.alloc(Vertex { point_id: pid }));
    }

    // Build annular face polygon.
    // Structure:
    //   [0]:               slit_boundary_vert (boundary[best_boundary_idx])
    //   [1..n_boundary-1]: rest of boundary walk
    //   [n_boundary]:      slit_boundary_vert_COPY (same position, for closing boundary)
    //   [n_boundary+1]:    first loop vert (loop[best_loop_idx])
    //   ...
    //   [n_boundary+n_loop]: last loop vert
    //
    // Edges:
    //   [n_boundary-1] → [n_boundary]: last boundary edge (closes the boundary back to slit pos)
    //   [n_boundary] → [n_boundary+1]: slit forward (slit_vert_copy → first_loop_vert)
    //   [last] → [0]: slit backward (last_loop_vert → slit_vert)
    //
    // After vertex merge, slit forward and backward become proper twins.
    let mut annular_verts: Vec<VertexIdx> = Vec::with_capacity(n_boundary + 1 + n_loop);

    // Boundary walk: allocate fresh vertices for all boundary points.
    for offset in 0..n_boundary {
        let idx = (best_boundary_idx + offset) % n_boundary;
        let pid = geom.add_point(boundary_positions[idx]);
        annular_verts.push(topo.vertices.alloc(Vertex { point_id: pid }));
    }

    // Duplicate slit boundary vertex (geometrically same as [0]).
    let slit_pos = boundary_positions[best_boundary_idx];
    let slit_copy_pid = geom.add_point(slit_pos);
    annular_verts.push(topo.vertices.alloc(Vertex { point_id: slit_copy_pid }));

    // Loop walk starting from best_loop_idx.
    for offset in 0..n_loop {
        let idx = (best_loop_idx + offset) % n_loop;
        annular_verts.push(loop_verts[idx]);
    }

    // Build interior face polygon: loop vertices in reverse winding.
    let mut interior_verts: Vec<VertexIdx> = Vec::with_capacity(n_loop);
    for offset in 0..n_loop {
        let idx = (best_loop_idx + n_loop - offset) % n_loop;
        let pid = geom.add_point(loop_points[idx]);
        interior_verts.push(topo.vertices.alloc(Vertex { point_id: pid }));
    }

    if annular_verts.len() < 3 || interior_verts.len() < 3 {
        return Err(SplitError::DegenerateSplit);
    }

    let face_a = build_face_from_verts(topo, geom, &annular_verts, surface_id, shell_idx);
    let face_b = build_face_from_verts(topo, geom, &interior_verts, surface_id, shell_idx);

    // Twin-match the slit edges within the annular face.
    // Slit forward: half-edge at index n_boundary (from slit_copy → first_loop_vert).
    // Slit backward: half-edge at index (n_boundary + n_loop) which is the closing edge
    //   (from last_loop_vert → annular_verts[0] = slit_vert).
    // These are A→B and B→A for the slit, so they should be twins.
    {
        let annular_loop = topo.faces.get(face_a).outer_loop;
        let first_he = topo.loops.get(annular_loop).half_edge;
        let total = annular_verts.len();

        // Collect half-edge indices by walking the loop.
        let mut he_indices = Vec::with_capacity(total);
        let mut curr = first_he;
        for _ in 0..total {
            he_indices.push(curr);
            curr = topo.half_edges.get(curr).next;
        }

        // Slit forward is at index n_boundary, slit backward is the last edge.
        let slit_fwd_idx = n_boundary;
        let slit_bwd_idx = total - 1;

        if slit_fwd_idx < he_indices.len() && slit_bwd_idx < he_indices.len() {
            let he_fwd = he_indices[slit_fwd_idx];
            let he_bwd = he_indices[slit_bwd_idx];

            topo.half_edges.get_mut(he_fwd).twin = Some(he_bwd);
            topo.half_edges.get_mut(he_bwd).twin = Some(he_fwd);

            // Share the edge entity.
            let shared_edge = topo.half_edges.get(he_fwd).edge;
            topo.half_edges.get_mut(he_bwd).edge = shared_edge;
            topo.edges.get_mut(shared_edge).half_edges[1] = he_bwd;
        }
    }

    Ok(SplitResult {
        face_a,
        face_b,
        vertex_start: loop_verts[best_loop_idx],
        vertex_end: annular_verts[0], // slit boundary vert
    })
}

/// Get the ordered boundary vertex positions for a face.
pub(crate) fn face_boundary_verts(
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    face: FaceIdx,
) -> Vec<Point3> {
    let loop_idx = topo.faces.get(face).outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;
    let mut pts = Vec::new();
    let mut he = start_he;
    loop {
        let vid = topo.half_edges.get(he).origin;
        let pid = topo.vertices.get(vid).point_id;
        pts.push(geom.point(pid));
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }
    pts
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_geom::AnalyticalGeomStore;
    use rustkernel_math::Vec3;
    use rustkernel_topology::face_util::face_boundary_points;
    use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
    use rustkernel_topology::store::TopoStore;

    #[test]
    fn test_split_square_face() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Find the top face (+Z at z=1).
        let shell = topo.solids.get(solid).outer_shell();
        let top_face = topo.shells.get(shell).faces.iter().find(|&&f| {
            let sid = topo.faces.get(f).surface_id;
            match geom.surface_kind(sid) {
                SurfaceKind::Plane { normal, .. } => normal.z > 0.5,
                _ => false,
            }
        }).copied().unwrap();

        // Split the top face with a line from (-1, 0, 1) to (1, 0, 1).
        let segment = TrimmedSegment {
            t_start: 0.0,
            t_end: 1.0,
            start_point: Point3::new(-1.0, 0.0, 1.0),
            end_point: Point3::new(1.0, 0.0, 1.0),
        };

        let result = split_face_along_segment(&mut topo, &mut geom, top_face, &segment).unwrap();

        // Each sub-face should have 4 vertices (rectangles).
        let pts_a = face_boundary_points(&topo, &geom, result.face_a);
        let pts_b = face_boundary_points(&topo, &geom, result.face_b);

        assert_eq!(pts_a.len(), 4, "Sub-face A should have 4 vertices");
        assert_eq!(pts_b.len(), 4, "Sub-face B should have 4 vertices");

        // Areas should sum to original (2x2 = 4).
        let area_a = polygon_area_3d(&pts_a);
        let area_b = polygon_area_3d(&pts_b);
        let total = area_a + area_b;
        assert!(
            (total - 4.0).abs() < 1e-6,
            "Areas should sum to 4.0, got {total} ({area_a} + {area_b})"
        );
    }

    #[test]
    fn test_split_snap_tolerance() {
        // Endpoint offset by ~1e-7 from edge → should succeed with relaxed tolerance.
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let shell = topo.solids.get(solid).outer_shell();
        let top_face = topo.shells.get(shell).faces.iter().find(|&&f| {
            let sid = topo.faces.get(f).surface_id;
            match geom.surface_kind(sid) {
                SurfaceKind::Plane { normal, .. } => normal.z > 0.5,
                _ => false,
            }
        }).copied().unwrap();

        // Offset start_point slightly off the edge (1e-7 in z).
        let segment = TrimmedSegment {
            t_start: 0.0,
            t_end: 1.0,
            start_point: Point3::new(-1.0, 0.0, 1.0 + 1e-7),
            end_point: Point3::new(1.0, 0.0, 1.0 + 1e-7),
        };

        let result = split_face_along_segment(&mut topo, &mut geom, top_face, &segment);
        assert!(result.is_ok(), "Split should succeed with 1e-7 drift");
    }

    #[test]
    fn test_split_degenerate_returns_error() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let shell = topo.solids.get(solid).outer_shell();
        let top_face = topo.shells.get(shell).faces.iter().find(|&&f| {
            let sid = topo.faces.get(f).surface_id;
            match geom.surface_kind(sid) {
                SurfaceKind::Plane { normal, .. } => normal.z > 0.5,
                _ => false,
            }
        }).copied().unwrap();

        // Both endpoints at same location on the same edge → degenerate.
        let segment = TrimmedSegment {
            t_start: 0.0,
            t_end: 0.0,
            start_point: Point3::new(-1.0, 0.0, 1.0),
            end_point: Point3::new(-1.0, 0.0, 1.0),
        };

        let result = split_face_along_segment(&mut topo, &mut geom, top_face, &segment);
        assert!(result.is_err(), "Degenerate split should return error");
    }

    #[test]
    fn test_split_far_endpoint_returns_error() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let shell = topo.solids.get(solid).outer_shell();
        let top_face = topo.shells.get(shell).faces.iter().find(|&&f| {
            let sid = topo.faces.get(f).surface_id;
            match geom.surface_kind(sid) {
                SurfaceKind::Plane { normal, .. } => normal.z > 0.5,
                _ => false,
            }
        }).copied().unwrap();

        // Endpoint far from any edge (0.1 away).
        let segment = TrimmedSegment {
            t_start: 0.0,
            t_end: 1.0,
            start_point: Point3::new(-1.0, 0.0, 1.1),
            end_point: Point3::new(1.0, 0.0, 1.1),
        };

        let result = split_face_along_segment(&mut topo, &mut geom, top_face, &segment);
        assert!(result.is_err(), "Far endpoint should return error");
    }

    /// Compute area of a planar polygon in 3D using cross products.
    fn polygon_area_3d(pts: &[Point3]) -> f64 {
        if pts.len() < 3 {
            return 0.0;
        }
        let mut area = Vec3::zeros();
        for i in 1..pts.len() - 1 {
            let ab = pts[i] - pts[0];
            let ac = pts[i + 1] - pts[0];
            area += ab.cross(&ac);
        }
        area.norm() * 0.5
    }
}

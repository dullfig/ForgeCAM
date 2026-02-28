use rustkernel_math::Point3;
use rustkernel_topology::arena::Idx;
use rustkernel_topology::geom_store::GeomAccess;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use crate::boolean::curve_trimming::TrimmedSegment;
use crate::geom::{AnalyticalGeomStore, LineSegment};

/// Result of splitting a face along an intersection segment.
pub struct SplitResult {
    /// The two new faces created by the split.
    pub face_a: FaceIdx,
    pub face_b: FaceIdx,
    /// The vertices created at the split points.
    pub vertex_start: VertexIdx,
    pub vertex_end: VertexIdx,
}

/// Split a planar face along a line segment that crosses it.
/// The original face is not modified (new topology is created).
///
/// Returns the two new faces. The segment's start and end points
/// must lie on the face boundary.
pub fn split_face_along_segment(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    face: FaceIdx,
    segment: &TrimmedSegment,
) -> SplitResult {
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

    // Find which edges the segment endpoints lie on.
    let enter_edge = find_edge_for_point(&vert_positions, &segment.start_point, 1e-8);
    let exit_edge = find_edge_for_point(&vert_positions, &segment.end_point, 1e-8);

    // Create vertices for the split points.
    let start_pid = geom.add_point(segment.start_point);
    let end_pid = geom.add_point(segment.end_point);
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

    // Build face topology for each sub-polygon.
    let face_a = build_face_from_verts(topo, geom, &poly_a_verts, surface_id, shell_idx);
    let face_b = build_face_from_verts(topo, geom, &poly_b_verts, surface_id, shell_idx);

    SplitResult {
        face_a,
        face_b,
        vertex_start: start_vert,
        vertex_end: end_vert,
    }
}

/// Find which edge of the polygon the point lies on.
/// Returns the edge index (edge i goes from vertex i to vertex (i+1) % n).
fn find_edge_for_point(verts: &[Point3], point: &Point3, tol: f64) -> usize {
    let n = verts.len();
    let mut best_edge = 0;
    let mut best_dist = f64::MAX;

    for i in 0..n {
        let j = (i + 1) % n;
        let dist = point_to_segment_distance(point, &verts[i], &verts[j]);
        if dist < best_dist {
            best_dist = dist;
            best_edge = i;
        }
    }

    assert!(
        best_dist < tol,
        "Segment endpoint not on face boundary (dist={best_dist})"
    );
    best_edge
}

/// Distance from a point to a line segment.
fn point_to_segment_distance(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = b - a;
    let ap = p - a;
    let len_sq = ab.dot(&ab);
    if len_sq < 1e-20 {
        return ap.norm();
    }
    let t = (ap.dot(&ab) / len_sq).clamp(0.0, 1.0);
    let proj = a + t * ab;
    (p - proj).norm()
}

/// Build a sub-polygon vertex list by walking from edge `from_edge` to edge `to_edge`.
/// The first vertex is `start_vert` (the split point on `from_edge`),
/// then the original vertices from `from_edge + 1` to `to_edge`,
/// then `end_vert` (the split point on `to_edge`).
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
    let mut verts = Vec::new();
    verts.push(start_vert);

    // Walk from (from_edge + 1) to (to_edge) inclusive, wrapping around.
    let mut i = (from_edge + 1) % n;
    loop {
        let pos = vert_positions[i];
        let pid = geom.add_point(pos);
        let vert = topo.vertices.alloc(Vertex { point_id: pid });
        verts.push(vert);

        if i == to_edge {
            break;
        }
        i = (i + 1) % n;
    }

    verts.push(end_vert);
    verts
}

/// Build a complete face (Face, Loop, HalfEdges, Edges) from a vertex sequence.
fn build_face_from_verts(
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
        let curve_id = geom.add_curve(LineSegment {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::Kernel;
    use rustkernel_math::Vec3;
    use rustkernel_topology::face_util::face_boundary_points;
    use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};

    #[test]
    fn test_split_square_face() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 2.0, 2.0);

        // Find the top face (+Z at z=1).
        let shell = k.topo.solids.get(solid).shell;
        let top_face = k.topo.shells.get(shell).faces.iter().find(|&&f| {
            let sid = k.topo.faces.get(f).surface_id;
            match k.geom.surface_kind(sid) {
                SurfaceKind::Plane { normal, .. } => normal.z > 0.5,
                _ => false,
            }
        }).copied().unwrap();

        // Split the top face with a line from (-1, 0, 1) to (1, 0, 1).
        // This cuts the square [(-1,-1), (1,-1), (1,1), (-1,1)] at z=1
        // horizontally through the middle.
        let segment = TrimmedSegment {
            t_start: 0.0,
            t_end: 1.0,
            start_point: Point3::new(-1.0, 0.0, 1.0),
            end_point: Point3::new(1.0, 0.0, 1.0),
        };

        let result = split_face_along_segment(&mut k.topo, &mut k.geom, top_face, &segment);

        // Each sub-face should have 4 vertices (rectangles).
        let pts_a = face_boundary_points(&k.topo, &k.geom, result.face_a);
        let pts_b = face_boundary_points(&k.topo, &k.geom, result.face_b);

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

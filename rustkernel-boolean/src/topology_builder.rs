use std::collections::HashMap;

use rustkernel_topology::arena::Idx;
use rustkernel_topology::geom_store::GeomAccess;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use crate::face_selector::SelectedFaces;
use rustkernel_geom::{AnalyticalGeomStore, LineSegment};

/// Errors during boolean result construction.
#[derive(Debug)]
pub enum BuildError {
    EulerViolation { v: usize, e: usize, f: usize },
    UnmatchedTwin(u32),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildError::EulerViolation { v, e, f: faces } => {
                write!(f, "Euler violation: V({v}) - E({e}) + F({faces}) != 2")
            }
            BuildError::UnmatchedTwin(he) => {
                write!(f, "Half-edge {he} has no twin")
            }
        }
    }
}

impl std::error::Error for BuildError {}

/// Build a result solid from selected faces.
///
/// Deep-copies all selected faces into fresh arena slots, optionally flips
/// normals, merges coincident vertices, matches twins, and validates.
pub fn build_result_solid(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    selected: &SelectedFaces,
) -> Result<SolidIdx, BuildError> {
    // Allocate result solid and shell.
    let result_solid = topo.solids.alloc(Solid {
        shell: Idx::from_raw(0),
        genus: 0,
    });
    let result_shell = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: result_solid,
    });
    topo.solids.get_mut(result_solid).shell = result_shell;

    // Deep-copy faces from A (normal orientation preserved).
    let mut all_new_faces = Vec::new();
    for &face_idx in &selected.keep_from_a {
        let new_face = deep_copy_face(topo, geom, face_idx, result_shell, false);
        all_new_faces.push(new_face);
    }

    // Deep-copy faces from B (normal orientation preserved).
    for &face_idx in &selected.keep_from_b {
        let new_face = deep_copy_face(topo, geom, face_idx, result_shell, false);
        all_new_faces.push(new_face);
    }

    // Deep-copy faces from B with flipped normals (for Cut operation).
    for &face_idx in &selected.flip_from_b {
        let new_face = deep_copy_face(topo, geom, face_idx, result_shell, true);
        all_new_faces.push(new_face);
    }

    // Add all new faces to the shell.
    for &nf in &all_new_faces {
        topo.shells.get_mut(result_shell).faces.push(nf);
    }

    // Merge coincident vertices.
    merge_coincident_vertices(topo, geom, &all_new_faces, 1e-8);

    // Match twins.
    match_twins(topo, geom, &all_new_faces);

    // Validate.
    validate(topo, &all_new_faces)?;

    Ok(result_solid)
}

/// Deep-copy a face's topology into fresh arena slots.
/// If `flip` is true, reverse the winding order and negate the surface normal.
fn deep_copy_face(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    face_idx: FaceIdx,
    shell: ShellIdx,
    flip: bool,
) -> FaceIdx {
    let old_face = topo.faces.get(face_idx);
    let old_surface_id = old_face.surface_id;
    let old_loop_idx = old_face.outer_loop;

    // Copy surface (possibly with flipped normal).
    let mut cloned_surface = geom.surfaces[old_surface_id as usize].clone();
    if flip {
        cloned_surface.flip_normal();
    }
    let new_surface_id = geom.add_surface(cloned_surface);

    // Collect old half-edge data.
    let start_he = topo.loops.get(old_loop_idx).half_edge;
    let mut old_he_data = Vec::new();
    let mut he = start_he;
    loop {
        let vert_idx = topo.half_edges.get(he).origin;
        let pid = topo.vertices.get(vert_idx).point_id;
        let pos = geom.point(pid);
        old_he_data.push(pos);
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }

    // If flipping, reverse the vertex order.
    if flip {
        old_he_data.reverse();
    }

    // Build new face topology.
    let n = old_he_data.len();

    let new_face = topo.faces.alloc(Face {
        outer_loop: Idx::from_raw(0),
        surface_id: new_surface_id,
        mesh_cache: None,
        shell,
    });
    let new_loop = topo.loops.alloc(Loop {
        half_edge: Idx::from_raw(0),
        face: new_face,
    });
    topo.faces.get_mut(new_face).outer_loop = new_loop;

    let mut new_verts = Vec::with_capacity(n);
    for pos in &old_he_data {
        let pid = geom.add_point(*pos);
        let vert = topo.vertices.alloc(Vertex { point_id: pid });
        new_verts.push(vert);
    }

    let mut he_idxs = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let start_pt = geom.point(topo.vertices.get(new_verts[i]).point_id);
        let end_pt = geom.point(topo.vertices.get(new_verts[j]).point_id);
        let curve_id = geom.add_line_segment(LineSegment {
            start: start_pt,
            end: end_pt,
        });

        let edge_idx = topo.edges.alloc(Edge {
            half_edges: [Idx::from_raw(0), Idx::from_raw(0)],
            curve_id,
        });

        let he_idx = topo.half_edges.alloc(HalfEdge {
            origin: new_verts[i],
            twin: None,
            next: Idx::from_raw(0),
            edge: edge_idx,
            loop_ref: new_loop,
        });

        topo.edges.get_mut(edge_idx).half_edges[0] = he_idx;
        he_idxs.push(he_idx);
    }

    for i in 0..n {
        let next = he_idxs[(i + 1) % n];
        topo.half_edges.get_mut(he_idxs[i]).next = next;
    }
    topo.loops.get_mut(new_loop).half_edge = he_idxs[0];

    new_face
}

/// Merge vertices that are geometrically coincident.
/// After merging, half-edges that pointed to merged-away vertices
/// are updated to point to the canonical vertex.
fn merge_coincident_vertices(
    topo: &mut TopoStore,
    geom: &AnalyticalGeomStore,
    faces: &[FaceIdx],
    tolerance: f64,
) {
    // Collect all vertices from these faces.
    let mut all_verts: Vec<VertexIdx> = Vec::new();
    for &face_idx in faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            all_verts.push(topo.half_edges.get(he).origin);
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }

    // Build merge map: vertex -> canonical vertex.
    let tol_sq = tolerance * tolerance;
    let mut canonical: HashMap<VertexIdx, VertexIdx> = HashMap::new();

    // Simple O(n^2) merge — fine for small vertex counts.
    let unique_verts: Vec<VertexIdx> = {
        let mut seen = Vec::new();
        for &v in &all_verts {
            if !seen.contains(&v) {
                seen.push(v);
            }
        }
        seen
    };

    for i in 0..unique_verts.len() {
        if canonical.contains_key(&unique_verts[i]) {
            continue;
        }
        let pi = geom.point(topo.vertices.get(unique_verts[i]).point_id);
        canonical.insert(unique_verts[i], unique_verts[i]);

        for j in (i + 1)..unique_verts.len() {
            if canonical.contains_key(&unique_verts[j]) {
                continue;
            }
            let pj = geom.point(topo.vertices.get(unique_verts[j]).point_id);
            if (pi - pj).norm_squared() < tol_sq {
                canonical.insert(unique_verts[j], unique_verts[i]);
            }
        }
    }

    // Update half-edge origins to canonical vertices.
    for &face_idx in faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let origin = topo.half_edges.get(he).origin;
            if let Some(&canon) = canonical.get(&origin) {
                if canon != origin {
                    topo.half_edges.get_mut(he).origin = canon;
                }
            }
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }
}

/// Match twins across all faces by looking for half-edge pairs (A→B) and (B→A)
/// that share the same canonical vertices.
fn match_twins(
    topo: &mut TopoStore,
    _geom: &AnalyticalGeomStore,
    faces: &[FaceIdx],
) {
    // Build map: (origin_vertex, dest_vertex) -> HalfEdgeIdx.
    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    for &face_idx in faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let origin = topo.half_edges.get(he).origin;
            let next_he = topo.half_edges.get(he).next;
            let _dest = topo.half_edges.get(next_he).origin;

            // Use vertex index (raw) as key since vertices are already merged.
            he_map.insert((origin.raw(), _dest.raw()), he);

            he = next_he;
            if he == start_he {
                break;
            }
        }
    }

    // Match twins.
    let keys: Vec<(u32, u32)> = he_map.keys().cloned().collect();
    for (a, b) in keys {
        if let (Some(&he_ab), Some(&he_ba)) = (he_map.get(&(a, b)), he_map.get(&(b, a))) {
            topo.half_edges.get_mut(he_ab).twin = Some(he_ba);
            topo.half_edges.get_mut(he_ba).twin = Some(he_ab);

            // Share edge entity.
            let shared_edge = topo.half_edges.get(he_ab).edge;
            topo.half_edges.get_mut(he_ba).edge = shared_edge;
            topo.edges.get_mut(shared_edge).half_edges[1] = he_ba;
        }
    }
}

/// Validate the result: Euler formula and all twins matched.
fn validate(topo: &TopoStore, faces: &[FaceIdx]) -> Result<(), BuildError> {
    let mut verts = std::collections::HashSet::new();
    let mut edges = std::collections::HashSet::new();
    let mut unmatched = Vec::new();

    for &face_idx in faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let he_data = topo.half_edges.get(he);
            verts.insert(he_data.origin.raw());
            edges.insert(he_data.edge.raw());
            if he_data.twin.is_none() {
                unmatched.push(he.raw());
            }
            he = he_data.next;
            if he == start_he {
                break;
            }
        }
    }

    if !unmatched.is_empty() {
        return Err(BuildError::UnmatchedTwin(unmatched[0]));
    }

    let v = verts.len();
    let e = edges.len();
    let f = faces.len();
    if v as i32 - e as i32 + f as i32 != 2 {
        return Err(BuildError::EulerViolation { v, e, f });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_geom::AnalyticalGeomStore;
    use rustkernel_math::Point3;
    use rustkernel_topology::store::TopoStore;

    #[test]
    fn test_sew_cube_from_standalone_faces() {
        // Create a box, then deep-copy all 6 faces as standalone faces
        // and reconstruct via build_result_solid.
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let shell_idx = topo.solids.get(solid).shell;
        let faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

        let selected = SelectedFaces {
            keep_from_a: faces,
            keep_from_b: Vec::new(),
            flip_from_b: Vec::new(),
        };

        let result = build_result_solid(&mut topo, &mut geom, &selected);
        assert!(result.is_ok(), "Should build valid solid: {:?}", result.err());
    }
}

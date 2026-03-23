use std::collections::HashMap;
use std::collections::HashSet;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::evolution::{FaceOrigin, ShapeEvolution};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{AnalyticalGeomStore, Plane, SurfaceDef};

use crate::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};
use crate::edge_analysis::{
    edge_adjacency, edge_convexity, edge_endpoints, face_centroid, plane_normal,
    EdgeConvexity,
};

/// Errors from the chamfer builder.
#[derive(Debug)]
pub enum ChamferError {
    EdgeAnalysis(crate::edge_analysis::EdgeAnalysisError),
    NotConvex(EdgeIdx),
    SharedVertex(VertexIdx),
}

impl std::fmt::Display for ChamferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChamferError::EdgeAnalysis(e) => write!(f, "Edge analysis error: {e}"),
            ChamferError::NotConvex(e) => write!(f, "Edge {:?} is not convex", e),
            ChamferError::SharedVertex(v) => {
                write!(f, "Vertex {:?} is shared by multiple chamfered edges", v)
            }
        }
    }
}

impl std::error::Error for ChamferError {}

impl From<crate::edge_analysis::EdgeAnalysisError> for ChamferError {
    fn from(e: crate::edge_analysis::EdgeAnalysisError) -> Self {
        ChamferError::EdgeAnalysis(e)
    }
}

/// Computed contact information for one chamfered edge.
struct ChamferContact {
    /// The edge being chamfered.
    edge: EdgeIdx,
    /// VertexIdx of endpoint A.
    vert_a: VertexIdx,
    /// VertexIdx of endpoint B.
    vert_b: VertexIdx,
    /// Face on the he_a side.
    face_a: FaceIdx,
    /// Face on the he_b side.
    face_b: FaceIdx,
    /// Contact point on face_a at endpoint A.
    c1_a: Point3,
    /// Contact point on face_a at endpoint B.
    c1_b: Point3,
    /// Contact point on face_b at endpoint A.
    c2_a: Point3,
    /// Contact point on face_b at endpoint B.
    c2_b: Point3,
}

/// Compute the inward direction from an edge toward a face's interior,
/// projected perpendicular to the edge direction and lying in the face's plane.
fn compute_inward(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    face_idx: FaceIdx,
    edge_midpoint: &Point3,
    edge_dir: &Vec3,
) -> Vec3 {
    let centroid = face_centroid(topo, geom, face_idx);
    let toward_face = centroid - edge_midpoint;
    // Remove component along edge direction.
    let perp = toward_face - edge_dir.scale(toward_face.dot(edge_dir));
    let len = perp.norm();
    if len < 1e-15 {
        // Degenerate: face centroid is on the edge line. Fall back to face normal cross edge.
        let n = plane_normal(geom, topo, face_idx);
        let fallback = n.cross(edge_dir);
        let fl = fallback.norm();
        if fl > 1e-15 { fallback / fl } else { Vec3::zeros() }
    } else {
        perp / len
    }
}

/// Build a chamfer on the specified edges of a solid.
///
/// Returns a new SolidIdx with chamfer faces inserted. The original solid's topology
/// is not modified (rebuild approach).
///
/// # Constraints
/// - All edges must be convex and have planar adjacent faces.
/// - No two edges may share a vertex.
pub fn chamfer_edges_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    edges: &[EdgeIdx],
    distance: f64,
) -> Result<(SolidIdx, ShapeEvolution), ChamferError> {
    let _span = info_span!("chamfer_edges", n_edges = edges.len(), distance).entered();

    // --- Phase 1: Validate ---
    let mut all_chamfer_verts: HashMap<VertexIdx, EdgeIdx> = HashMap::new();
    let mut contacts: Vec<ChamferContact> = Vec::with_capacity(edges.len());

    for &edge_idx in edges {
        let conv = edge_convexity(topo, geom, edge_idx)?;
        if conv != EdgeConvexity::Convex {
            return Err(ChamferError::NotConvex(edge_idx));
        }

        let adj = edge_adjacency(topo, edge_idx)?;
        let (pt_a, pt_b) = edge_endpoints(topo, geom, edge_idx);

        // Get the VertexIdx for A and B.
        let he_a_data = topo.half_edges.get(adj.he_a);
        let he_b_data = topo.half_edges.get(adj.he_b);
        let vert_a = he_a_data.origin;
        let vert_b = he_b_data.origin;

        // Check no shared vertex with another chamfered edge.
        for &v in &[vert_a, vert_b] {
            if let Some(&prev_edge) = all_chamfer_verts.get(&v) {
                if prev_edge != edge_idx {
                    return Err(ChamferError::SharedVertex(v));
                }
            }
            all_chamfer_verts.insert(v, edge_idx);
        }

        // Compute contact points.
        let edge_vec = pt_b - pt_a;
        let edge_len = edge_vec.norm();
        let edge_dir = if edge_len > 1e-15 { edge_vec / edge_len } else { Vec3::zeros() };
        let mid = Point3::from((pt_a.coords + pt_b.coords) * 0.5);

        let inward_1 = compute_inward(topo, geom, adj.face_a, &mid, &edge_dir);
        let inward_2 = compute_inward(topo, geom, adj.face_b, &mid, &edge_dir);

        contacts.push(ChamferContact {
            edge: edge_idx,
            vert_a,
            vert_b,
            face_a: adj.face_a,
            face_b: adj.face_b,
            c1_a: pt_a + distance * inward_1,
            c1_b: pt_b + distance * inward_1,
            c2_a: pt_a + distance * inward_2,
            c2_b: pt_b + distance * inward_2,
        });
    }

    // --- Phase 2: Build lookup maps ---
    // For each vertex on a chamfered edge, map (vertex, face) → replacement point(s).
    // If the face is face_a of the edge → replace vertex with c1 contact.
    // If the face is face_b → replace with c2 contact.
    // If the face is neither (an "end face") → replace with two contacts [c1, c2].

    // Map: (VertexIdx, FaceIdx) → Vec<Point3> replacement points.
    let mut vertex_face_replacements: HashMap<(u32, u32), Vec<Point3>> = HashMap::new();

    // Map: EdgeIdx → chamfer face polygon points + surface.
    let mut chamfer_faces: Vec<(Vec<Point3>, SurfaceDef)> = Vec::new();

    for contact in &contacts {
        // face_a: vertex → c1 contact
        vertex_face_replacements.insert(
            (contact.vert_a.raw(), contact.face_a.raw()),
            vec![contact.c1_a],
        );
        vertex_face_replacements.insert(
            (contact.vert_b.raw(), contact.face_a.raw()),
            vec![contact.c1_b],
        );

        // face_b: vertex → c2 contact
        vertex_face_replacements.insert(
            (contact.vert_a.raw(), contact.face_b.raw()),
            vec![contact.c2_a],
        );
        vertex_face_replacements.insert(
            (contact.vert_b.raw(), contact.face_b.raw()),
            vec![contact.c2_b],
        );

        // Chamfer face: [c1_A, c1_B, c2_B, c2_A] with outward winding.
        // The normal should point outward. We determine correct winding by checking
        // that the chamfer face normal is consistent with the adjacent faces.
        let chamfer_normal = (contact.c1_b - contact.c1_a)
            .cross(&(contact.c2_a - contact.c1_a));
        let cn_len = chamfer_normal.norm();
        let chamfer_normal = if cn_len > 1e-15 { chamfer_normal / cn_len } else { Vec3::zeros() };

        // Check: the chamfer normal should point away from face_a's interior.
        // If dot(chamfer_normal, inward_to_face_a_from_edge) < 0, it's correct.
        // We use face_a normal as a proxy — chamfer normal should have a positive component along face_a's normal.
        let n_a = plane_normal(geom, topo, contact.face_a);
        let n_b = plane_normal(geom, topo, contact.face_b);
        let expected_dir = (n_a + n_b).normalize();

        let poly = if chamfer_normal.dot(&expected_dir) > 0.0 {
            vec![contact.c1_a, contact.c1_b, contact.c2_b, contact.c2_a]
        } else {
            vec![contact.c2_a, contact.c2_b, contact.c1_b, contact.c1_a]
        };

        let surface = SurfaceDef::Plane(Plane {
            origin: contact.c1_a,
            normal: if chamfer_normal.dot(&expected_dir) > 0.0 {
                chamfer_normal
            } else {
                -chamfer_normal
            },
        });

        chamfer_faces.push((poly, surface));
    }

    // --- Phase 3: Handle "end faces" ---
    // An end face is a face that contains a chamfered vertex but is neither face_a nor face_b
    // of that vertex's chamfered edge. For end faces, the vertex gets replaced by two points.
    let shell_idx = topo.solids.get(solid).outer_shell();
    let face_list: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

    for &face_idx in &face_list {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let vert = topo.half_edges.get(he).origin;
            if all_chamfer_verts.contains_key(&vert) {
                let key = (vert.raw(), face_idx.raw());
                if !vertex_face_replacements.contains_key(&key) {
                    // This is an end face. Find which chamfered edge this vertex belongs to.
                    let &edge_idx = all_chamfer_verts.get(&vert).unwrap();
                    let contact = contacts.iter().find(|c| c.edge == edge_idx).unwrap();

                    // Determine the two replacement points (c1, c2) for this vertex.
                    let (c1, c2) = if vert == contact.vert_a {
                        (contact.c1_a, contact.c2_a)
                    } else {
                        (contact.c1_b, contact.c2_b)
                    };

                    // Determine ordering: the incoming half-edge's twin tells us which
                    // adjacent face we came from. If it's face_a, c1 comes first; if face_b, c2 first.
                    // We need to determine which order preserves correct winding for the end face.
                    //
                    // Strategy: the vertex is replaced by [c1, c2] or [c2, c1].
                    // The previous vertex in the loop → current vertex edge is on some face
                    // via the twin. If that neighboring face is face_a of the chamfer edge,
                    // then the first replacement should be c1 (near face_a), otherwise c2.
                    let prev_he = find_prev_he(topo, loop_idx, he);
                    let prev_twin_face = topo.half_edges.get(prev_he).twin
                        .map(|twin| {
                            let twin_loop = topo.half_edges.get(twin).loop_ref;
                            topo.loops.get(twin_loop).face
                        });

                    let ordered = if prev_twin_face == Some(contact.face_a) {
                        vec![c1, c2]
                    } else {
                        vec![c2, c1]
                    };

                    vertex_face_replacements.insert(key, ordered);
                }
            }
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }

    // --- Phase 4: Rebuild solid ---
    let mut evo = ShapeEvolution::new();

    // Collect faces adjacent to chamfered edges (they get Modified provenance).
    let mut adjacent_faces: HashSet<u32> = HashSet::new();
    for contact in &contacts {
        adjacent_faces.insert(contact.face_a.raw());
        adjacent_faces.insert(contact.face_b.raw());
    }

    // Record deleted edges (the chamfered edges).
    for contact in &contacts {
        evo.record_deleted_edge(contact.edge);
    }

    let new_solid_idx = topo.solids.alloc(Solid {
        shells: vec![Idx::from_raw(0)],
        genus: 0,
    });
    let new_shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: new_solid_idx,
    });
    topo.solids.get_mut(new_solid_idx).shells = vec![new_shell_idx];

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();
    let mut point_to_vert: HashMap<u64, VertexIdx> = HashMap::new();

    // Helper: register a point, returning its VertexIdx. Deduplicates by bit pattern.
    let mut get_or_create_vert = |geom: &mut AnalyticalGeomStore,
                                   topo: &mut TopoStore,
                                   pt: Point3,
                                   point_to_vert: &mut HashMap<u64, VertexIdx>|
     -> VertexIdx {
        let key = point_hash(&pt);
        if let Some(&vi) = point_to_vert.get(&key) {
            return vi;
        }
        let pid = geom.add_point(pt);
        let vi = topo.vertices.alloc(Vertex { point_id: pid });
        point_to_vert.insert(key, vi);
        vi
    };

    // Rebuild each original face with modified vertex polygons.
    for &face_idx in &face_list {
        let surface_id = topo.faces.get(face_idx).surface_id;
        // Clone the surface for the new face.
        let new_surface_id = geom.add_surface(geom.surfaces[surface_id as usize].clone());

        // Walk the loop and build the new polygon.
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut new_polygon: Vec<Point3> = Vec::new();
        let mut he = start_he;
        loop {
            let vert = topo.half_edges.get(he).origin;
            let key = (vert.raw(), face_idx.raw());
            if let Some(replacements) = vertex_face_replacements.get(&key) {
                new_polygon.extend_from_slice(replacements);
            } else {
                let pid = topo.vertices.get(vert).point_id;
                new_polygon.push(geom.points[pid as usize]);
            }
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }

        // Create vertex indices.
        let vert_idxs: Vec<VertexIdx> = new_polygon
            .iter()
            .map(|pt| get_or_create_vert(geom, topo, *pt, &mut point_to_vert))
            .collect();

        let new_face = build_face_from_vert_idxs(topo, geom, &vert_idxs, new_surface_id, new_shell_idx, &mut he_map);
        // Track face provenance.
        if adjacent_faces.contains(&face_idx.raw()) {
            evo.record_face(new_face, FaceOrigin::Modified(face_idx));
        } else {
            evo.record_face(new_face, FaceOrigin::CopiedFrom(face_idx));
        }
    }

    // Add chamfer faces.
    for (i, (poly, surface)) in chamfer_faces.into_iter().enumerate() {
        let new_surface_id = geom.add_surface(surface);
        let vert_idxs: Vec<VertexIdx> = poly
            .iter()
            .map(|pt| get_or_create_vert(geom, topo, *pt, &mut point_to_vert))
            .collect();
        let new_face = build_face_from_vert_idxs(topo, geom, &vert_idxs, new_surface_id, new_shell_idx, &mut he_map);
        evo.record_face(new_face, FaceOrigin::FromEdge(contacts[i].edge));
    }

    // Match twins.
    match_twins_from_map(topo, &he_map);

    Ok((new_solid_idx, evo))
}

/// Find the half-edge before `target_he` in the loop.
/// Delegates to the canonical implementation in `euler`.
fn find_prev_he(topo: &TopoStore, loop_idx: LoopIdx, target_he: HalfEdgeIdx) -> HalfEdgeIdx {
    rustkernel_topology::euler::find_prev_he(topo, loop_idx, target_he)
}

/// Hash a point by its bit representation for deduplication.
pub fn point_hash(pt: &Point3) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    pt.x.to_bits().hash(&mut hasher);
    pt.y.to_bits().hash(&mut hasher);
    pt.z.to_bits().hash(&mut hasher);
    hasher.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box_builder::make_box_into;
    use crate::edge_analysis::solid_edges;
    use rustkernel_topology::evolution::FaceOrigin;

    fn verify_euler(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).outer_shell();
        let faces = &topo.shells.get(shell_idx).faces;
        let genus = topo.solids.get(solid).genus;

        let mut verts = HashSet::new();
        let mut edges = HashSet::new();
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                verts.insert(topo.half_edges.get(he).origin.raw());
                edges.insert(topo.half_edges.get(he).edge.raw());
                he = topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }

        let v = verts.len() as i32;
        let e = edges.len() as i32;
        let f = faces.len() as i32;
        let expected = 2 - 2 * genus as i32;
        assert_eq!(
            v - e + f, expected,
            "Euler: V({v}) - E({e}) + F({f}) = {} != {expected}",
            v - e + f
        );
    }

    fn verify_all_twins(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).outer_shell();
        let faces = &topo.shells.get(shell_idx).faces;
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                assert!(
                    topo.half_edges.get(he).twin.is_some(),
                    "Half-edge {} has no twin",
                    he.raw()
                );
                he = topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }
    }

    fn face_count(topo: &TopoStore, solid: SolidIdx) -> usize {
        let shell_idx = topo.solids.get(solid).outer_shell();
        topo.shells.get(shell_idx).faces.len()
    }

    /// Pick one edge of a box between two specific face normals.
    fn find_box_edge_between(
        topo: &TopoStore,
        geom: &AnalyticalGeomStore,
        solid: SolidIdx,
        n1: Vec3,
        n2: Vec3,
    ) -> EdgeIdx {
        let all_edges = solid_edges(topo, solid);
        for &edge_idx in &all_edges {
            let adj = edge_adjacency(topo, edge_idx).unwrap();
            let fn_a = plane_normal(geom, topo, adj.face_a);
            let fn_b = plane_normal(geom, topo, adj.face_b);
            if (fn_a - n1).norm() < 0.01 && (fn_b - n2).norm() < 0.01 {
                return edge_idx;
            }
            if (fn_a - n2).norm() < 0.01 && (fn_b - n1).norm() < 0.01 {
                return edge_idx;
            }
        }
        panic!("No edge found between normals {:?} and {:?}", n1, n2);
    }

    #[test]
    fn test_chamfer_one_edge_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Chamfer one edge: between bottom (-Z) and front (-Y).
        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (chamfered, _evo) = chamfer_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();
        verify_euler(&topo, chamfered);
    }

    #[test]
    fn test_chamfer_one_edge_twins() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (chamfered, _evo) = chamfer_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();
        verify_all_twins(&topo, chamfered);
    }

    #[test]
    fn test_chamfer_one_edge_face_count() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (chamfered, _evo) = chamfer_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();
        assert_eq!(face_count(&topo, chamfered), 7, "Chamfer 1 edge → 7 faces");
    }

    #[test]
    fn test_chamfer_two_opposite_edges() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Pick two opposite edges (bottom-front and top-back) — these share no vertices.
        let e1 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );
        let e2 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
        );

        let (chamfered, _evo) = chamfer_edges_into(&mut topo, &mut geom, solid, &[e1, e2], 0.3).unwrap();
        verify_euler(&topo, chamfered);
        verify_all_twins(&topo, chamfered);
        assert_eq!(face_count(&topo, chamfered), 8, "Chamfer 2 edges → 8 faces");
    }

    #[test]
    fn test_chamfer_one_edge_evolution() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (_chamfered, evo) = chamfer_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();

        // 7 faces in provenance: 6 original + 1 chamfer.
        assert_eq!(evo.face_provenance.len(), 7, "Should track 7 faces");

        // Count provenance categories.
        let modified = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::Modified(_)))
            .count();
        let copied = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::CopiedFrom(_)))
            .count();
        let from_edge = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::FromEdge(_)))
            .count();

        assert_eq!(modified, 2, "2 faces adjacent to chamfered edge → Modified");
        assert_eq!(copied, 4, "4 faces not adjacent → CopiedFrom");
        assert_eq!(from_edge, 1, "1 chamfer face → FromEdge");

        // 1 deleted edge (the chamfered edge).
        assert_eq!(evo.deleted_edges.len(), 1);
        assert_eq!(evo.deleted_edges[0], edge);
    }
}

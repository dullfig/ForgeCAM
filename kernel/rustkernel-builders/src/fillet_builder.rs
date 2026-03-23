use std::collections::HashMap;
use std::collections::HashSet;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::evolution::{EdgeOrigin, FaceOrigin, ShapeEvolution, VertexOrigin};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{AnalyticalGeomStore, CylinderSurface, SurfaceDef};

use crate::chamfer_builder::point_hash;
use crate::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};
use crate::edge_analysis::{
    dihedral_angle, edge_adjacency, edge_convexity, edge_endpoints, face_centroid, plane_normal,
    EdgeConvexity,
};

/// Errors from the fillet builder.
#[derive(Debug)]
pub enum FilletError {
    EdgeAnalysis(crate::edge_analysis::EdgeAnalysisError),
    NotConvex(EdgeIdx),
    SharedVertex(VertexIdx),
}

impl std::fmt::Display for FilletError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FilletError::EdgeAnalysis(e) => write!(f, "Edge analysis error: {e}"),
            FilletError::NotConvex(e) => write!(f, "Edge {:?} is not convex", e),
            FilletError::SharedVertex(v) => {
                write!(f, "Vertex {:?} is shared by multiple filleted edges", v)
            }
        }
    }
}

impl std::error::Error for FilletError {}

impl From<crate::edge_analysis::EdgeAnalysisError> for FilletError {
    fn from(e: crate::edge_analysis::EdgeAnalysisError) -> Self {
        FilletError::EdgeAnalysis(e)
    }
}

/// Computed contact information for one filleted edge.
struct FilletContact {
    edge: EdgeIdx,
    vert_a: VertexIdx,
    vert_b: VertexIdx,
    face_a: FaceIdx,
    face_b: FaceIdx,
    /// Contact point on face_a at endpoint A.
    c1_a: Point3,
    /// Contact point on face_a at endpoint B.
    c1_b: Point3,
    /// Contact point on face_b at endpoint A.
    c2_a: Point3,
    /// Contact point on face_b at endpoint B.
    c2_b: Point3,
    /// Arc intermediate points at endpoint A (not including c1_a and c2_a).
    arc_a: Vec<Point3>,
    /// Arc intermediate points at endpoint B (not including c1_b and c2_b).
    arc_b: Vec<Point3>,
    /// Fillet cylinder surface.
    surface: SurfaceDef,
}

/// Compute the inward direction from an edge toward a face's interior,
/// projected perpendicular to the edge direction.
fn compute_inward(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    face_idx: FaceIdx,
    edge_midpoint: &Point3,
    edge_dir: &Vec3,
) -> Vec3 {
    let centroid = face_centroid(topo, geom, face_idx);
    let toward_face = centroid - edge_midpoint;
    let perp = toward_face - edge_dir.scale(toward_face.dot(edge_dir));
    let len = perp.norm();
    if len < 1e-15 {
        let n = plane_normal(geom, topo, face_idx);
        let fallback = n.cross(edge_dir);
        let fl = fallback.norm();
        if fl > 1e-15 { fallback / fl } else { Vec3::zeros() }
    } else {
        perp / len
    }
}

/// Generate arc points from `start` to `end` by rotating `(start - center)` around `axis`.
/// Returns `n_segments - 1` intermediate points (not including start and end).
fn generate_arc_points(
    center: &Point3,
    start: &Point3,
    axis: &Vec3,
    subtended_angle: f64,
    n_segments: usize,
) -> Vec<Point3> {
    let mut result = Vec::with_capacity(n_segments.saturating_sub(1));
    let r = start - center;

    for i in 1..n_segments {
        let t = i as f64 / n_segments as f64;
        let angle = subtended_angle * t;
        // Rodrigues rotation: r_rot = r*cos(θ) + (axis × r)*sin(θ) + axis*(axis·r)*(1-cos(θ))
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let r_rot = r.scale(cos_a)
            + axis.cross(&r).scale(sin_a)
            + axis.scale(axis.dot(&r) * (1.0 - cos_a));
        result.push(center + r_rot);
    }
    result
}

/// Build fillets on the specified edges of a solid.
///
/// Returns a new SolidIdx with fillet (cylindrical) faces inserted. The original solid's
/// topology is not modified (rebuild approach).
///
/// # Constraints
/// - All edges must be convex and have planar adjacent faces.
/// - No two edges may share a vertex.
/// - `arc_segments` controls the number of segments in the fillet arc (minimum 2).
pub fn fillet_edges_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    edges: &[EdgeIdx],
    radius: f64,
    arc_segments: usize,
) -> Result<(SolidIdx, ShapeEvolution), FilletError> {
    let _span = info_span!("fillet_edges", n_edges = edges.len(), radius, arc_segments).entered();
    let arc_segments = arc_segments.max(2);

    // --- Phase 1: Validate and compute contacts ---
    let mut all_fillet_verts: HashMap<VertexIdx, EdgeIdx> = HashMap::new();
    let mut contacts: Vec<FilletContact> = Vec::with_capacity(edges.len());

    for &edge_idx in edges {
        let conv = edge_convexity(topo, geom, edge_idx)?;
        if conv != EdgeConvexity::Convex {
            return Err(FilletError::NotConvex(edge_idx));
        }

        let adj = edge_adjacency(topo, edge_idx)?;
        let (pt_a, pt_b) = edge_endpoints(topo, geom, edge_idx);
        let alpha = dihedral_angle(topo, geom, edge_idx)?;

        let he_a_data = topo.half_edges.get(adj.he_a);
        let he_b_data = topo.half_edges.get(adj.he_b);
        let vert_a = he_a_data.origin;
        let vert_b = he_b_data.origin;

        for &v in &[vert_a, vert_b] {
            if let Some(&prev_edge) = all_fillet_verts.get(&v) {
                if prev_edge != edge_idx {
                    return Err(FilletError::SharedVertex(v));
                }
            }
            all_fillet_verts.insert(v, edge_idx);
        }

        // Offset distance: d = R / tan(α/2)
        let d = radius / (alpha / 2.0).tan();

        let edge_vec = pt_b - pt_a;
        let edge_len = edge_vec.norm();
        let edge_dir = if edge_len > 1e-15 { edge_vec / edge_len } else { Vec3::zeros() };
        let mid = Point3::from((pt_a.coords + pt_b.coords) * 0.5);

        let inward_1 = compute_inward(topo, geom, adj.face_a, &mid, &edge_dir);
        let inward_2 = compute_inward(topo, geom, adj.face_b, &mid, &edge_dir);

        let c1_a = pt_a + d * inward_1;
        let c1_b = pt_b + d * inward_1;
        let c2_a = pt_a + d * inward_2;
        let c2_b = pt_b + d * inward_2;

        // Fillet center at endpoint A: center = A + d * inward_1 - R * n1
        // This puts the center at distance R from both face planes.
        let n1 = plane_normal(geom, topo, adj.face_a);
        let center_a = pt_a + d * inward_1 - radius * n1;

        // The arc subtends angle (π - α) from c1 to c2 around the edge direction.
        let subtended = std::f64::consts::PI - alpha;

        // Generate arc points at endpoints A and B.
        let arc_a = generate_arc_points(&center_a, &c1_a, &edge_dir, subtended, arc_segments);
        let center_b = pt_b + d * inward_1 - radius * n1;
        let arc_b = generate_arc_points(&center_b, &c1_b, &edge_dir, subtended, arc_segments);

        // Fillet surface: cylinder along the edge through the centers.
        let surface = SurfaceDef::Cylinder(CylinderSurface {
            origin: center_a,
            axis: edge_dir,
            radius,
        });

        contacts.push(FilletContact {
            edge: edge_idx,
            vert_a,
            vert_b,
            face_a: adj.face_a,
            face_b: adj.face_b,
            c1_a,
            c1_b,
            c2_a,
            c2_b,
            arc_a,
            arc_b,
            surface,
        });
    }

    // --- Phase 2: Build vertex replacement map ---
    let mut vertex_face_replacements: HashMap<(u32, u32), Vec<Point3>> = HashMap::new();
    let mut fillet_faces: Vec<(Vec<Point3>, SurfaceDef)> = Vec::new();

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

        // Fillet face polygon:
        // [c1_A, c1_B, arc_B[last], ..., arc_B[0], c2_B, c2_A, arc_A[0], ..., arc_A[last]]
        // Wait — the winding needs to be outward. Let me construct it carefully.
        //
        // The polygon goes along face_a side (c1_A → c1_B), then down through arc at B
        // (c1_B → arc_B → c2_B), then along face_b side (c2_B → c2_A), then up through arc at A
        // (c2_A → arc_A_reversed → c1_A).
        let mut fillet_poly = Vec::new();
        fillet_poly.push(contact.c1_a);
        fillet_poly.push(contact.c1_b);
        // Arc at B goes from c1_b → c2_b.
        for p in &contact.arc_b {
            fillet_poly.push(*p);
        }
        fillet_poly.push(contact.c2_b);
        fillet_poly.push(contact.c2_a);
        // Arc at A goes from c2_a back to c1_a (reverse of c1_a → c2_a).
        for p in contact.arc_a.iter().rev() {
            fillet_poly.push(*p);
        }

        // Check winding: the normal of this polygon should be consistent with outward direction.
        let n_a = plane_normal(geom, topo, contact.face_a);
        let n_b = plane_normal(geom, topo, contact.face_b);
        let expected_dir = (n_a + n_b).normalize();

        // Compute polygon normal from first three points.
        if fillet_poly.len() >= 3 {
            let e1 = fillet_poly[1] - fillet_poly[0];
            let e2 = fillet_poly[2] - fillet_poly[0];
            let poly_n = e1.cross(&e2);
            if poly_n.dot(&expected_dir) < 0.0 {
                fillet_poly.reverse();
            }
        }

        fillet_faces.push((fillet_poly, contact.surface.clone()));
    }

    // --- Phase 3: Handle end faces ---
    let shell_idx = topo.solids.get(solid).outer_shell();
    let face_list: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

    for &face_idx in &face_list {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let vert = topo.half_edges.get(he).origin;
            if all_fillet_verts.contains_key(&vert) {
                let key = (vert.raw(), face_idx.raw());
                if !vertex_face_replacements.contains_key(&key) {
                    let &edge_idx = all_fillet_verts.get(&vert).unwrap();
                    let contact = contacts.iter().find(|c| c.edge == edge_idx).unwrap();

                    let (c1, arc, c2) = if vert == contact.vert_a {
                        (contact.c1_a, &contact.arc_a, contact.c2_a)
                    } else {
                        (contact.c1_b, &contact.arc_b, contact.c2_b)
                    };

                    // Determine ordering via the previous half-edge's twin face.
                    let prev_he = find_prev_he(topo, loop_idx, he);
                    let prev_twin_face = topo.half_edges.get(prev_he).twin
                        .map(|twin| {
                            let twin_loop = topo.half_edges.get(twin).loop_ref;
                            topo.loops.get(twin_loop).face
                        });

                    let ordered = if prev_twin_face == Some(contact.face_a) {
                        // Coming from face_a side: c1 first, then arc, then c2.
                        let mut pts = vec![c1];
                        pts.extend_from_slice(arc);
                        pts.push(c2);
                        pts
                    } else {
                        // Coming from face_b side: c2 first, then arc reversed, then c1.
                        let mut pts = vec![c2];
                        pts.extend(arc.iter().rev().copied());
                        pts.push(c1);
                        pts
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

    // Collect faces adjacent to filleted edges.
    let mut adjacent_faces: HashSet<u32> = HashSet::new();
    for contact in &contacts {
        adjacent_faces.insert(contact.face_a.raw());
        adjacent_faces.insert(contact.face_b.raw());
    }

    // Record deleted edges.
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

    // Rebuild each original face.
    for &face_idx in &face_list {
        let surface_id = topo.faces.get(face_idx).surface_id;
        let new_surface_id = geom.add_surface(geom.surfaces[surface_id as usize].clone());

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

        let vert_idxs: Vec<VertexIdx> = new_polygon
            .iter()
            .map(|pt| get_or_create_vert(geom, topo, *pt, &mut point_to_vert))
            .collect();

        let new_face = build_face_from_vert_idxs(topo, geom, &vert_idxs, new_surface_id, new_shell_idx, &mut he_map);
        if adjacent_faces.contains(&face_idx.raw()) {
            evo.record_face(new_face, FaceOrigin::Modified(face_idx));
        } else {
            evo.record_face(new_face, FaceOrigin::CopiedFrom(face_idx));
        }
    }

    // Add fillet faces.
    for (i, (poly, surface)) in fillet_faces.into_iter().enumerate() {
        let new_surface_id = geom.add_surface(surface);
        let vert_idxs: Vec<VertexIdx> = poly
            .iter()
            .map(|pt| get_or_create_vert(geom, topo, *pt, &mut point_to_vert))
            .collect();
        let new_face = build_face_from_vert_idxs(topo, geom, &vert_idxs, new_surface_id, new_shell_idx, &mut he_map);
        evo.record_face(new_face, FaceOrigin::FromEdge(contacts[i].edge));
    }

    match_twins_from_map(topo, &he_map);

    // Edge and vertex provenance — rebuild creates entirely new topology.
    let (_, new_edges, new_verts) = crate::euler_chamfer::collect_solid_entities(topo, new_solid_idx);
    for edge in new_edges {
        evo.record_edge(edge, EdgeOrigin::Primitive);
    }
    for vert in new_verts {
        evo.record_vertex(vert, VertexOrigin::Primitive);
    }

    Ok((new_solid_idx, evo))
}

/// Find the half-edge before `target_he` in the loop.
/// Delegates to the canonical implementation in `euler`.
fn find_prev_he(topo: &TopoStore, loop_idx: LoopIdx, target_he: HalfEdgeIdx) -> HalfEdgeIdx {
    rustkernel_topology::euler::find_prev_he(topo, loop_idx, target_he)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box_builder::make_box_into;
    use crate::edge_analysis::solid_edges;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

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
    fn test_fillet_one_edge_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (filleted, _evo) = fillet_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();
        verify_euler(&topo, filleted);
    }

    #[test]
    fn test_fillet_one_edge_twins() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (filleted, _evo) = fillet_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();
        verify_all_twins(&topo, filleted);
    }

    #[test]
    fn test_fillet_one_edge_face_count() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (filleted, _evo) = fillet_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();
        assert_eq!(face_count(&topo, filleted), 7, "Fillet 1 edge → 7 faces");
    }

    #[test]
    fn test_fillet_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (filleted, _evo) = fillet_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();
        let shell = topo.solids.get(filleted).outer_shell();
        tessellate_shell(&mut topo, shell, &geom);

        let mut total_tris = 0;
        for &face_idx in &topo.shells.get(shell).faces.clone() {
            let mesh = topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert!(total_tris > 0, "Should have triangles");
    }

    #[test]
    fn test_fillet_two_opposite_edges() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

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

        let (filleted, _evo) = fillet_edges_into(&mut topo, &mut geom, solid, &[e1, e2], 0.3, 4).unwrap();
        verify_euler(&topo, filleted);
        verify_all_twins(&topo, filleted);
        assert_eq!(face_count(&topo, filleted), 8, "Fillet 2 edges → 8 faces");
    }

    #[test]
    fn test_fillet_one_edge_evolution() {
        use rustkernel_topology::evolution::FaceOrigin;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (_filleted, evo) = fillet_edges_into(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();

        // 7 faces: 6 original + 1 fillet.
        assert_eq!(evo.face_provenance.len(), 7, "Should track 7 faces");

        let modified = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::Modified(_)))
            .count();
        let copied = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::CopiedFrom(_)))
            .count();
        let from_edge = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::FromEdge(_)))
            .count();

        assert_eq!(modified, 2, "2 faces adjacent to filleted edge → Modified");
        assert_eq!(copied, 4, "4 faces not adjacent → CopiedFrom");
        assert_eq!(from_edge, 1, "1 fillet face → FromEdge");

        assert_eq!(evo.deleted_edges.len(), 1);
        assert_eq!(evo.deleted_edges[0], edge);
    }
}

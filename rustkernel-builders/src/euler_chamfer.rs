//! Euler-operator-based chamfer builder.
//!
//! Operates in-place on the solid's topology using a 9-op sequence per edge:
//! 4× SEMV + 2× MEF + 1× KEF + 2× KEV.
//!
//! Scope: straight edges between planar faces, convex only, no shared vertices.

use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::euler::{
    find_he_from_vertex, find_prev_he, kef, kev, mef, semv, EulerError,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{AnalyticalGeomStore, CurveDef, LineSegment, Plane, SurfaceDef};

use crate::edge_analysis::{
    edge_adjacency, edge_convexity, edge_endpoints, face_centroid, plane_normal, EdgeConvexity,
};

// ─── Error type ──────────────────────────────────────────────────────────────

#[derive(Debug)]
pub enum EulerChamferError {
    EdgeAnalysis(crate::edge_analysis::EdgeAnalysisError),
    NotConvex(EdgeIdx),
    SharedVertex(VertexIdx),
    Euler(EulerError),
    /// A half-edge could not be found in a face loop.
    HalfEdgeNotFound { face: FaceIdx, vertex: VertexIdx },
}

impl std::fmt::Display for EulerChamferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EulerChamferError::EdgeAnalysis(e) => write!(f, "Edge analysis error: {e}"),
            EulerChamferError::NotConvex(e) => write!(f, "Edge {:?} is not convex", e),
            EulerChamferError::SharedVertex(v) => {
                write!(f, "Vertex {:?} is shared by multiple chamfered edges", v)
            }
            EulerChamferError::Euler(e) => write!(f, "Euler operator error: {e}"),
            EulerChamferError::HalfEdgeNotFound { face, vertex } => {
                write!(
                    f,
                    "No half-edge from vertex {:?} in face {:?}",
                    vertex, face
                )
            }
        }
    }
}

impl std::error::Error for EulerChamferError {}

impl From<crate::edge_analysis::EdgeAnalysisError> for EulerChamferError {
    fn from(e: crate::edge_analysis::EdgeAnalysisError) -> Self {
        EulerChamferError::EdgeAnalysis(e)
    }
}

impl From<EulerError> for EulerChamferError {
    fn from(e: EulerError) -> Self {
        EulerChamferError::Euler(e)
    }
}

// ─── Contact geometry ────────────────────────────────────────────────────────

/// Computed contact information for one chamfered edge.
struct ChamferContact {
    edge: EdgeIdx,
    vert_a: VertexIdx,
    vert_b: VertexIdx,
    face_a: FaceIdx,
    face_b: FaceIdx,
    he_a: HalfEdgeIdx,
    he_b: HalfEdgeIdx,
    c1_a: Point3,
    c1_b: Point3,
    c2_a: Point3,
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
    let perp = toward_face - edge_dir.scale(toward_face.dot(edge_dir));
    let len = perp.norm();
    if len < 1e-15 {
        let n = plane_normal(geom, topo, face_idx);
        let fallback = n.cross(edge_dir);
        let fl = fallback.norm();
        if fl > 1e-15 {
            fallback / fl
        } else {
            Vec3::zeros()
        }
    } else {
        perp / len
    }
}

/// Compute contact info for one edge.
fn compute_contact(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
    distance: f64,
) -> Result<ChamferContact, EulerChamferError> {
    let conv = edge_convexity(topo, geom, edge_idx)?;
    if conv != EdgeConvexity::Convex {
        return Err(EulerChamferError::NotConvex(edge_idx));
    }

    let adj = edge_adjacency(topo, edge_idx)?;
    let (pt_a, pt_b) = edge_endpoints(topo, geom, edge_idx);

    let vert_a = topo.half_edges.get(adj.he_a).origin;
    let vert_b = topo.half_edges.get(adj.he_b).origin;

    let edge_vec = pt_b - pt_a;
    let edge_len = edge_vec.norm();
    let edge_dir = if edge_len > 1e-15 {
        edge_vec / edge_len
    } else {
        Vec3::zeros()
    };
    let mid = Point3::from((pt_a.coords + pt_b.coords) * 0.5);

    let inward_1 = compute_inward(topo, geom, adj.face_a, &mid, &edge_dir);
    let inward_2 = compute_inward(topo, geom, adj.face_b, &mid, &edge_dir);

    Ok(ChamferContact {
        edge: edge_idx,
        vert_a,
        vert_b,
        face_a: adj.face_a,
        face_b: adj.face_b,
        he_a: adj.he_a,
        he_b: adj.he_b,
        c1_a: pt_a + distance * inward_1,
        c1_b: pt_b + distance * inward_1,
        c2_a: pt_a + distance * inward_2,
        c2_b: pt_b + distance * inward_2,
    })
}

// ─── Single-edge Euler chamfer ───────────────────────────────────────────────

/// Apply a single-edge chamfer using the 9-op Euler sequence.
///
/// Returns the FaceIdx of the new chamfer face (quad).
fn euler_chamfer_single_edge(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    contact: &ChamferContact,
) -> Result<FaceIdx, EulerChamferError> {
    let _span = info_span!("euler_chamfer_single_edge", edge = contact.edge.raw()).entered();

    // Alias key topology from the contact
    let face_a = contact.face_a;
    let face_b = contact.face_b;
    let vert_a = contact.vert_a;
    let vert_b = contact.vert_b;

    // ── Phase A: 4× SEMV — Create contact vertices on adjacent edges ──

    // We need to find the half-edges around V_a and V_b in both faces.
    // h_a goes V_a → V_b in face_a. So:
    //   h_prev_a: ...→V_a (the HE ending at V_a in face_a, i.e. the HE before h_a)
    //   h_next_a: V_b→... (the HE after h_a)
    // h_b goes V_b → V_a in face_b. So:
    //   h_prev_b: ...→V_b (the HE before h_b)
    //   h_next_b: V_a→... (the HE after h_b)

    let h_a = contact.he_a; // V_a → V_b in face_a
    let h_b = contact.he_b; // V_b → V_a in face_b

    let loop_a = topo.half_edges.get(h_a).loop_ref;
    let loop_b = topo.half_edges.get(h_b).loop_ref;

    // Find predecessors and successors
    let h_prev_a = find_prev_he(topo, loop_a, h_a); // ...→V_a in face_a
    let h_next_a = topo.half_edges.get(h_a).next; // V_b→... in face_a

    let h_prev_b = find_prev_he(topo, loop_b, h_b); // ...→V_b in face_b
    let h_next_b = topo.half_edges.get(h_b).next; // V_a→... in face_b

    // Register contact point geometry
    let pid_c1a = geom.add_point(contact.c1_a);
    let pid_c1b = geom.add_point(contact.c1_b);
    let pid_c2a = geom.add_point(contact.c2_a);
    let pid_c2b = geom.add_point(contact.c2_b);

    // For SEMV curve_ids, we create line segments for the split edges.
    // The original edge's endpoints are known from geom.
    let pt_a = geom.points[topo.vertices.get(vert_a).point_id as usize];
    let pt_b = geom.points[topo.vertices.get(vert_b).point_id as usize];

    // 1. SEMV(h_prev_a): splits edge ending at V_a in face_a → inserts c1_a
    //    h_prev_a goes from some V_prev to V_a. Split creates c1_a between V_prev and V_a.
    let prev_a_origin_pid = topo.vertices.get(topo.half_edges.get(h_prev_a).origin).point_id;
    let prev_a_origin_pt = geom.points[prev_a_origin_pid as usize];
    let curve_prev_a_1 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: prev_a_origin_pt,
        end: contact.c1_a,
    }));
    let curve_prev_a_2 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c1_a,
        end: pt_a,
    }));
    let semv1 = semv(topo, h_prev_a, pid_c1a, curve_prev_a_1, curve_prev_a_2);
    // c1_a vertex created. face_a loop now has: ...→ c1_a → V_a → V_b → ...

    // 2. SEMV(h_next_a): splits edge starting at V_b in face_a → inserts c1_b
    //    h_next_a goes V_b → V_next. Split creates c1_b between V_b and V_next.
    //    After SEMV1, h_next_a may still be valid (it's on a different edge).
    //    But we must re-fetch h_next_a after mutations. Actually h_next_a was read before
    //    mutations and points at the original HE. SEMV1 split h_prev_a's edge, so h_next_a is unaffected.
    let _next_a_dest_he = topo.half_edges.get(h_next_a).next;
    // Actually for h_next_a, origin=V_b, dest=V_next. We need V_next's point.
    // The destination of h_next_a is the origin of h_next_a.next (the HE after h_next_a).
    // But after SEMV1 this chain may have changed. Let's just read the endpoints of h_next_a's edge.
    let h_next_a_origin_pid = topo.vertices.get(topo.half_edges.get(h_next_a).origin).point_id;
    let h_next_a_origin_pt = geom.points[h_next_a_origin_pid as usize];
    // Destination of h_next_a: walk to next HE in the loop to get its origin.
    let h_next_a_next = topo.half_edges.get(h_next_a).next;
    let h_next_a_dest_pid = topo.vertices.get(topo.half_edges.get(h_next_a_next).origin).point_id;
    let h_next_a_dest_pt = geom.points[h_next_a_dest_pid as usize];
    let curve_next_a_1 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: h_next_a_origin_pt,
        end: contact.c1_b,
    }));
    let curve_next_a_2 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c1_b,
        end: h_next_a_dest_pt,
    }));
    let semv2 = semv(topo, h_next_a, pid_c1b, curve_next_a_1, curve_next_a_2);
    // c1_b vertex created. face_a loop now: ...→ c1_a → V_a → V_b → c1_b → ...

    // 3. SEMV(h_next_b): splits edge starting at V_a in face_b → inserts c2_a
    //    h_next_b goes V_a → V_next_b. Split creates c2_a between V_a and V_next_b.
    let h_next_b_origin_pid = topo.vertices.get(topo.half_edges.get(h_next_b).origin).point_id;
    let h_next_b_origin_pt = geom.points[h_next_b_origin_pid as usize];
    let h_next_b_next = topo.half_edges.get(h_next_b).next;
    let h_next_b_dest_pid =
        topo.vertices.get(topo.half_edges.get(h_next_b_next).origin).point_id;
    let h_next_b_dest_pt = geom.points[h_next_b_dest_pid as usize];
    let curve_next_b_1 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: h_next_b_origin_pt,
        end: contact.c2_a,
    }));
    let curve_next_b_2 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c2_a,
        end: h_next_b_dest_pt,
    }));
    let semv3 = semv(topo, h_next_b, pid_c2a, curve_next_b_1, curve_next_b_2);
    // c2_a vertex created. face_b loop: ...→ c2_b? no → V_b → V_a → c2_a → ...

    // 4. SEMV(h_prev_b): splits edge ending at V_b in face_b → inserts c2_b
    //    h_prev_b goes V_prev_b → V_b. Split creates c2_b between V_prev_b and V_b.
    let prev_b_origin_pid = topo.vertices.get(topo.half_edges.get(h_prev_b).origin).point_id;
    let prev_b_origin_pt = geom.points[prev_b_origin_pid as usize];
    let curve_prev_b_1 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: prev_b_origin_pt,
        end: contact.c2_b,
    }));
    let curve_prev_b_2 = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c2_b,
        end: pt_b,
    }));
    let semv4 = semv(topo, h_prev_b, pid_c2b, curve_prev_b_1, curve_prev_b_2);
    // c2_b vertex created. face_b loop: ...→ c2_b → V_b → V_a → c2_a → ...

    // ── Phase B: 2× MEF — Split faces along contact lines ──

    // After the 4 SEMVs:
    // face_a loop: ...→ c1_a → V_a → V_b → c1_b → ...
    // face_b loop: ...→ c2_b → V_b → V_a → c2_a → ...
    //
    // MEF(he@c1_a_in_face_a, he@c1_b_in_face_a) splits face_a:
    //   strip_a = [c1_a, V_a, V_b, c1_b]
    //   trimmed face_a = remainder
    //
    // We need to find the HE with origin=c1_a and origin=c1_b in face_a's current loop.
    // After mutations, we must re-walk.

    let v_c1a = semv1.new_vertex;
    let v_c1b = semv2.new_vertex;
    let v_c2a = semv3.new_vertex;
    let v_c2b = semv4.new_vertex;

    // Find HEs in face_a for c1_a and c1_b
    // After SEMV1, there's a HE from c1_a in face_a (the new_he_a from semv1 or the
    // continuation). Let's walk the loop to find them.
    let he_c1a_in_fa = find_he_from_vertex(topo, face_a, v_c1a)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_a,
            vertex: v_c1a,
        })?;

    // For MEF, we need h1=HE@c1_a and h2=HE@c1_b.
    // But we need to find the HE that starts the strip c1_b → ... (the one AFTER the strip).
    // Actually MEF(h1, h2) connects h1.origin to h2.origin.
    // We want to cut face_a along the line c1_a ↔ c1_b.
    // The strip that contains V_a and V_b should go: c1_a → V_a → V_b → c1_b.
    // We need h1 = HE@c1_b (after c1_b in the loop) and h2 = HE@c1_a.
    // Wait, let me think about MEF semantics:
    //
    // MEF(h1, h2): connects h1.origin (V1) and h2.origin (V2).
    // Old loop keeps: h1 → ... → prev(h2) → he_a(V2→V1) → h1
    // New loop gets: h2 → ... → prev(h1) → he_b(V1→V2) → h2
    //
    // face_a loop order: ... → c1_a → V_a → V_b → c1_b → rest → ...
    //
    // We want strip_a = new face = [c1_a, V_a, V_b, c1_b].
    // That means the new loop should contain: c1_a → V_a → V_b → c1_b → (back via he_b to c1_a).
    // So new loop = h2 → ... → prev(h1) → he_b → h2 where h2.origin = c1_a.
    // The path h2 → ... → prev(h1) = c1_a → V_a → V_b.
    // So h1.origin = c1_b (to have prev(h1) end at c1_b's predecessor = V_b's HE to c1_b).
    // Wait, prev(h1) should be the HE from V_b to c1_b. And h1 starts at c1_b.
    //
    // Let me re-read the MEF code. After MEF(h1=he@c1_b, h2=he@c1_a):
    //   Old loop (face_a remains): c1_b → rest → ... → (back to prev(c1_a)) → he_a(c1_a→c1_b) → c1_b
    //   New loop (strip_a):        c1_a → V_a → V_b → (reaches prev(c1_b)) → he_b(c1_b→c1_a) → c1_a
    //
    // That gives strip_a = [c1_a, V_a, V_b, he_b_back_to_c1_a]. This is the strip containing
    // V_a and V_b between the contact points. That's what we want.

    let he_c1b_in_fa = find_he_from_vertex(topo, face_a, v_c1b)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_a,
            vertex: v_c1b,
        })?;

    // Curve for the MEF edge (c1_a ↔ c1_b line)
    let curve_mef_a = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c1_a,
        end: contact.c1_b,
    }));

    // The surface_id for the NEW face (strip_a) — use face_a's surface since it's the same plane
    let surface_a_id = topo.faces.get(face_a).surface_id;

    let mef1 = mef(topo, he_c1b_in_fa, he_c1a_in_fa, surface_a_id, curve_mef_a)?;
    // mef1.new_face = strip_a (contains c1_a, V_a, V_b, c1_b)
    let strip_a = mef1.new_face;

    // Now do the same for face_b.
    // face_b loop order after SEMVs: ... → c2_b → V_b → V_a → c2_a → rest → ...
    // We want strip_b = new face = [c2_b, V_b, V_a, c2_a].
    // By the same logic: MEF(he@c2_a, he@c2_b) creates:
    //   New loop (strip_b): c2_b → V_b → V_a → he_b(c2_a→c2_b) → c2_b

    let he_c2a_in_fb = find_he_from_vertex(topo, face_b, v_c2a)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_b,
            vertex: v_c2a,
        })?;
    let he_c2b_in_fb = find_he_from_vertex(topo, face_b, v_c2b)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_b,
            vertex: v_c2b,
        })?;

    let curve_mef_b = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c2_b,
        end: contact.c2_a,
    }));

    let surface_b_id = topo.faces.get(face_b).surface_id;

    let mef2 = mef(topo, he_c2a_in_fb, he_c2b_in_fb, surface_b_id, curve_mef_b)?;
    let _strip_b = mef2.new_face;

    // ── Phase C: 1× KEF — Merge strips by removing edge E ──
    //
    // Edge E (the original chamfered edge) is shared between strip_a and strip_b.
    // KEF(E, strip_a) merges strip_b into strip_a, creating a hexagon:
    //   [c1_a, V_a, c2_a, c2_b, V_b, c1_b]

    kef(topo, contact.edge, strip_a)?;
    // strip_a is now the hexagon. strip_b is dead.
    let chamfer_face = strip_a;

    // ── Phase D: 2× KEV — Remove corner vertices ──
    //
    // V_a has edges: (c1_a ↔ V_a) and (V_a ↔ c2_a) in the hexagon, plus possibly
    // connections on end faces. But after the MEFs, V_a only connects to c1_a and c2_a
    // within the chamfer face. However, V_a is also on the end faces (the faces at
    // each end of the original edge, perpendicular to it).
    //
    // The edge to kill for V_a: the edge between c1_a and V_a (from SEMV1's leftover).
    // After SEMV1, the original edge (h_prev_a's edge) became V_prev → c1_a, and
    // the new edge (semv1.new_edge) became c1_a → V_a.
    // Similarly for V_b with SEMV4: semv4.new_edge is c2_b → V_b.
    //
    // Wait — we need to find the edges connecting V_a to its neighbors c1_a and c2_a.
    // After KEF, the hexagon has edges: c1_a-V_a, V_a-c2_a, c2_a-c2_b, c2_b-V_b, V_b-c1_b, c1_b-c1_a.
    // Plus the MEF edges. Let me think about which edge connects V_a to c1_a.
    //
    // From SEMV1: split h_prev_a's edge. The new edge (semv1.new_edge) goes c1_a → V_a.
    // But that edge is between the end face and face_a (now strip_a/chamfer). After MEF + KEF
    // the strip absorbed the V_a vertex. The edge c1_a↔V_a is shared between the chamfer
    // face and an end face.
    //
    // For KEV(edge, v_keep=c1_a): removes V_a, absorbs it into c1_a.
    // This removes the c1_a↔V_a edge from both the chamfer face and the end face.

    // Find the edge between c1_a and V_a. After SEMV1, semv1.new_edge connects c1_a to V_a.
    // Actually, SEMV splits h_prev_a. After split:
    //   - original edge (h_prev_a's edge): V_prev → c1_a
    //   - new_edge: c1_a → V_a
    // The new_edge's half_edges connect c1_a and V_a. That's what we need.
    let edge_c1a_va = semv1.new_edge;
    kev(topo, edge_c1a_va, v_c1a)?;

    // Similarly, find the edge between c2_b and V_b.
    // SEMV4 split h_prev_b. After split:
    //   - original edge: V_prev_b → c2_b
    //   - new_edge: c2_b → V_b
    let edge_c2b_vb = semv4.new_edge;
    kev(topo, edge_c2b_vb, v_c2b)?;

    // ── Phase E: Update geometry ──

    // Set chamfer face surface to a new plane through the contact points.
    let chamfer_normal = (contact.c1_b - contact.c1_a).cross(&(contact.c2_a - contact.c1_a));
    let cn_len = chamfer_normal.norm();
    let chamfer_normal_unit = if cn_len > 1e-15 {
        chamfer_normal / cn_len
    } else {
        Vec3::zeros()
    };

    // Ensure normal points outward (same direction as average of face_a and face_b normals)
    let n_a = plane_normal(geom, topo, face_a);
    let n_b = plane_normal(geom, topo, face_b);
    let expected_dir = (n_a + n_b).normalize();
    let final_normal = if chamfer_normal_unit.dot(&expected_dir) > 0.0 {
        chamfer_normal_unit
    } else {
        -chamfer_normal_unit
    };

    let chamfer_surface_id = geom.add_surface(SurfaceDef::Plane(Plane {
        origin: contact.c1_a,
        normal: final_normal,
    }));
    topo.faces.get_mut(chamfer_face).surface_id = chamfer_surface_id;

    // Update curve_ids on the chamfer face edges. Walk the loop and set line segments.
    let chamfer_loop = topo.faces.get(chamfer_face).outer_loop;
    let start = topo.loops.get(chamfer_loop).half_edge;
    let mut he = start;
    loop {
        let origin_vid = topo.half_edges.get(he).origin;
        let origin_pt = geom.points[topo.vertices.get(origin_vid).point_id as usize];
        let next_he = topo.half_edges.get(he).next;
        let dest_vid = topo.half_edges.get(next_he).origin;
        let dest_pt = geom.points[topo.vertices.get(dest_vid).point_id as usize];
        let curve_id = geom.add_curve(CurveDef::LineSegment(LineSegment {
            start: origin_pt,
            end: dest_pt,
        }));
        let edge_idx = topo.half_edges.get(he).edge;
        topo.edges.get_mut(edge_idx).curve_id = curve_id;
        he = next_he;
        if he == start {
            break;
        }
    }

    // Invalidate mesh_cache on all affected faces
    topo.faces.get_mut(face_a).mesh_cache = None;
    topo.faces.get_mut(face_b).mesh_cache = None;
    topo.faces.get_mut(chamfer_face).mesh_cache = None;

    // Also invalidate end faces (faces that share edges with V_a's and V_b's old positions).
    // Since KEV rewrites vertex origins across the shell, end faces automatically get
    // the new contact vertices. We invalidate all faces on the shell to be safe.
    let shell_idx = topo.faces.get(chamfer_face).shell;
    let all_faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for &f in &all_faces {
        topo.faces.get_mut(f).mesh_cache = None;
    }

    Ok(chamfer_face)
}

// ─── Public API ──────────────────────────────────────────────────────────────

/// Apply a constant-distance Euler-based chamfer to straight edges between planar faces.
///
/// Operates in-place on the solid (no rebuild). All edges must be convex and no two
/// edges may share a vertex.
///
/// Returns the same SolidIdx (modified in-place).
pub fn euler_chamfer_edges(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    edges: &[EdgeIdx],
    distance: f64,
) -> Result<SolidIdx, EulerChamferError> {
    let _span = info_span!("euler_chamfer_edges", n_edges = edges.len(), distance).entered();

    // Validate all edges and compute contacts up front
    let mut all_chamfer_verts: HashMap<VertexIdx, EdgeIdx> = HashMap::new();
    let mut contacts: Vec<ChamferContact> = Vec::with_capacity(edges.len());

    for &edge_idx in edges {
        let contact = compute_contact(topo, geom, edge_idx, distance)?;

        // SharedVertex check
        for &v in &[contact.vert_a, contact.vert_b] {
            if let Some(&prev_edge) = all_chamfer_verts.get(&v) {
                if prev_edge != edge_idx {
                    return Err(EulerChamferError::SharedVertex(v));
                }
            }
            all_chamfer_verts.insert(v, edge_idx);
        }

        contacts.push(contact);
    }

    // Process each edge sequentially
    for contact in &contacts {
        euler_chamfer_single_edge(topo, geom, contact)?;
    }

    Ok(solid)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box_builder::make_box_into;
    use crate::edge_analysis::{plane_normal, solid_edges};
    use rustkernel_topology::diagnostics::validate_solid;
    use std::collections::HashSet;

    /// Count V, E, F by walking the shell.
    fn count_vef(topo: &TopoStore, solid: SolidIdx) -> (usize, usize, usize) {
        let shell = topo.solids.get(solid).shell;
        let faces = &topo.shells.get(shell).faces;
        let mut verts = HashSet::new();
        let mut edges = HashSet::new();
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                verts.insert(topo.half_edges.get(he).origin.raw());
                edges.insert(topo.half_edges.get(he).edge.raw());
                he = topo.half_edges.get(he).next;
                if he == start {
                    break;
                }
            }
        }
        (verts.len(), edges.len(), faces.len())
    }

    fn verify_all_twins(topo: &TopoStore, solid: SolidIdx) {
        let shell = topo.solids.get(solid).shell;
        let faces = &topo.shells.get(shell).faces;
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                let he_data = topo.half_edges.get(he);
                assert!(
                    he_data.twin.is_some(),
                    "Half-edge {} has no twin",
                    he.raw()
                );
                let twin = he_data.twin.unwrap();
                assert_eq!(
                    topo.half_edges.get(twin).twin,
                    Some(he),
                    "Twin asymmetry at HE {}",
                    he.raw()
                );
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
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
    fn test_euler_chamfer_one_edge() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Before: box = 8V, 12E, 6F
        assert_eq!(count_vef(&topo, solid), (8, 12, 6));

        let edge = find_box_edge_between(
            &topo,
            &geom,
            solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let result = euler_chamfer_edges(&mut topo, &mut geom, solid, &[edge], 0.3);
        assert!(result.is_ok(), "Euler chamfer failed: {:?}", result.err());

        // After: V=10, E=15, F=7, Euler=2
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (10, 15, 7), "One chamfer: V=10, E=15, F=7");
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");

        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_euler_chamfer_one_edge_geometry() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        // Box from origin: corners at (0,0,0) and (2,2,2)
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Chamfer the bottom-front edge (between -Z and -Y faces)
        let edge = find_box_edge_between(
            &topo,
            &geom,
            solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let dist = 0.5;
        euler_chamfer_edges(&mut topo, &mut geom, solid, &[edge], dist).unwrap();

        // Collect all vertex positions
        let shell = topo.solids.get(solid).shell;
        let faces = &topo.shells.get(shell).faces;
        let mut all_pts: HashSet<[i64; 3]> = HashSet::new();
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                let vid = topo.half_edges.get(he).origin;
                let pid = topo.vertices.get(vid).point_id;
                let pt = geom.points[pid as usize];
                // Quantize to check positions
                let key = [
                    (pt.x * 100.0).round() as i64,
                    (pt.y * 100.0).round() as i64,
                    (pt.z * 100.0).round() as i64,
                ];
                all_pts.insert(key);
                he = topo.half_edges.get(he).next;
                if he == start {
                    break;
                }
            }
        }

        // Should have 10 unique vertices (8 original - 2 removed + 4 contact)
        assert_eq!(all_pts.len(), 10);

        // Original 6 vertices that were NOT on the chamfered edge should still be present
        // (the box has 8 vertices, 2 on the chamfered edge get removed, leaving 6 original + 4 new)
    }

    #[test]
    fn test_euler_chamfer_two_opposite_edges() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Pick two opposite edges: bottom-front and top-back (no shared vertices)
        let e1 = find_box_edge_between(
            &topo,
            &geom,
            solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );
        let e2 = find_box_edge_between(
            &topo,
            &geom,
            solid,
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
        );

        let result = euler_chamfer_edges(&mut topo, &mut geom, solid, &[e1, e2], 0.3);
        assert!(result.is_ok(), "Two-edge chamfer failed: {:?}", result.err());

        // After: V=12, E=18, F=8
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (12, 18, 8), "Two chamfers: V=12, E=18, F=8");
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");

        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_euler_chamfer_validate_solid() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo,
            &geom,
            solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        euler_chamfer_edges(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();

        let report = validate_solid(&topo, solid);
        assert!(
            report.is_valid(),
            "validate_solid failed: {:?}",
            report.errors()
        );
    }
}

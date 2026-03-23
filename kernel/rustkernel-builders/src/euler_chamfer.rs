//! Euler-operator-based chamfer builder.
//!
//! Operates in-place on the solid's topology using a 9-op sequence per edge:
//! 4× SEMV + 2× MEF + 1× KEF + 2× KEV.
//!
//! Supports chamfering multiple edges including edges meeting at a corner vertex.

use std::collections::{HashMap, HashSet};

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::euler::{
    find_he_from_vertex, find_prev_he, kef, kev, kev_same_loop, mef, semv, EulerError, SemvResult,
};
use rustkernel_topology::evolution::{
    EdgeOrigin, FaceOrigin, ShapeEvolution, VertexOrigin,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{AnalyticalGeomStore, CurveDef, LineSegment, Plane, SurfaceDef};

use crate::edge_analysis::{
    edge_adjacency, edge_convexity, edge_endpoints, face_centroid, plane_normal,
    surface_normal_at_point, EdgeConvexity,
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
pub(crate) struct ChamferContact {
    pub(crate) edge: EdgeIdx,
    pub(crate) vert_a: VertexIdx,
    pub(crate) vert_b: VertexIdx,
    pub(crate) face_a: FaceIdx,
    pub(crate) face_b: FaceIdx,
    pub(crate) he_a: HalfEdgeIdx,
    pub(crate) he_b: HalfEdgeIdx,
    pub(crate) c1_a: Point3,
    pub(crate) c1_b: Point3,
    pub(crate) c2_a: Point3,
    pub(crate) c2_b: Point3,
}

/// Geometry-only contact data (positions). Computed once, survives topology edits.
#[derive(Clone)]
pub(crate) struct ContactGeometry {
    pub(crate) edge: EdgeIdx,
    pub(crate) pt_a: Point3,
    pub(crate) pt_b: Point3,
    pub(crate) c1_a: Point3,
    pub(crate) c1_b: Point3,
    pub(crate) c2_a: Point3,
    pub(crate) c2_b: Point3,
}

/// Topology-dependent contact data (live HE/face/vertex queries).
pub(crate) struct ContactTopology {
    edge: EdgeIdx,
    vert_a: VertexIdx,
    vert_b: VertexIdx,
    face_a: FaceIdx,
    face_b: FaceIdx,
    he_a: HalfEdgeIdx,
    he_b: HalfEdgeIdx,
}

/// Compute the inward direction from an edge toward a face's interior,
/// projected perpendicular to the edge direction and lying in the face's plane.
pub(crate) fn compute_inward(
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
        let n = surface_normal_at_point(geom, topo, face_idx, edge_midpoint);
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

/// Compute geometry-only contact positions for one edge (idempotent, topology-independent).
pub(crate) fn compute_contact_geometry(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
    distance: f64,
) -> Result<ContactGeometry, EulerChamferError> {
    let conv = edge_convexity(topo, geom, edge_idx)?;
    if conv != EdgeConvexity::Convex {
        return Err(EulerChamferError::NotConvex(edge_idx));
    }

    let adj = edge_adjacency(topo, edge_idx)?;
    let (pt_a, pt_b) = edge_endpoints(topo, geom, edge_idx);

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

    Ok(ContactGeometry {
        edge: edge_idx,
        pt_a,
        pt_b,
        c1_a: pt_a + distance * inward_1,
        c1_b: pt_b + distance * inward_1,
        c2_a: pt_a + distance * inward_2,
        c2_b: pt_b + distance * inward_2,
    })
}

/// Query live topology for an edge (must be called fresh after any topology mutation).
pub(crate) fn query_contact_topology(
    topo: &TopoStore,
    edge_idx: EdgeIdx,
) -> Result<ContactTopology, EulerChamferError> {
    let adj = edge_adjacency(topo, edge_idx)?;
    let vert_a = topo.half_edges.get(adj.he_a).origin;
    let vert_b = topo.half_edges.get(adj.he_b).origin;
    Ok(ContactTopology {
        edge: edge_idx,
        vert_a,
        vert_b,
        face_a: adj.face_a,
        face_b: adj.face_b,
        he_a: adj.he_a,
        he_b: adj.he_b,
    })
}

/// Assemble a ChamferContact from pre-computed geometry and live topology.
pub(crate) fn assemble_contact(cg: &ContactGeometry, ct: &ContactTopology) -> ChamferContact {
    ChamferContact {
        edge: ct.edge,
        vert_a: ct.vert_a,
        vert_b: ct.vert_b,
        face_a: ct.face_a,
        face_b: ct.face_b,
        he_a: ct.he_a,
        he_b: ct.he_b,
        c1_a: cg.c1_a,
        c1_b: cg.c1_b,
        c2_a: cg.c2_a,
        c2_b: cg.c2_b,
    }
}

// ─── Edge classification ────────────────────────────────────────────────────

/// A group of edges sharing a common corner vertex.
pub(crate) struct CornerGroup {
    pub(crate) vertex: VertexIdx,
    pub(crate) edge_indices: Vec<usize>,
}

/// Classification of edges into independent edges and corner groups.
pub(crate) struct EdgeClassification {
    pub(crate) independent: Vec<usize>,
    pub(crate) corners: Vec<CornerGroup>,
}

/// Classify edges into independent (no shared vertices) and corner groups.
fn classify_edges(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edges: &[EdgeIdx],
    distance: f64,
) -> Result<(Vec<ContactGeometry>, EdgeClassification), EulerChamferError> {
    let mut geometries = Vec::with_capacity(edges.len());
    let mut vert_to_edges: HashMap<VertexIdx, Vec<usize>> = HashMap::new();

    for (i, &edge_idx) in edges.iter().enumerate() {
        let cg = compute_contact_geometry(topo, geom, edge_idx, distance)?;
        geometries.push(cg);

        let adj = edge_adjacency(topo, edge_idx)?;
        let va = topo.half_edges.get(adj.he_a).origin;
        let vb = topo.half_edges.get(adj.he_b).origin;
        vert_to_edges.entry(va).or_default().push(i);
        vert_to_edges.entry(vb).or_default().push(i);
    }

    // Edges at a shared vertex form a corner group
    let mut corner_edge_set = std::collections::HashSet::new();
    let mut corners = Vec::new();
    for (vertex, edge_list) in &vert_to_edges {
        if edge_list.len() >= 2 {
            for &idx in edge_list {
                corner_edge_set.insert(idx);
            }
            corners.push(CornerGroup {
                vertex: *vertex,
                edge_indices: edge_list.clone(),
            });
        }
    }

    let independent: Vec<usize> = (0..edges.len())
        .filter(|i| !corner_edge_set.contains(i))
        .collect();

    Ok((geometries, EdgeClassification { independent, corners }))
}

// ─── Phase result types ─────────────────────────────────────────────────────

pub(crate) struct SemvPhaseResult {
    pub(crate) v_c1a: VertexIdx,
    pub(crate) v_c1b: VertexIdx,
    pub(crate) v_c2a: VertexIdx,
    pub(crate) v_c2b: VertexIdx,
    pub(crate) semv1: SemvResult,
    pub(crate) semv4: SemvResult,
}

pub(crate) struct MefPhaseResult {
    pub(crate) strip_a: FaceIdx,
}

// ─── Phase sub-functions ────────────────────────────────────────────────────

/// Phase A: 4× SEMV — Create contact vertices on adjacent edges.
pub(crate) fn phase_a_semv(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    contact: &ChamferContact,
) -> SemvPhaseResult {
    let vert_a = contact.vert_a;
    let vert_b = contact.vert_b;
    let h_a = contact.he_a;
    let h_b = contact.he_b;

    let loop_a = topo.half_edges.get(h_a).loop_ref;
    let loop_b = topo.half_edges.get(h_b).loop_ref;

    let h_prev_a = find_prev_he(topo, loop_a, h_a);
    let h_next_a = topo.half_edges.get(h_a).next;
    let h_prev_b = find_prev_he(topo, loop_b, h_b);
    let h_next_b = topo.half_edges.get(h_b).next;

    let pid_c1a = geom.add_point(contact.c1_a);
    let pid_c1b = geom.add_point(contact.c1_b);
    let pid_c2a = geom.add_point(contact.c2_a);
    let pid_c2b = geom.add_point(contact.c2_b);

    let pt_a = geom.points[topo.vertices.get(vert_a).point_id as usize];
    let pt_b = geom.points[topo.vertices.get(vert_b).point_id as usize];

    // 1. SEMV(h_prev_a): splits edge ending at V_a in face_a → inserts c1_a
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

    // 2. SEMV(h_next_a): splits edge starting at V_b in face_a → inserts c1_b
    let h_next_a_origin_pid = topo.vertices.get(topo.half_edges.get(h_next_a).origin).point_id;
    let h_next_a_origin_pt = geom.points[h_next_a_origin_pid as usize];
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

    // 3. SEMV(h_next_b): splits edge starting at V_a in face_b → inserts c2_a
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

    // 4. SEMV(h_prev_b): splits edge ending at V_b in face_b → inserts c2_b
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

    SemvPhaseResult {
        v_c1a: semv1.new_vertex,
        v_c1b: semv2.new_vertex,
        v_c2a: semv3.new_vertex,
        v_c2b: semv4.new_vertex,
        semv1,
        semv4,
    }
}

/// Phase B: 2× MEF — Split faces along contact lines.
pub(crate) fn phase_b_mef(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    contact: &ChamferContact,
    semv_res: &SemvPhaseResult,
) -> Result<MefPhaseResult, EulerChamferError> {
    let face_a = contact.face_a;
    let face_b = contact.face_b;

    // Find HEs in face_a for c1_a and c1_b
    let he_c1a_in_fa = find_he_from_vertex(topo, face_a, semv_res.v_c1a)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_a,
            vertex: semv_res.v_c1a,
        })?;
    let he_c1b_in_fa = find_he_from_vertex(topo, face_a, semv_res.v_c1b)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_a,
            vertex: semv_res.v_c1b,
        })?;

    let curve_mef_a = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c1_a,
        end: contact.c1_b,
    }));
    let surface_a_id = topo.faces.get(face_a).surface_id;
    let mef1 = mef(topo, he_c1b_in_fa, he_c1a_in_fa, surface_a_id, curve_mef_a)?;
    let strip_a = mef1.new_face;

    // Find HEs in face_b for c2_a and c2_b
    let he_c2a_in_fb = find_he_from_vertex(topo, face_b, semv_res.v_c2a)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_b,
            vertex: semv_res.v_c2a,
        })?;
    let he_c2b_in_fb = find_he_from_vertex(topo, face_b, semv_res.v_c2b)
        .ok_or(EulerChamferError::HalfEdgeNotFound {
            face: face_b,
            vertex: semv_res.v_c2b,
        })?;

    let curve_mef_b = geom.add_curve(CurveDef::LineSegment(LineSegment {
        start: contact.c2_b,
        end: contact.c2_a,
    }));
    let surface_b_id = topo.faces.get(face_b).surface_id;
    mef(topo, he_c2a_in_fb, he_c2b_in_fb, surface_b_id, curve_mef_b)?;

    Ok(MefPhaseResult { strip_a })
}

/// Phase C: 1× KEF — Merge strips by removing the original edge.
pub(crate) fn phase_c_kef(
    topo: &mut TopoStore,
    contact: &ChamferContact,
    strip_a: FaceIdx,
) -> Result<FaceIdx, EulerChamferError> {
    kef(topo, contact.edge, strip_a)?;
    Ok(strip_a)
}

/// Phase D: 2× KEV — Remove corner vertices V_a and V_b.
pub(crate) fn phase_d_kev(
    topo: &mut TopoStore,
    semv_res: &SemvPhaseResult,
) -> Result<(), EulerChamferError> {
    kev(topo, semv_res.semv1.new_edge, semv_res.v_c1a)?;
    kev(topo, semv_res.semv4.new_edge, semv_res.v_c2b)?;
    Ok(())
}

/// Phase E: Update geometry — set chamfer face surface and edge curves.
fn phase_e_geom(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    contact: &ChamferContact,
    chamfer_face: FaceIdx,
) {
    // Set chamfer face surface to a new plane through the contact points.
    let chamfer_normal = (contact.c1_b - contact.c1_a).cross(&(contact.c2_a - contact.c1_a));
    let cn_len = chamfer_normal.norm();
    let chamfer_normal_unit = if cn_len > 1e-15 {
        chamfer_normal / cn_len
    } else {
        Vec3::zeros()
    };

    let n_a = plane_normal(geom, topo, contact.face_a);
    let n_b = plane_normal(geom, topo, contact.face_b);
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

    // Update curve_ids on the chamfer face edges.
    update_face_curves(topo, geom, chamfer_face);

    // Invalidate mesh_cache on the entire shell.
    let shell_idx = topo.faces.get(chamfer_face).shell;
    let all_faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for &f in &all_faces {
        topo.faces.get_mut(f).mesh_cache = None;
    }
}

/// Walk a face loop and set each edge's curve_id to a LineSegment.
pub(crate) fn update_face_curves(topo: &mut TopoStore, geom: &mut AnalyticalGeomStore, face: FaceIdx) {
    let loop_idx = topo.faces.get(face).outer_loop;
    let start = topo.loops.get(loop_idx).half_edge;
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
}

// ─── Single-edge Euler chamfer ───────────────────────────────────────────────

/// Apply a single-edge chamfer using the 9-op Euler sequence.
///
/// Returns the FaceIdx of the new chamfer face (quad).
pub(crate) fn euler_chamfer_single_edge(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    contact: &ChamferContact,
) -> Result<FaceIdx, EulerChamferError> {
    let _span = info_span!("euler_chamfer_single_edge", edge = contact.edge.raw()).entered();

    let semv_res = phase_a_semv(topo, geom, contact);
    let mef_res = phase_b_mef(topo, geom, contact, &semv_res)?;
    let chamfer_face = phase_c_kef(topo, contact, mef_res.strip_a)?;
    phase_d_kev(topo, &semv_res)?;
    phase_e_geom(topo, geom, contact, chamfer_face);

    Ok(chamfer_face)
}

// ─── Corner processing ──────────────────────────────────────────────────────

/// Result of corner topology processing: far edges + the corner face itself.
pub(crate) struct CornerProcessResult {
    pub(crate) far_edges: Vec<EdgeIdx>,
    pub(crate) corner_face: FaceIdx,
}

/// Info about one edge at a corner vertex, from the corner vertex's perspective.
pub(crate) struct CornerEdgeInfo {
    pub(crate) edge: EdgeIdx,
    /// Contact point position (at distance d from the corner vertex along the edge).
    pub(crate) contact_pt: Point3,
    /// The far endpoint (the other vertex of the edge, not the corner vertex).
    pub(crate) far_vertex: VertexIdx,
}

/// Process the topology of a corner vertex where N chamfered edges meet.
///
/// Performs corner SEMV → corner MEF → corner KEF+KEV → geometry assignment.
/// Returns the `far_edges` (sub-edges on the far side of each contact point)
/// so the caller can apply its own single-edge operation (chamfer or fillet).
pub(crate) fn process_corner_topology(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    corner: &CornerGroup,
    edges: &[EdgeIdx],
    geometries: &[ContactGeometry],
    edge_distances: &[f64],
) -> Result<CornerProcessResult, EulerChamferError> {
    let _span = info_span!("process_corner_topology", vertex = corner.vertex.raw(), n_edges = corner.edge_indices.len()).entered();
    let v_corner = corner.vertex;

    // Build per-edge info from the corner vertex's perspective
    let mut edge_infos: Vec<CornerEdgeInfo> = Vec::new();
    for (local_i, &idx) in corner.edge_indices.iter().enumerate() {
        let d_i = edge_distances[local_i];
        let cg = &geometries[idx];
        let edge_idx = edges[idx];
        let adj = edge_adjacency(topo, edge_idx)?;
        let va = topo.half_edges.get(adj.he_a).origin;
        let vb = topo.half_edges.get(adj.he_b).origin;

        // Determine which end is the corner vertex
        let (contact_pt, far_vertex) = if va == v_corner {
            let dir = (cg.pt_b - cg.pt_a).normalize();
            (cg.pt_a + dir * d_i, vb)
        } else {
            let dir = (cg.pt_a - cg.pt_b).normalize();
            (cg.pt_b + dir * d_i, va)
        };

        edge_infos.push(CornerEdgeInfo {
            edge: edge_idx,
            contact_pt,
            far_vertex,
        });
    }

    let n = edge_infos.len();

    // ── 4a. Corner SEMV: insert contact vertices near V on each edge ──
    // For each edge, find the HE from V_corner in that edge, and SEMV to split.
    let mut contact_vertices: Vec<VertexIdx> = Vec::with_capacity(n);
    let mut near_edges: Vec<EdgeIdx> = Vec::with_capacity(n);
    let mut far_edges: Vec<EdgeIdx> = Vec::with_capacity(n);

    for info in &edge_infos {
        let adj = edge_adjacency(topo, info.edge)?;
        let va = topo.half_edges.get(adj.he_a).origin;

        // Find the HE from V_corner on this edge
        let he_from_corner = if va == v_corner { adj.he_a } else { adj.he_b };

        // Register the contact point
        let pid = geom.add_point(info.contact_pt);
        let v_corner_pt = geom.points[topo.vertices.get(v_corner).point_id as usize];
        let far_pt = geom.points[topo.vertices.get(info.far_vertex).point_id as usize];

        let curve_near = geom.add_curve(CurveDef::LineSegment(LineSegment {
            start: v_corner_pt,
            end: info.contact_pt,
        }));
        let curve_far = geom.add_curve(CurveDef::LineSegment(LineSegment {
            start: info.contact_pt,
            end: far_pt,
        }));

        // SEMV: split the edge at he_from_corner. The original edge becomes V_corner↔P_i
        // (near), and new_edge becomes P_i↔V_far (far).
        let sr = semv(topo, he_from_corner, pid, curve_near, curve_far);

        contact_vertices.push(sr.new_vertex);
        // Original edge (he_from_corner's edge) = near edge (V_corner↔P_i)
        near_edges.push(topo.half_edges.get(he_from_corner).edge);
        // New edge = far edge (P_i↔V_far)
        far_edges.push(sr.new_edge);
    }

    // ── 4b. Corner MEF: on each face meeting at V, connect P_left and P_right ──
    // Identify which faces meet at V_corner, and which two contact vertices appear on each face.
    // Walk each face loop to find: ...P_left → V_corner → P_right...
    // Then MEF(he@P_right, he@P_left) to split off triangle [P_left, V, P_right].

    // Collect faces meeting at V_corner from the shell
    let shell_idx = {
        let he0 = topo.edges.get(near_edges[0]).half_edges[0];
        let l0 = topo.half_edges.get(he0).loop_ref;
        topo.faces.get(topo.loops.get(l0).face).shell
    };
    let all_faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

    // For each face, find if V_corner is in its loop and which 2 contact vertices bracket it.
    let contact_set: std::collections::HashSet<VertexIdx> = contact_vertices.iter().copied().collect();
    let mut strip_faces: Vec<FaceIdx> = Vec::new();
    let mut mef_edges: Vec<EdgeIdx> = Vec::new();

    for &face_idx in &all_faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start = topo.loops.get(loop_idx).half_edge;
        let mut he = start;
        let mut found_v = false;
        loop {
            if topo.half_edges.get(he).origin == v_corner {
                found_v = true;
                break;
            }
            he = topo.half_edges.get(he).next;
            if he == start {
                break;
            }
        }
        if !found_v {
            continue;
        }

        // he has origin=V_corner. Walk backward to find P_left (predecessor of V).
        let he_at_v = he;
        let he_prev = find_prev_he(topo, loop_idx, he_at_v);
        let p_left = topo.half_edges.get(he_prev).origin;
        // P_right is the destination of he_at_v: origin of he_at_v.next
        let he_next = topo.half_edges.get(he_at_v).next;
        let p_right = topo.half_edges.get(he_next).origin;

        // Both P_left and P_right should be contact vertices
        if !contact_set.contains(&p_left) || !contact_set.contains(&p_right) {
            continue;
        }

        // MEF(he@P_right, he@P_left) to split off triangle [P_left, V_corner, P_right].
        // After MEF(h1=he@P_right, h2=he@P_left):
        //   new loop: P_left → V_corner → P_right → he_b(P_right→P_left)
        //   = triangle strip

        // Find HE with origin=P_right and origin=P_left in this face
        let he_p_right = find_he_from_vertex(topo, face_idx, p_right)
            .ok_or(EulerChamferError::HalfEdgeNotFound {
                face: face_idx,
                vertex: p_right,
            })?;
        let he_p_left = find_he_from_vertex(topo, face_idx, p_left)
            .ok_or(EulerChamferError::HalfEdgeNotFound {
                face: face_idx,
                vertex: p_left,
            })?;

        let pt_left = geom.points[topo.vertices.get(p_left).point_id as usize];
        let pt_right = geom.points[topo.vertices.get(p_right).point_id as usize];
        let curve_id = geom.add_curve(CurveDef::LineSegment(LineSegment {
            start: pt_right,
            end: pt_left,
        }));
        let surface_id = topo.faces.get(face_idx).surface_id;

        let mef_r = mef(topo, he_p_right, he_p_left, surface_id, curve_id)?;
        strip_faces.push(mef_r.new_face);
        mef_edges.push(mef_r.new_edge);
    }

    // ── 4c. Corner KEF+KEV: merge strips, remove V_corner ──
    // Kill near edges via KEF, merging strips. After (N-1) KEFs, V_corner is in one face
    // with spike topology. Then kev_same_loop removes V_corner.

    if strip_faces.is_empty() {
        // Degenerate case: no faces to strip. Return far edges with a dummy face.
        let fallback_face = {
            let he0 = topo.edges.get(near_edges[0]).half_edges[0];
            let l0 = topo.half_edges.get(he0).loop_ref;
            topo.loops.get(l0).face
        };
        return Ok(CornerProcessResult {
            far_edges,
            corner_face: fallback_face,
        });
    }

    // KEF on each near edge. Each near edge connects V_corner to a contact vertex.
    // After MEFs, each near edge is shared between two strip faces (or a strip and another face).
    let mut surviving_face = strip_faces[0];

    for i in 0..near_edges.len() {
        let ne = near_edges[i];
        // Check if this edge is adjacent to surviving_face
        let [he0, he1] = topo.edges.get(ne).half_edges;
        let l0 = topo.half_edges.get(he0).loop_ref;
        let l1 = topo.half_edges.get(he1).loop_ref;
        let f0 = topo.loops.get(l0).face;
        let f1 = topo.loops.get(l1).face;

        if f0 == f1 {
            // Both HEs in the same face — this near edge is a spike already.
            // Skip KEF, will be handled by kev_same_loop below.
            continue;
        }

        let survive = if f0 == surviving_face {
            f0
        } else if f1 == surviving_face {
            f1
        } else {
            // surviving_face not adjacent; pick whichever is a strip face
            if strip_faces.contains(&f0) { f0 } else { f1 }
        };

        kef(topo, ne, survive)?;
        surviving_face = survive;
    }

    // Now V_corner should be connected by at most one near edge still in the surviving face.
    // Use kev_same_loop to remove V_corner.
    // Find the remaining edge connecting V_corner to a contact vertex.
    // Walk the surviving face's loop to find HE with origin=V_corner.
    let sf_loop = topo.faces.get(surviving_face).outer_loop;
    let sf_start = topo.loops.get(sf_loop).half_edge;
    let mut he = sf_start;
    let mut v_corner_edge = None;
    loop {
        if topo.half_edges.get(he).origin == v_corner {
            v_corner_edge = Some(topo.half_edges.get(he).edge);
            break;
        }
        he = topo.half_edges.get(he).next;
        if he == sf_start {
            break;
        }
    }

    if let Some(remaining_edge) = v_corner_edge {
        // Find which contact vertex to keep
        let [re0, re1] = topo.edges.get(remaining_edge).half_edges;
        let ro0 = topo.half_edges.get(re0).origin;
        let ro1 = topo.half_edges.get(re1).origin;
        let v_keep = if ro0 == v_corner { ro1 } else { ro0 };
        kev_same_loop(topo, remaining_edge, v_keep)?;
    }

    // ── 4d. Corner geometry: assign Plane surface to corner face ──
    // The corner face is the surviving face. Compute its plane from 3 contact points.
    if contact_vertices.len() >= 3 {
        let p0 = geom.points[topo.vertices.get(contact_vertices[0]).point_id as usize];
        let p1 = geom.points[topo.vertices.get(contact_vertices[1]).point_id as usize];
        let p2 = geom.points[topo.vertices.get(contact_vertices[2]).point_id as usize];

        let corner_normal = (p1 - p0).cross(&(p2 - p0));
        let cn_len = corner_normal.norm();
        if cn_len > 1e-15 {
            let corner_normal_unit = corner_normal / cn_len;
            // Orient outward: the corner vertex was "inside" relative to the corner face.
            // The normal should point away from the solid interior.
            let v_corner_pt = geom.points[topo.vertices.get(v_corner).point_id as usize];
            let to_corner = v_corner_pt - p0;
            let final_normal = if corner_normal_unit.dot(&to_corner) > 0.0 {
                -corner_normal_unit
            } else {
                corner_normal_unit
            };
            let sid = geom.add_surface(SurfaceDef::Plane(Plane {
                origin: p0,
                normal: final_normal,
            }));
            topo.faces.get_mut(surviving_face).surface_id = sid;
        }
    }

    // Update edge curves on the corner face
    update_face_curves(topo, geom, surviving_face);

    // Invalidate mesh caches
    let all_faces2: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for &f in &all_faces2 {
        topo.faces.get_mut(f).mesh_cache = None;
    }

    Ok(CornerProcessResult {
        far_edges,
        corner_face: surviving_face,
    })
}

// ─── Entity collection helper ────────────────────────────────────────────────

/// Collect all face, edge, and vertex indices from a solid by walking its shell.
pub(crate) fn collect_solid_entities(
    topo: &TopoStore,
    solid: SolidIdx,
) -> (HashSet<FaceIdx>, HashSet<EdgeIdx>, HashSet<VertexIdx>) {
    let shell = topo.solids.get(solid).outer_shell();
    let faces_list = &topo.shells.get(shell).faces;
    let mut faces = HashSet::new();
    let mut edges = HashSet::new();
    let mut vertices = HashSet::new();

    for &face_idx in faces_list {
        faces.insert(face_idx);
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start = topo.loops.get(loop_idx).half_edge;
        let mut he = start;
        loop {
            let he_data = topo.half_edges.get(he);
            vertices.insert(he_data.origin);
            edges.insert(he_data.edge);
            he = he_data.next;
            if he == start {
                break;
            }
        }
    }

    (faces, edges, vertices)
}

// ─── Public API ──────────────────────────────────────────────────────────────

/// Apply a constant-distance Euler-based chamfer to straight edges between planar faces.
///
/// Operates in-place on the solid (no rebuild). Supports independent edges and
/// edges meeting at corner vertices.
///
/// Returns the same SolidIdx (modified in-place) and a ShapeEvolution record.
pub fn euler_chamfer_edges(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    edges: &[EdgeIdx],
    distance: f64,
) -> Result<(SolidIdx, ShapeEvolution), EulerChamferError> {
    let _span = info_span!("euler_chamfer_edges", n_edges = edges.len(), distance).entered();

    // Snapshot entities before modification.
    let (faces_before, edges_before, verts_before) = collect_solid_entities(topo, solid);

    // Collect faces adjacent to chamfered edges (they'll be Modified).
    let mut adjacent_faces: HashSet<FaceIdx> = HashSet::new();
    for &edge_idx in edges {
        let adj = edge_adjacency(topo, edge_idx)?;
        adjacent_faces.insert(adj.face_a);
        adjacent_faces.insert(adj.face_b);
    }

    let (geometries, classification) = classify_edges(topo, geom, edges, distance)?;

    // Track new faces created by chamfer/corner operations.
    let mut chamfer_faces: Vec<(FaceIdx, EdgeIdx)> = Vec::new();
    let mut corner_faces: Vec<FaceIdx> = Vec::new();

    // 1. Process corners first (creates far sub-edges)
    for corner in &classification.corners {
        let edge_distances: Vec<f64> = vec![distance; corner.edge_indices.len()];
        let result = process_corner_topology(topo, geom, corner, edges, &geometries, &edge_distances)?;
        corner_faces.push(result.corner_face);
        // Far-end processing: chamfer each far sub-edge as an independent edge
        for (local_i, &far_edge) in result.far_edges.iter().enumerate() {
            let cg = compute_contact_geometry(topo, geom, far_edge, distance)?;
            let ct = query_contact_topology(topo, far_edge)?;
            let contact = assemble_contact(&cg, &ct);
            let chamfer_face = euler_chamfer_single_edge(topo, geom, &contact)?;
            let original_edge = edges[corner.edge_indices[local_i]];
            chamfer_faces.push((chamfer_face, original_edge));
        }
    }

    // 2. Process independent edges (recompute topology per edge)
    for &idx in &classification.independent {
        let topo_data = query_contact_topology(topo, edges[idx])?;
        let contact = assemble_contact(&geometries[idx], &topo_data);
        let chamfer_face = euler_chamfer_single_edge(topo, geom, &contact)?;
        chamfer_faces.push((chamfer_face, edges[idx]));
    }

    // Build evolution by diffing before/after.
    let (faces_after, edges_after, verts_after) = collect_solid_entities(topo, solid);
    let mut evo = ShapeEvolution::new();

    // Track explicitly-created chamfer faces.
    let corner_face_set: HashSet<FaceIdx> = corner_faces.iter().copied().collect();

    for &face in &faces_after {
        if let Some(&(_, edge)) = chamfer_faces.iter().find(|(f, _)| *f == face) {
            evo.record_face(face, FaceOrigin::FromEdge(edge));
        } else if corner_face_set.contains(&face) {
            // Corner faces are new faces generated at corner vertices.
            // Use FromEdge with the first edge of the corner group as provenance.
            evo.record_face(face, FaceOrigin::Primitive);
        } else if faces_before.contains(&face) {
            // Face existed before — was it adjacent to a chamfered edge?
            if adjacent_faces.contains(&face) {
                evo.record_face(face, FaceOrigin::Modified(face));
            } else {
                evo.record_face(face, FaceOrigin::CopiedFrom(face));
            }
        } else {
            // New face not explicitly tracked — created by corner topology.
            evo.record_face(face, FaceOrigin::Primitive);
        }
    }

    // Edge provenance.
    let chamfered_edge_set: HashSet<EdgeIdx> = edges.iter().copied().collect();
    // Edges adjacent to chamfered edges (share a face with a chamfered edge).
    let mut adjacent_edges: HashSet<EdgeIdx> = HashSet::new();
    for &face in &adjacent_faces {
        if faces_before.contains(&face) {
            let loop_idx = topo.faces.get(face).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                let e = topo.half_edges.get(he).edge;
                if !chamfered_edge_set.contains(&e) {
                    adjacent_edges.insert(e);
                }
                he = topo.half_edges.get(he).next;
                if he == start { break; }
            }
        }
    }

    for &edge in &edges_after {
        if edges_before.contains(&edge) {
            if adjacent_edges.contains(&edge) {
                evo.record_edge(edge, EdgeOrigin::Modified(edge));
            } else {
                evo.record_edge(edge, EdgeOrigin::CopiedFrom(edge));
            }
        } else {
            // New edge — created by chamfer/corner Euler ops.
            // Check if it borders a chamfer face (FromIntersection of adjacent faces)
            // or is entirely new (Primitive for corner edges).
            evo.record_edge(edge, EdgeOrigin::Primitive);
        }
    }

    // Vertex provenance.
    for &vert in &verts_after {
        if verts_before.contains(&vert) {
            evo.record_vertex(vert, VertexOrigin::CopiedFrom(vert));
        } else {
            // New vertex — chamfer contact point or corner vertex.
            evo.record_vertex(vert, VertexOrigin::Primitive);
        }
    }

    // Deleted faces/edges/vertices.
    for &face in &faces_before {
        if !faces_after.contains(&face) {
            evo.record_deleted_face(face);
        }
    }
    for &edge in &edges_before {
        if !edges_after.contains(&edge) {
            evo.record_deleted_edge(edge);
        }
    }
    for &vert in &verts_before {
        if !verts_after.contains(&vert) {
            evo.deleted_vertices.push(vert);
        }
    }

    Ok((solid, evo))
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box_builder::make_box_into;
    use crate::edge_analysis::{edge_adjacency, plane_normal, solid_edges};
    use rustkernel_topology::diagnostics::validate_solid;
    use std::collections::HashSet;

    /// Count V, E, F by walking the shell.
    fn count_vef(topo: &TopoStore, solid: SolidIdx) -> (usize, usize, usize) {
        let shell = topo.solids.get(solid).outer_shell();
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
        let shell = topo.solids.get(solid).outer_shell();
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
        let (_solid, _evo) = result.unwrap();

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
        let (_solid, _evo) = euler_chamfer_edges(&mut topo, &mut geom, solid, &[edge], dist).unwrap();

        // Collect all vertex positions
        let shell = topo.solids.get(solid).outer_shell();
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
        let (_solid, _evo) = result.unwrap();

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

        let (_solid, _evo) = euler_chamfer_edges(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();

        let report = validate_solid(&topo, solid);
        assert!(
            report.is_valid(),
            "validate_solid failed: {:?}",
            report.errors()
        );
    }

    #[test]
    fn test_euler_chamfer_three_non_adjacent_edges() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Pick 3 non-adjacent edges (no shared vertices):
        // bottom-front (-Z/-Y), top-back (+Z/+Y), left-right (+X/-X... actually pick 3 that share no vertices)
        // A box has 8 vertices. 3 non-adjacent edges means 6 distinct vertices.
        // bottom-front: between -Z and -Y
        // top-back: between +Z and +Y
        // left-middle: between -X and +Z — shares a vertex with top-back? Let's be careful.
        // Use: bottom-front (-Z/-Y), top-back (+Z/+Y), left-right... pick one that is fully disjoint.
        // bottom-front edge touches (0,0,0)-(2,0,0)
        // top-back edge touches (0,2,2)-(2,2,2)
        // The remaining vertices are (0,0,2),(2,0,2),(0,2,0),(2,2,0)
        // An edge between -X and +Y goes from (0,2,0) to (0,2,2) — shares with top-back at (0,2,2).
        // An edge between +X and -Y goes from (2,0,0) to (2,0,2) — shares with bottom-front at (2,0,0).
        // An edge between -X and -Z: (0,0,0) to (0,2,0) — shares with bottom-front at (0,0,0).
        // So on a 2x2x2 box, any 3 edges either share a vertex or are parallel pairs.
        // Let's use a larger box or pick 2 opposite edges + verify with 2 non-adjacent.
        // Actually pick: bottom-front, top-back — those share no vertices. That's only 2.
        // For 3 non-adjacent: we need to pick from the 4 "middle" edges (connecting top/bottom).
        // middle-left-front between -X and -Y: (0,0,0)-(0,0,2) — shares (0,0,0) with bottom-front.
        // Hmm, on a unit box every vertex touches 3 edges. Let's just pick 2 opposite + verify it works.
        // Actually for testing 3 independent edges: pick bottom-front, top-back, and middle-right-back.
        // bottom-front: (0,0,0)-(2,0,0) on -Z/-Y
        // top-back: (0,2,2)-(2,2,2) on +Z/+Y
        // middle-right-back: (2,2,0)-(2,2,2)? That shares (2,2,2) with top-back.
        // It seems like on a cube, any 3 edges will share at least one vertex unless they're the
        // 3 parallel edges of one direction. Let's check:
        // The 4 edges parallel to X: (0,0,0)-(2,0,0), (0,0,2)-(2,0,2), (0,2,0)-(2,2,0), (0,2,2)-(2,2,2)
        // Pick any 3 of these 4 — none share vertices? No: e.g. none of these share vertices among the group.
        // Each has 2 vertices from {(0,y,z),(2,y,z)} and no two share a y,z pair.
        // So pick 3 edges parallel to X: bottom-front(-Z,-Y), bottom-back(-Z,+Y), top-front(+Z,-Y).
        let e1 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, -1.0, 0.0),
        );
        let e2 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, 1.0, 0.0),
        );
        let e3 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, -1.0, 0.0),
        );

        let result = euler_chamfer_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3);
        assert!(result.is_ok(), "Three-edge chamfer failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        // After: 3 chamfers on a box: V = 8 + 3*2 = 14, E = 12 + 3*3 = 21, F = 6 + 3 = 9
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (14, 21, 9), "Three chamfers: V=14, E=21, F=9");
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");

        verify_all_twins(&topo, solid);

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "validate_solid failed: {:?}", report.errors());
    }

    #[test]
    fn test_euler_chamfer_three_edges_one_corner() {
        // Pick 3 edges that share one vertex of a box.
        // On a 2x2x2 box at origin, vertex (0,0,0) has 3 edges:
        //   bottom-front (-Z/-Y): connects (0,0,0)-(2,0,0)
        //   bottom-left (-Z/-X): connects (0,0,0)-(0,2,0)... wait, need to check normals.
        // Edges at vertex (0,0,0):
        //   1. between -Z and -Y faces (bottom-front edge along X)
        //   2. between -Z and -X faces (bottom-left edge along Y)
        //   3. between -Y and -X faces (front-left edge along Z)
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        assert_eq!(count_vef(&topo, solid), (8, 12, 6));

        let e1 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, -1.0, 0.0),
        );
        let e2 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(-1.0, 0.0, 0.0),
        );
        let e3 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, -1.0, 0.0), Vec3::new(-1.0, 0.0, 0.0),
        );

        // Verify these 3 edges share a vertex
        let adj1 = edge_adjacency(&topo, e1).unwrap();
        let adj2 = edge_adjacency(&topo, e2).unwrap();
        let adj3 = edge_adjacency(&topo, e3).unwrap();
        let verts1: HashSet<_> = [topo.half_edges.get(adj1.he_a).origin, topo.half_edges.get(adj1.he_b).origin].into();
        let verts2: HashSet<_> = [topo.half_edges.get(adj2.he_a).origin, topo.half_edges.get(adj2.he_b).origin].into();
        let verts3: HashSet<_> = [topo.half_edges.get(adj3.he_a).origin, topo.half_edges.get(adj3.he_b).origin].into();
        let shared: Vec<_> = verts1.intersection(&verts2).copied().filter(|v| verts3.contains(v)).collect();
        assert_eq!(shared.len(), 1, "3 edges should share exactly 1 vertex");

        let result = euler_chamfer_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3);
        assert!(result.is_ok(), "Corner chamfer failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        // After chamfering 3 edges at a corner:
        // The corner adds 1 triangle face + 3 chamfer faces (from far-end processing).
        // Corner SEMV: +3V, +3E (split each edge at contact point)
        // Corner MEF: +3E, +3F (split each face at corner)
        // Corner KEF×2 + KEV: -(2E+2F) -(1V+1E) = -3E, -2F, -1V
        // Corner net: +2V, +3E, +1F (corner face)
        // Far-end processing: 3 × single-edge chamfer = 3 × (+2V, +3E, +1F)
        // Total: 8 + 2 + 6 = 16V, 12 + 3 + 9 = 24E, 6 + 1 + 3 = 10F
        // Euler: 16 - 24 + 10 = 2 ✓
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");

        verify_all_twins(&topo, solid);

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "validate_solid failed: {:?}", report.errors());
    }

    #[test]
    fn test_euler_chamfer_mixed_corner_and_independent() {
        // 3 edges at a corner + 1 independent edge (on opposite side of box)
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // 3 corner edges at (0,0,0)
        let e1 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, -1.0, 0.0),
        );
        let e2 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(-1.0, 0.0, 0.0),
        );
        let e3 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, -1.0, 0.0), Vec3::new(-1.0, 0.0, 0.0),
        );

        // 1 independent edge: top-back-right (between +Z and +Y, no vertex shared with corner)
        let e4 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0),
        );

        let result = euler_chamfer_edges(&mut topo, &mut geom, solid, &[e1, e2, e3, e4], 0.3);
        assert!(result.is_ok(), "Mixed chamfer failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");

        verify_all_twins(&topo, solid);

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "validate_solid failed: {:?}", report.errors());
    }

    #[test]
    fn test_euler_chamfer_one_edge_evolution() {
        use rustkernel_topology::evolution::{EdgeOrigin, FaceOrigin, VertexOrigin};

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (_solid, evo) = euler_chamfer_edges(&mut topo, &mut geom, solid, &[edge], 0.3).unwrap();

        // After: 7 faces total.
        assert_eq!(evo.face_provenance.len(), 7, "Should track 7 faces");

        let from_edge = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::FromEdge(_)))
            .count();
        let modified = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::Modified(_)))
            .count();
        let copied = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::CopiedFrom(_)))
            .count();

        assert_eq!(from_edge, 1, "1 chamfer face → FromEdge");
        assert_eq!(modified, 2, "2 adjacent faces → Modified");
        assert_eq!(copied, 4, "4 untouched faces → CopiedFrom");

        // The chamfered edge should be in deleted_edges.
        assert!(evo.deleted_edges.contains(&edge), "Chamfered edge should be deleted");

        // Edge provenance: all edges in result should be tracked.
        assert!(!evo.edge_provenance.is_empty(), "Should have edge provenance");
        let copied_edges = evo.edge_provenance.values()
            .filter(|o| matches!(o, EdgeOrigin::CopiedFrom(_)))
            .count();
        let modified_edges = evo.edge_provenance.values()
            .filter(|o| matches!(o, EdgeOrigin::Modified(_)))
            .count();
        let new_edges = evo.edge_provenance.values()
            .filter(|o| matches!(o, EdgeOrigin::Primitive))
            .count();
        // Chamfer replaces 1 edge with 2 new edges (chamfer face borders).
        // Adjacent edges on the 2 modified faces get Modified provenance.
        assert!(copied_edges > 0, "Some edges should be CopiedFrom");
        assert!(modified_edges > 0, "Adjacent edges should be Modified");
        assert!(new_edges > 0, "Chamfer should create new edges");

        // Vertex provenance: all vertices in result should be tracked.
        assert!(!evo.vertex_provenance.is_empty(), "Should have vertex provenance");
        let copied_verts = evo.vertex_provenance.values()
            .filter(|o| matches!(o, VertexOrigin::CopiedFrom(_)))
            .count();
        let new_verts = evo.vertex_provenance.values()
            .filter(|o| matches!(o, VertexOrigin::Primitive))
            .count();
        assert!(copied_verts > 0, "Original vertices should be CopiedFrom");
        assert!(new_verts > 0, "Chamfer contact points should be new vertices");
    }

    #[test]
    fn test_euler_chamfer_corner_evolution() {
        use rustkernel_topology::evolution::FaceOrigin;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // Pick 3 edges that share vertex (0,0,0) — same pattern as existing corner test.
        let e1 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, -1.0, 0.0),
        );
        let e2 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0), Vec3::new(-1.0, 0.0, 0.0),
        );
        let e3 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, -1.0, 0.0), Vec3::new(-1.0, 0.0, 0.0),
        );

        let (_solid, evo) = euler_chamfer_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3).unwrap();

        // Should have face provenance entries and tracked deletions.
        assert!(!evo.face_provenance.is_empty(), "Should have face provenance");
        assert!(!evo.deleted_edges.is_empty(), "Should have deleted edges");

        // Should have chamfer faces (FromEdge) and at least one corner face (Primitive).
        let from_edge = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::FromEdge(_)))
            .count();
        let primitive = evo.face_provenance.values()
            .filter(|o| matches!(o, FaceOrigin::Primitive))
            .count();
        assert!(from_edge >= 3, "Should have at least 3 chamfer faces from edges, got {}", from_edge);
        assert!(primitive >= 1, "Should have at least 1 corner face (Primitive), got {}", primitive);
    }
}

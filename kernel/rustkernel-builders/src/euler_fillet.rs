//! Euler-operator-based fillet builder.
//!
//! Reuses phases A-D from `euler_chamfer` (identical topology sequence), then adds:
//! - Phase E: SEMV arc intermediate vertices on transverse edges.
//! - Phase F: Assign CylinderSurface + LineSegment chord curves.
//!
//! Supports filleting multiple edges including edges meeting at a corner vertex.

use std::collections::{HashMap, HashSet};

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::euler::{find_he_from_vertex, semv, EulerError};
use rustkernel_topology::evolution::{
    FaceOrigin, ShapeEvolution,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{
    AnalyticalGeomStore, CurveDef, CylinderSurface, LineSegment, SphereSurface, SurfaceDef,
};

use crate::edge_analysis::{
    dihedral_angle, edge_adjacency, edge_endpoints, plane_normal,
};
use crate::euler_chamfer::{
    assemble_contact, collect_solid_entities, compute_contact_geometry, compute_inward,
    phase_a_semv, phase_b_mef, phase_c_kef, phase_d_kev, process_corner_topology,
    query_contact_topology, update_face_curves, ContactGeometry, EulerChamferError,
};

// ─── Error type ──────────────────────────────────────────────────────────────

#[derive(Debug)]
pub enum EulerFilletError {
    EdgeAnalysis(crate::edge_analysis::EdgeAnalysisError),
    NotConvex(EdgeIdx),
    Euler(EulerError),
    Chamfer(EulerChamferError),
    HalfEdgeNotFound { face: FaceIdx, vertex: VertexIdx },
}

impl std::fmt::Display for EulerFilletError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EulerFilletError::EdgeAnalysis(e) => write!(f, "Edge analysis error: {e}"),
            EulerFilletError::NotConvex(e) => write!(f, "Edge {:?} is not convex", e),
            EulerFilletError::Euler(e) => write!(f, "Euler operator error: {e}"),
            EulerFilletError::Chamfer(e) => write!(f, "Chamfer phase error: {e}"),
            EulerFilletError::HalfEdgeNotFound { face, vertex } => {
                write!(
                    f,
                    "No half-edge from vertex {:?} in face {:?}",
                    vertex, face
                )
            }
        }
    }
}

impl std::error::Error for EulerFilletError {}

impl From<crate::edge_analysis::EdgeAnalysisError> for EulerFilletError {
    fn from(e: crate::edge_analysis::EdgeAnalysisError) -> Self {
        EulerFilletError::EdgeAnalysis(e)
    }
}

impl From<EulerError> for EulerFilletError {
    fn from(e: EulerError) -> Self {
        EulerFilletError::Euler(e)
    }
}

impl From<EulerChamferError> for EulerFilletError {
    fn from(e: EulerChamferError) -> Self {
        EulerFilletError::Chamfer(e)
    }
}

// ─── Fillet geometry ────────────────────────────────────────────────────────

/// Fillet-specific geometry for one edge. Includes the base contact geometry
/// (reused by chamfer phases A-D) plus arc and fillet surface data.
struct FilletGeometry {
    /// Base contact geometry: c1_a, c1_b, c2_a, c2_b positions.
    base: ContactGeometry,
    /// Intermediate arc points at vertex A (not including c1_a and c2_a).
    arc_a: Vec<Point3>,
    /// Intermediate arc points at vertex B (not including c1_b and c2_b).
    arc_b: Vec<Point3>,
    /// Fillet center at the vertex-A end (on the rolling ball axis).
    center_a: Point3,
    /// Edge direction (unit vector from A to B).
    edge_dir: Vec3,
    /// Fillet radius.
    radius: f64,
    /// Arc subtended angle (π - dihedral_angle).
    subtended_angle: f64,
    /// The fillet surface (Cylinder for Plane-Plane, Torus for Plane-Cylinder, etc.)
    fillet_surface: SurfaceDef,
}

/// Generate arc points by Rodrigues rotation of `(start - center)` around `axis`.
///
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
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let r_rot = r.scale(cos_a)
            + axis.cross(&r).scale(sin_a)
            + axis.scale(axis.dot(&r) * (1.0 - cos_a));
        result.push(center + r_rot);
    }
    result
}

/// Compute fillet geometry for one edge: contact positions, arc points, fillet surface.
///
/// Works for any surface pair — uses surface normals at the edge midpoint to
/// compute the rolling-ball center, and determines the fillet surface type
/// from the adjacent face surface types.
fn compute_fillet_geometry(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
    radius: f64,
    arc_segments: usize,
) -> Result<FilletGeometry, EulerFilletError> {
    use crate::edge_analysis::surface_normal_at_point;
    use crate::fillet_geometry::determine_fillet_surface;

    let alpha = dihedral_angle(topo, geom, edge_idx)?;
    let d = radius / (alpha / 2.0).tan();

    // Reuse chamfer contact geometry computation with d as the offset distance.
    let base = compute_contact_geometry(topo, geom, edge_idx, d)
        .map_err(EulerFilletError::Chamfer)?;

    let (pt_a, pt_b) = edge_endpoints(topo, geom, edge_idx);
    let edge_vec = pt_b - pt_a;
    let edge_len = edge_vec.norm();
    let edge_dir = if edge_len > 1e-15 {
        edge_vec / edge_len
    } else {
        Vec3::zeros()
    };
    let mid = Point3::from((pt_a.coords + pt_b.coords) * 0.5);

    let adj = edge_adjacency(topo, edge_idx)?;
    let inward_1 = compute_inward(topo, geom, adj.face_a, &mid, &edge_dir);

    // Surface normal at face_a — works for any surface type
    let n1 = surface_normal_at_point(geom, topo, adj.face_a, &mid);

    // Fillet center: vertex + d*inward_face_a - R*normal_face_a
    let center_a = pt_a + d * inward_1 - radius * n1;
    let center_b = pt_b + d * inward_1 - radius * n1;

    let subtended_angle = std::f64::consts::PI - alpha;

    let arc_a = generate_arc_points(&center_a, &base.c1_a, &edge_dir, subtended_angle, arc_segments);
    let arc_b = generate_arc_points(&center_b, &base.c1_b, &edge_dir, subtended_angle, arc_segments);

    // Determine fillet surface type from the adjacent face surfaces
    let surf_a = &geom.surfaces[topo.faces.get(adj.face_a).surface_id as usize];
    let surf_b = &geom.surfaces[topo.faces.get(adj.face_b).surface_id as usize];
    let fillet_surface = determine_fillet_surface(surf_a, surf_b, &center_a, &edge_dir, radius);

    Ok(FilletGeometry {
        base,
        arc_a,
        arc_b,
        center_a,
        edge_dir,
        radius,
        subtended_angle,
        fillet_surface,
    })
}

// ─── Phase E: Arc SEMV ──────────────────────────────────────────────────────

/// Phase E: Insert arc intermediate vertices on the two transverse edges of the
/// quad fillet face produced by phases A-D.
///
/// The quad has vertices c1_a, c1_b, c2_b, c2_a (in some loop order).
/// Transverse edges connect c1↔c2 at each endpoint (A and B).
/// We SEMV `arc_segments-1` times on each transverse edge.
fn fillet_phase_e_arc(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    fillet_face: FaceIdx,
    v_c1a: VertexIdx,
    v_c1b: VertexIdx,
    v_c2a: VertexIdx,
    v_c2b: VertexIdx,
    fillet_geom: &FilletGeometry,
) -> Result<(), EulerFilletError> {
    let arc_a = &fillet_geom.arc_a;
    let arc_b = &fillet_geom.arc_b;

    // ── Transverse edge at vertex A: c2_a → c1_a ──
    // The arc goes from c1_a → (arc_a points) → c2_a.
    // The edge direction in the loop is c2_a → c1_a.
    // We insert arc points from c2_a side toward c1_a, so in reverse order of arc_a.
    //
    // Find the HE from c2_a in the fillet face.
    if !arc_a.is_empty() {
        let mut he_from_c2a = find_he_from_vertex(topo, fillet_face, v_c2a)
            .ok_or(EulerFilletError::HalfEdgeNotFound {
                face: fillet_face,
                vertex: v_c2a,
            })?;

        // Walk to find the HE from c2_a that leads to c1_a (the transverse edge).
        let loop_idx = topo.faces.get(fillet_face).outer_loop;
        let start = topo.loops.get(loop_idx).half_edge;
        let mut he = start;
        let mut found = false;
        loop {
            let origin = topo.half_edges.get(he).origin;
            let next_origin = topo.half_edges.get(topo.half_edges.get(he).next).origin;
            if origin == v_c2a && next_origin == v_c1a {
                he_from_c2a = he;
                found = true;
                break;
            }
            he = topo.half_edges.get(he).next;
            if he == start {
                break;
            }
        }
        if !found {
            // The transverse edge might go c1_a → c2_a instead. Try that direction.
            he = start;
            loop {
                let origin = topo.half_edges.get(he).origin;
                let next_origin = topo.half_edges.get(topo.half_edges.get(he).next).origin;
                if origin == v_c1a && next_origin == v_c2a {
                    // Insert arc points in forward order on this edge.
                    he_from_c2a = he; // Actually this is he_from_c1a
                    found = true;
                    break;
                }
                he = topo.half_edges.get(he).next;
                if he == start {
                    break;
                }
            }
            if found {
                // Edge goes c1_a → c2_a. Insert arc_a points in forward order.
                let mut current_he = he_from_c2a;
                for pt in arc_a.iter() {
                    let pid = geom.add_point(*pt);
                    let origin_vid = topo.half_edges.get(current_he).origin;
                    let origin_pt = geom.points[topo.vertices.get(origin_vid).point_id as usize];
                    let next_he = topo.half_edges.get(current_he).next;
                    let dest_vid = topo.half_edges.get(next_he).origin;
                    let dest_pt = geom.points[topo.vertices.get(dest_vid).point_id as usize];
                    let curve_a = geom.add_curve(CurveDef::LineSegment(LineSegment {
                        start: origin_pt,
                        end: *pt,
                    }));
                    let curve_b = geom.add_curve(CurveDef::LineSegment(LineSegment {
                        start: *pt,
                        end: dest_pt,
                    }));
                    let sr = semv(topo, current_he, pid, curve_a, curve_b);
                    // Next SEMV operates on the new trailing half-edge.
                    current_he = sr.new_he_a;
                }
            }
        } else {
            // Edge goes c2_a → c1_a. Insert arc_a points in reverse order.
            let mut current_he = he_from_c2a;
            for pt in arc_a.iter().rev() {
                let pid = geom.add_point(*pt);
                let origin_vid = topo.half_edges.get(current_he).origin;
                let origin_pt = geom.points[topo.vertices.get(origin_vid).point_id as usize];
                let next_he = topo.half_edges.get(current_he).next;
                let dest_vid = topo.half_edges.get(next_he).origin;
                let dest_pt = geom.points[topo.vertices.get(dest_vid).point_id as usize];
                let curve_a = geom.add_curve(CurveDef::LineSegment(LineSegment {
                    start: origin_pt,
                    end: *pt,
                }));
                let curve_b = geom.add_curve(CurveDef::LineSegment(LineSegment {
                    start: *pt,
                    end: dest_pt,
                }));
                let sr = semv(topo, current_he, pid, curve_a, curve_b);
                current_he = sr.new_he_a;
            }
        }
    }

    // ── Transverse edge at vertex B: c1_b → c2_b ──
    // The arc goes from c1_b → (arc_b points) → c2_b.
    // Insert arc_b points in forward order.
    if !arc_b.is_empty() {
        let loop_idx = topo.faces.get(fillet_face).outer_loop;
        let start = topo.loops.get(loop_idx).half_edge;
        let mut he = start;
        let mut he_transverse_b = None;
        let mut direction_forward = true; // c1_b → c2_b

        // Find edge c1_b → c2_b
        loop {
            let origin = topo.half_edges.get(he).origin;
            let next_origin = topo.half_edges.get(topo.half_edges.get(he).next).origin;
            if origin == v_c1b && next_origin == v_c2b {
                he_transverse_b = Some(he);
                direction_forward = true;
                break;
            }
            if origin == v_c2b && next_origin == v_c1b {
                he_transverse_b = Some(he);
                direction_forward = false;
                break;
            }
            he = topo.half_edges.get(he).next;
            if he == start {
                break;
            }
        }

        if let Some(he_start) = he_transverse_b {
            let points: Vec<Point3> = if direction_forward {
                arc_b.iter().copied().collect()
            } else {
                arc_b.iter().rev().copied().collect()
            };

            let mut current_he = he_start;
            for pt in &points {
                let pid = geom.add_point(*pt);
                let origin_vid = topo.half_edges.get(current_he).origin;
                let origin_pt = geom.points[topo.vertices.get(origin_vid).point_id as usize];
                let next_he = topo.half_edges.get(current_he).next;
                let dest_vid = topo.half_edges.get(next_he).origin;
                let dest_pt = geom.points[topo.vertices.get(dest_vid).point_id as usize];
                let curve_a = geom.add_curve(CurveDef::LineSegment(LineSegment {
                    start: origin_pt,
                    end: *pt,
                }));
                let curve_b = geom.add_curve(CurveDef::LineSegment(LineSegment {
                    start: *pt,
                    end: dest_pt,
                }));
                let sr = semv(topo, current_he, pid, curve_a, curve_b);
                current_he = sr.new_he_a;
            }
        }
    }

    Ok(())
}

// ─── Phase F: CylinderSurface geometry ──────────────────────────────────────

/// Phase F: Assign fillet surface to the fillet face and update edge curves.
fn fillet_phase_f_geom(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    fillet_face: FaceIdx,
    fillet_geom: &FilletGeometry,
) {
    let surface_id = geom.add_surface(fillet_geom.fillet_surface.clone());
    topo.faces.get_mut(fillet_face).surface_id = surface_id;

    // Update all edge curves on the fillet face to LineSegment chords.
    update_face_curves(topo, geom, fillet_face);

    // Invalidate mesh caches on the entire shell.
    let shell_idx = topo.faces.get(fillet_face).shell;
    let all_faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for &f in &all_faces {
        topo.faces.get_mut(f).mesh_cache = None;
    }
}

// ─── Sphere corner eligibility ──────────────────────────────────────────────

/// Check if a corner is eligible for a SphereSurface patch.
///
/// Requirements: exactly 3 edges, all dihedral angles within `tol_deg` of 90°,
/// all adjacent faces are Plane.
fn is_sphere_eligible(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    corner: &crate::euler_chamfer::CornerGroup,
    edges: &[EdgeIdx],
    tol_deg: f64,
) -> bool {
    if corner.edge_indices.len() != 3 {
        return false;
    }
    let half_pi = std::f64::consts::FRAC_PI_2;
    let tol_rad = tol_deg.to_radians();
    for &idx in &corner.edge_indices {
        // Check dihedral angle is ~90°
        let alpha = match dihedral_angle(topo, geom, edges[idx]) {
            Ok(a) => a,
            Err(_) => return false,
        };
        if (alpha - half_pi).abs() > tol_rad {
            return false;
        }
        // Check both faces are Plane
        let adj = match edge_adjacency(topo, edges[idx]) {
            Ok(a) => a,
            Err(_) => return false,
        };
        if !matches!(geom.surfaces[topo.faces.get(adj.face_a).surface_id as usize], SurfaceDef::Plane(_)) {
            return false;
        }
        if !matches!(geom.surfaces[topo.faces.get(adj.face_b).surface_id as usize], SurfaceDef::Plane(_)) {
            return false;
        }
    }
    true
}

/// Collect unique outward face normals meeting at a corner vertex.
///
/// Must be called BEFORE topology surgery (which removes V_corner).
fn collect_corner_face_normals(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    corner: &crate::euler_chamfer::CornerGroup,
    edges: &[EdgeIdx],
) -> Vec<Vec3> {
    let mut normals = Vec::new();
    let mut seen_faces = std::collections::HashSet::new();
    for &idx in &corner.edge_indices {
        if let Ok(adj) = edge_adjacency(topo, edges[idx]) {
            for &face in &[adj.face_a, adj.face_b] {
                if seen_faces.insert(face) {
                    if matches!(geom.surfaces[topo.faces.get(face).surface_id as usize], SurfaceDef::Plane(_)) {
                        normals.push(plane_normal(geom, topo, face));
                    }
                }
            }
        }
    }
    normals
}

// ─── NURBS corner patch (Tier 3) ────────────────────────────────────────────

/// Compute a best-fit blend sphere center for any corner geometry.
/// For orthogonal corners this is exact; for non-orthogonal it's an approximation.
fn compute_blend_sphere(
    corner_pt: &Point3,
    face_normals: &[Vec3],
    radius: f64,
) -> (Point3, f64) {
    let normal_sum: Vec3 = face_normals.iter().copied().fold(Vec3::zeros(), |a, b| a + b);
    let center = corner_pt - radius * normal_sum;
    (center, radius)
}

/// Compute a weighted blend sphere for a corner with dissimilar radii.
///
/// Each face normal is offset by its corresponding radius (from the edge between
/// the other two faces). Returns a best-fit center and average radius.
fn compute_blend_sphere_weighted(
    corner_pt: &Point3,
    face_normals: &[Vec3],
    radii: &[f64],
) -> (Point3, f64) {
    // For 3 edges at a corner, there are 3 face normals. Each edge i sits between
    // two faces; the radius R_i offsets along the normal of the third face.
    // As an approximation, offset each normal by the average of the two adjacent edge radii.
    let avg_r = radii.iter().sum::<f64>() / radii.len() as f64;
    let mut offset = Vec3::zeros();
    for (i, n) in face_normals.iter().enumerate() {
        // Use the radius of the edge opposite to this face in the corner.
        // For a 3-edge corner: face i is opposite edge i (approximately).
        let r = if i < radii.len() { radii[i] } else { avg_r };
        offset += n.scale(r);
    }
    let center = corner_pt - offset;
    (center, avg_r)
}

/// Build a NURBS corner patch for dissimilar radii using a weighted blend sphere.
///
/// Similar to `build_corner_nurbs_patch` but uses per-edge radius weighting:
/// boundary points near a given edge are projected onto a sphere whose radius
/// is biased toward that edge's fillet radius, producing a smooth transition
/// between the three different-radius cylinders.
fn build_dissimilar_corner_patch(
    boundary_pts: &[Point3],
    sphere_center: &Point3,
    avg_radius: f64,
    corner_radii: &[f64],
    face_normals: &[Vec3],
) -> Result<curvo::prelude::NurbsSurface3D<f64>, EulerFilletError> {
    use curvo::prelude::{Interpolation, NurbsCurve3D, NurbsSurface3D};

    let make_err = || EulerFilletError::Chamfer(EulerChamferError::HalfEdgeNotFound {
        face: FaceIdx::from_raw(0),
        vertex: VertexIdx::from_raw(0),
    });

    if boundary_pts.len() < 3 {
        return Err(make_err());
    }

    let n = boundary_pts.len();

    // Compute centroid and project onto weighted sphere for apex.
    let centroid = {
        let sum: Vec3 = boundary_pts.iter().map(|p| p.coords).fold(Vec3::zeros(), |a, b| a + b);
        Point3::from(sum / n as f64)
    };
    let to_centroid = centroid - sphere_center;
    let to_centroid_len = to_centroid.norm();
    let apex = if to_centroid_len > 1e-15 {
        sphere_center + to_centroid.scale(avg_radius / to_centroid_len)
    } else {
        centroid
    };

    // For each boundary point, compute a local blend radius based on proximity
    // to each edge's region. We use the face normals to determine which edge
    // a boundary point is closest to, then weight the radius accordingly.
    let local_radius = |pt: &Point3| -> f64 {
        if corner_radii.is_empty() || face_normals.is_empty() {
            return avg_radius;
        }
        // Weight by inverse distance to each face plane (higher weight = closer to that face).
        let mut total_weight = 0.0_f64;
        let mut weighted_r = 0.0_f64;
        for (i, n) in face_normals.iter().enumerate() {
            let r = if i < corner_radii.len() { corner_radii[i] } else { avg_radius };
            // Distance from point to the face plane through sphere_center
            let dist = (pt - sphere_center).dot(n).abs().max(1e-10);
            let w = 1.0 / dist;
            total_weight += w;
            weighted_r += w * r;
        }
        if total_weight > 1e-15 {
            weighted_r / total_weight
        } else {
            avg_radius
        }
    };

    // Build rings at various fractions toward the apex, projected onto locally-radiused sphere.
    let make_ring = |t: f64| -> Result<NurbsCurve3D<f64>, EulerFilletError> {
        let ring_pts: Vec<Point3> = boundary_pts.iter().map(|p| {
            let interp = Point3::from(p.coords * (1.0 - t) + apex.coords * t);
            let r = local_radius(&interp);
            let to_interp = interp - sphere_center;
            let len = to_interp.norm();
            if len > 1e-15 {
                sphere_center + to_interp.scale(r / len)
            } else {
                interp
            }
        }).collect();
        let mut closed = ring_pts;
        closed.push(closed[0]);
        let degree = if closed.len() >= 4 { 3 } else { 1 };
        NurbsCurve3D::<f64>::interpolate(&closed, degree).map_err(|_| make_err())
    };

    // 4 rings for better quality with dissimilar radii.
    let curves = vec![
        make_ring(0.95)?,
        make_ring(0.65)?,
        make_ring(0.35)?,
        make_ring(0.0)?,
    ];
    let surface = NurbsSurface3D::<f64>::try_loft(&curves, Some(2))
        .map_err(|_| make_err())?;

    Ok(surface)
}

/// Project a point onto the closest position on a CylinderSurface.
fn project_onto_cylinder(pt: &Point3, cyl: &CylinderSurface) -> Point3 {
    let to_pt = pt - cyl.origin;
    let along = cyl.axis.scale(to_pt.dot(&cyl.axis));
    let radial = to_pt - along;
    let radial_len = radial.norm();
    if radial_len < 1e-15 {
        // Point is on the axis; return a point at radius along an arbitrary perpendicular.
        return *pt;
    }
    let radial_dir = radial / radial_len;
    cyl.origin + along + radial_dir.scale(cyl.radius.abs())
}

/// Subdivide each boundary edge of the corner face by inserting intermediate
/// vertices projected onto the adjacent surface (cylinder or plane).
///
/// Each edge gets `n_subdiv` intermediate points via SEMV.
fn subdivide_corner_boundary(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    corner_face: FaceIdx,
    n_subdiv: usize,
) -> Result<(), EulerFilletError> {
    use rustkernel_topology::euler::semv;

    if n_subdiv == 0 {
        return Ok(());
    }

    // Collect edges of the corner face and their twin's surface info.
    let loop_idx = topo.faces.get(corner_face).outer_loop;
    let start = topo.loops.get(loop_idx).half_edge;
    let mut he_list = Vec::new();
    let mut he = start;
    loop {
        he_list.push(he);
        he = topo.half_edges.get(he).next;
        if he == start {
            break;
        }
    }

    // Process each edge: collect (he, start_pt, end_pt, adjacent_surface_id)
    let mut edge_data: Vec<(HalfEdgeIdx, Point3, Point3, Option<CylinderSurface>)> = Vec::new();
    for &he_idx in &he_list {
        let origin_vid = topo.half_edges.get(he_idx).origin;
        let next_he = topo.half_edges.get(he_idx).next;
        let dest_vid = topo.half_edges.get(next_he).origin;
        let origin_pt = geom.points[topo.vertices.get(origin_vid).point_id as usize];
        let dest_pt = geom.points[topo.vertices.get(dest_vid).point_id as usize];

        // Find the adjacent face via twin
        let twin_opt = topo.half_edges.get(he_idx).twin;
        let adj_cyl = if let Some(twin) = twin_opt {
            let twin_loop = topo.half_edges.get(twin).loop_ref;
            let twin_face = topo.loops.get(twin_loop).face;
            let sid = topo.faces.get(twin_face).surface_id;
            if let SurfaceDef::Cylinder(ref c) = geom.surfaces[sid as usize] {
                Some(c.clone())
            } else {
                None
            }
        } else {
            None
        };

        edge_data.push((he_idx, origin_pt, dest_pt, adj_cyl));
    }

    // Now SEMV intermediate points on each edge (process in order).
    for (he_idx, origin_pt, dest_pt, adj_cyl) in &edge_data {
        let mut current_he = *he_idx;
        for i in 1..=n_subdiv {
            let t = i as f64 / (n_subdiv + 1) as f64;
            let mut interp = Point3::from(origin_pt.coords * (1.0 - t) + dest_pt.coords * t);

            // Project onto adjacent cylinder if available
            if let Some(cyl) = adj_cyl {
                interp = project_onto_cylinder(&interp, cyl);
            }

            let pid = geom.add_point(interp);
            let cur_origin_vid = topo.half_edges.get(current_he).origin;
            let cur_origin_pt = geom.points[topo.vertices.get(cur_origin_vid).point_id as usize];
            let cur_next = topo.half_edges.get(current_he).next;
            let cur_dest_vid = topo.half_edges.get(cur_next).origin;
            let cur_dest_pt = geom.points[topo.vertices.get(cur_dest_vid).point_id as usize];

            let curve_a = geom.add_curve(CurveDef::LineSegment(LineSegment {
                start: cur_origin_pt,
                end: interp,
            }));
            let curve_b = geom.add_curve(CurveDef::LineSegment(LineSegment {
                start: interp,
                end: cur_dest_pt,
            }));

            let sr = semv(topo, current_he, pid, curve_a, curve_b);
            current_he = sr.new_he_a;
        }
    }

    Ok(())
}

/// Collect boundary vertex positions of a face in loop order.
fn collect_face_boundary_points(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    face: FaceIdx,
) -> Vec<Point3> {
    let loop_idx = topo.faces.get(face).outer_loop;
    let start = topo.loops.get(loop_idx).half_edge;
    let mut pts = Vec::new();
    let mut he = start;
    loop {
        let vid = topo.half_edges.get(he).origin;
        pts.push(geom.points[topo.vertices.get(vid).point_id as usize]);
        he = topo.half_edges.get(he).next;
        if he == start {
            break;
        }
    }
    pts
}

/// Build a NURBS surface patch for a corner face using concentric ring loft.
///
/// Strategy: create concentric rings at 5%, 50%, and 100% of boundary size,
/// all projected onto the blend sphere. Loft between them.
fn build_corner_nurbs_patch(
    boundary_pts: &[Point3],
    sphere_center: &Point3,
    sphere_radius: f64,
) -> Result<curvo::prelude::NurbsSurface3D<f64>, EulerFilletError> {
    use curvo::prelude::{Interpolation, NurbsCurve3D, NurbsSurface3D};

    let make_err = || EulerFilletError::Chamfer(EulerChamferError::HalfEdgeNotFound {
        face: FaceIdx::from_raw(0),
        vertex: VertexIdx::from_raw(0),
    });

    if boundary_pts.len() < 3 {
        return Err(make_err());
    }

    // Compute centroid of boundary → project onto sphere for apex direction.
    let centroid = {
        let sum: Vec3 = boundary_pts.iter().map(|p| p.coords).fold(Vec3::zeros(), |a, b| a + b);
        Point3::from(sum / boundary_pts.len() as f64)
    };
    let to_centroid = centroid - sphere_center;
    let to_centroid_len = to_centroid.norm();
    let apex = if to_centroid_len > 1e-15 {
        sphere_center + to_centroid.scale(sphere_radius / to_centroid_len)
    } else {
        centroid
    };

    // Build a closed ring at a given fraction t toward the apex, projected onto sphere.
    let make_ring = |t: f64| -> Result<NurbsCurve3D<f64>, EulerFilletError> {
        let ring_pts: Vec<Point3> = boundary_pts.iter().map(|p| {
            let interp = Point3::from(p.coords * (1.0 - t) + apex.coords * t);
            let to_interp = interp - sphere_center;
            let len = to_interp.norm();
            if len > 1e-15 {
                sphere_center + to_interp.scale(sphere_radius / len)
            } else {
                interp
            }
        }).collect();
        let mut closed = ring_pts;
        closed.push(closed[0]);
        let degree = if closed.len() >= 4 { 3 } else { 1 };
        NurbsCurve3D::<f64>::interpolate(&closed, degree).map_err(|_| make_err())
    };

    // 3 concentric rings: inner (near apex), middle, outer (boundary).
    let inner_curve = make_ring(0.9)?;
    let mid_curve = make_ring(0.5)?;
    let outer_curve = make_ring(0.0)?;

    // Loft: inner → mid → outer
    let curves = vec![inner_curve, mid_curve, outer_curve];
    let surface = NurbsSurface3D::<f64>::try_loft(&curves, Some(2))
        .map_err(|_| make_err())?;

    Ok(surface)
}

/// Assign a NURBS surface patch to a corner face, checking normal orientation.
fn assign_nurbs_corner_surface(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    corner_face: FaceIdx,
    nurbs_patch: curvo::prelude::NurbsSurface3D<f64>,
    v_corner_pt: &Point3,
) {
    use rustkernel_geom::FlippableNurbs;

    // Get corner face centroid to check normal orientation.
    let boundary_pts = collect_face_boundary_points(topo, geom, corner_face);
    let centroid = {
        let sum: Vec3 = boundary_pts.iter().map(|p| p.coords).fold(Vec3::zeros(), |a, b| a + b);
        Point3::from(sum / boundary_pts.len() as f64)
    };

    // Sample NURBS normal at center parameter via finite differences.
    let surface_normal = {
        let pt = nurbs_patch.point_at(0.5, 0.5);
        let du = nurbs_patch.point_at(0.5 + 1e-6, 0.5) - pt;
        let dv = nurbs_patch.point_at(0.5, 0.5 + 1e-6) - pt;
        du.cross(&dv)
    };

    // Normal should point away from the original corner vertex (outward).
    let to_corner = v_corner_pt - centroid;
    let flipped = surface_normal.dot(&to_corner) > 0.0;

    let sid = geom.add_surface(SurfaceDef::Nurbs(FlippableNurbs {
        surface: nurbs_patch,
        flipped,
    }));
    topo.faces.get_mut(corner_face).surface_id = sid;
    update_face_curves(topo, geom, corner_face);

    // Invalidate mesh caches.
    let shell_idx = topo.faces.get(corner_face).shell;
    let all_faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for &f in &all_faces {
        topo.faces.get_mut(f).mesh_cache = None;
    }
}

// ─── Single-edge Euler fillet ───────────────────────────────────────────────

/// Apply a single-edge fillet using phases A-D (from euler_chamfer) + E + F.
///
/// Returns the FaceIdx of the new fillet face.
fn euler_fillet_single_edge(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    fillet_geom: &FilletGeometry,
    arc_segments: usize,
) -> Result<FaceIdx, EulerFilletError> {
    let edge_idx = fillet_geom.base.edge;
    let _span = info_span!("euler_fillet_single_edge", edge = edge_idx.raw()).entered();

    // Query fresh topology and assemble contact for phases A-D.
    let ct = query_contact_topology(topo, edge_idx)?;
    let contact = assemble_contact(&fillet_geom.base, &ct);

    // Phase A: 4× SEMV
    let semv_res = phase_a_semv(topo, geom, &contact);

    // Phase B: 2× MEF
    let mef_res = phase_b_mef(topo, geom, &contact, &semv_res)?;

    // Phase C: 1× KEF
    let fillet_face = phase_c_kef(topo, &contact, mef_res.strip_a)?;

    // Phase D: 2× KEV
    phase_d_kev(topo, &semv_res)?;

    // Phase E: Insert arc intermediate vertices on transverse edges
    if arc_segments > 2 {
        fillet_phase_e_arc(
            topo, geom, fillet_face,
            semv_res.v_c1a, semv_res.v_c1b, semv_res.v_c2a, semv_res.v_c2b,
            fillet_geom,
        )?;
    }

    // Phase F: Assign CylinderSurface + update curves
    fillet_phase_f_geom(topo, geom, fillet_face, fillet_geom);

    Ok(fillet_face)
}

// ─── Edge classification ────────────────────────────────────────────────────

/// Classify fillet edges and compute per-edge geometry with per-edge radii.
fn classify_fillet_edges_with_radii(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edges: &[EdgeIdx],
    radii: &[f64],
    arc_segments: usize,
) -> Result<(Vec<FilletGeometry>, Vec<ContactGeometry>, crate::euler_chamfer::EdgeClassification), EulerFilletError> {
    use crate::euler_chamfer::{CornerGroup, EdgeClassification};

    let mut fillet_geoms = Vec::with_capacity(edges.len());
    let mut contact_geoms = Vec::with_capacity(edges.len());
    let mut vert_to_edges: HashMap<VertexIdx, Vec<usize>> = HashMap::new();

    for (i, &edge_idx) in edges.iter().enumerate() {
        let fg = compute_fillet_geometry(topo, geom, edge_idx, radii[i], arc_segments)?;
        contact_geoms.push(fg.base.clone());
        fillet_geoms.push(fg);

        let adj = edge_adjacency(topo, edge_idx)?;
        let va = topo.half_edges.get(adj.he_a).origin;
        let vb = topo.half_edges.get(adj.he_b).origin;
        vert_to_edges.entry(va).or_default().push(i);
        vert_to_edges.entry(vb).or_default().push(i);
    }

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

    Ok((fillet_geoms, contact_geoms, EdgeClassification { independent, corners }))
}

// ─── Public API ──────────────────────────────────────────────────────────────

/// Apply a constant-radius Euler-based fillet to straight edges between planar faces.
///
/// Operates in-place on the solid (no rebuild). Supports independent edges and
/// edges meeting at corner vertices.
///
/// Returns the same SolidIdx (modified in-place).
pub fn euler_fillet_edges(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    edges: &[EdgeIdx],
    radius: f64,
    arc_segments: usize,
) -> Result<(SolidIdx, ShapeEvolution), EulerFilletError> {
    let radii: Vec<f64> = vec![radius; edges.len()];
    euler_fillet_edges_with_radii(topo, geom, solid, edges, &radii, arc_segments)
}

/// Apply per-edge-radius Euler-based fillets to straight edges between planar faces.
///
/// `radii[i]` is the fillet radius for `edges[i]`. When edges meet at a corner vertex,
/// dissimilar radii produce a NURBS blend patch; equal radii use the faster sphere or
/// concentric-ring NURBS patch.
///
/// Operates in-place on the solid (no rebuild).
pub fn euler_fillet_edges_with_radii(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    edges: &[EdgeIdx],
    radii: &[f64],
    arc_segments: usize,
) -> Result<(SolidIdx, ShapeEvolution), EulerFilletError> {
    let _span = info_span!("euler_fillet_edges_with_radii", n_edges = edges.len(), arc_segments).entered();
    assert_eq!(edges.len(), radii.len(), "edges and radii must have the same length");
    let arc_segments = arc_segments.max(2);

    // Snapshot entities before modification.
    let (faces_before, edges_before, verts_before) = collect_solid_entities(topo, solid);

    // Collect faces adjacent to filleted edges.
    let mut adjacent_faces: HashSet<FaceIdx> = HashSet::new();
    for &edge_idx in edges {
        let adj = edge_adjacency(topo, edge_idx)?;
        adjacent_faces.insert(adj.face_a);
        adjacent_faces.insert(adj.face_b);
    }

    let (fillet_geoms, contact_geoms, classification) =
        classify_fillet_edges_with_radii(topo, geom, edges, radii, arc_segments)?;

    // Track new faces created by fillet/corner operations.
    let mut fillet_faces: Vec<(FaceIdx, EdgeIdx)> = Vec::new();
    let mut corner_faces: Vec<FaceIdx> = Vec::new();

    // 1. Process corners first (topology surgery, returns far sub-edges).
    for corner in &classification.corners {
        // Per-edge radii and distances at this corner.
        let corner_radii: Vec<f64> = corner.edge_indices.iter()
            .map(|&idx| radii[idx])
            .collect();
        let all_radii_equal = corner_radii.windows(2)
            .all(|w| (w[0] - w[1]).abs() < 1e-12);

        // Capture face normals BEFORE topology surgery (V_corner removed after).
        let sphere_eligible = all_radii_equal
            && is_sphere_eligible(topo, geom, corner, edges, 5.0);
        let corner_normals = collect_corner_face_normals(topo, geom, corner, edges);
        let v_corner_pt = geom.points[topo.vertices.get(corner.vertex).point_id as usize];

        // Per-edge distance: d_i = R_i / tan(α_i / 2)
        let edge_distances: Vec<f64> = corner.edge_indices.iter()
            .enumerate()
            .map(|(local_i, &idx)| {
                let alpha = dihedral_angle(topo, geom, edges[idx])
                    .unwrap_or(std::f64::consts::FRAC_PI_2);
                corner_radii[local_i] / (alpha / 2.0).tan()
            })
            .collect();

        let result = process_corner_topology(
            topo, geom, corner, edges, &contact_geoms, &edge_distances,
        )?;

        corner_faces.push(result.corner_face);

        // Far-end processing: fillet each far sub-edge with its own radius.
        for (local_i, &far_edge) in result.far_edges.iter().enumerate() {
            let r_i = corner_radii[local_i];
            let fg = compute_fillet_geometry(topo, geom, far_edge, r_i, arc_segments)?;
            let fillet_face = euler_fillet_single_edge(topo, geom, &fg, arc_segments)?;
            let original_edge = edges[corner.edge_indices[local_i]];
            fillet_faces.push((fillet_face, original_edge));
        }

        // Assign curved surface to the corner face.
        if sphere_eligible && corner_normals.len() == 3 {
            // Tier 2: SphereSurface for orthogonal 3-edge equal-radius corners.
            let r = corner_radii[0]; // all equal
            let normal_sum: Vec3 = corner_normals.iter().copied().fold(Vec3::zeros(), |a, b| a + b);
            let sphere_center = v_corner_pt - r * normal_sum;
            let sid = geom.add_surface(SurfaceDef::Sphere(SphereSurface {
                center: sphere_center,
                radius: r,
            }));
            topo.faces.get_mut(result.corner_face).surface_id = sid;
            update_face_curves(topo, geom, result.corner_face);

            let shell_idx = topo.faces.get(result.corner_face).shell;
            let all_faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
            for &f in &all_faces {
                topo.faces.get_mut(f).mesh_cache = None;
            }
        } else if corner.edge_indices.len() >= 3 && corner_normals.len() >= 3 {
            if all_radii_equal {
                // Tier 3: NURBS corner patch for non-orthogonal equal-radius corners.
                let r = corner_radii[0];
                let (sphere_center, sphere_r) =
                    compute_blend_sphere(&v_corner_pt, &corner_normals, r);

                subdivide_corner_boundary(topo, geom, result.corner_face, 2)?;

                let boundary_pts = collect_face_boundary_points(topo, geom, result.corner_face);
                match build_corner_nurbs_patch(&boundary_pts, &sphere_center, sphere_r) {
                    Ok(nurbs_patch) => {
                        assign_nurbs_corner_surface(
                            topo, geom, result.corner_face, nurbs_patch, &v_corner_pt,
                        );
                    }
                    Err(_) => {
                        tracing::warn!(
                            n_boundary = boundary_pts.len(),
                            "NURBS corner patch construction failed; keeping Plane fallback"
                        );
                    }
                }
            } else {
                // Tier 4: Dissimilar-radii NURBS corner patch.
                let (sphere_center, sphere_r) =
                    compute_blend_sphere_weighted(&v_corner_pt, &corner_normals, &corner_radii);

                subdivide_corner_boundary(topo, geom, result.corner_face, 3)?;

                let boundary_pts = collect_face_boundary_points(topo, geom, result.corner_face);
                match build_dissimilar_corner_patch(&boundary_pts, &sphere_center, sphere_r, &corner_radii, &corner_normals) {
                    Ok(nurbs_patch) => {
                        assign_nurbs_corner_surface(
                            topo, geom, result.corner_face, nurbs_patch, &v_corner_pt,
                        );
                    }
                    Err(_) => {
                        tracing::warn!(
                            corner_radii = ?corner_radii,
                            "Dissimilar corner patch failed; falling back to averaged-radius"
                        );
                        // Fallback: use averaged radius with the Tier 3 approach.
                        let avg_r = corner_radii.iter().sum::<f64>() / corner_radii.len() as f64;
                        let (sc, sr) = compute_blend_sphere(&v_corner_pt, &corner_normals, avg_r);
                        let boundary_pts = collect_face_boundary_points(topo, geom, result.corner_face);
                        if let Ok(patch) = build_corner_nurbs_patch(&boundary_pts, &sc, sr) {
                            assign_nurbs_corner_surface(
                                topo, geom, result.corner_face, patch, &v_corner_pt,
                            );
                        }
                    }
                }
            }
        }
    }

    // 2. Process independent edges.
    for &idx in &classification.independent {
        let fillet_face = euler_fillet_single_edge(topo, geom, &fillet_geoms[idx], arc_segments)?;
        fillet_faces.push((fillet_face, edges[idx]));
    }

    // Build evolution by diffing before/after.
    let (faces_after, edges_after, verts_after) = collect_solid_entities(topo, solid);
    let mut evo = ShapeEvolution::new();

    let corner_face_set: HashSet<FaceIdx> = corner_faces.iter().copied().collect();

    for &face in &faces_after {
        if let Some(&(_, edge)) = fillet_faces.iter().find(|(f, _)| *f == face) {
            evo.record_face(face, FaceOrigin::FromEdge(edge));
        } else if corner_face_set.contains(&face) {
            evo.record_face(face, FaceOrigin::Primitive);
        } else if faces_before.contains(&face) {
            if adjacent_faces.contains(&face) {
                evo.record_face(face, FaceOrigin::Modified(face));
            } else {
                evo.record_face(face, FaceOrigin::CopiedFrom(face));
            }
        } else {
            evo.record_face(face, FaceOrigin::Primitive);
        }
    }

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
    use crate::edge_analysis::{edge_adjacency, edge_convexity, plane_normal, solid_edges, EdgeConvexity};
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

    // ── Step 2 verification: geometry computation ──

    #[test]
    fn test_fillet_geometry_90deg() {
        // For a box edge (90° dihedral), d = R/tan(45°) = R.
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let radius = 0.3;
        let fg = compute_fillet_geometry(&topo, &geom, edge, radius, 4).unwrap();

        // d should equal R for 90° dihedral
        let (pt_a, pt_b) = edge_endpoints(&topo, &geom, edge);
        let d = radius; // R / tan(45°) = R

        // Contact points should be at distance d from vertices along inward directions
        let dist_c1a = (fg.base.c1_a - pt_a).norm();
        let dist_c2a = (fg.base.c2_a - pt_a).norm();
        assert!((dist_c1a - d).abs() < 1e-10, "c1_a distance: {dist_c1a} != {d}");
        assert!((dist_c2a - d).abs() < 1e-10, "c2_a distance: {dist_c2a} != {d}");

        // Arc subtended angle should be 90° (π - π/2)
        assert!((fg.subtended_angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10);

        // Arc points should have 3 intermediate points for 4 segments
        assert_eq!(fg.arc_a.len(), 3, "4 arc segments → 3 intermediate points");
        assert_eq!(fg.arc_b.len(), 3);

        // All arc points should be at distance R from cylinder center
        for pt in &fg.arc_a {
            let to_center = pt - fg.center_a;
            let along_axis = fg.edge_dir.scale(to_center.dot(&fg.edge_dir));
            let radial = to_center - along_axis;
            assert!(
                (radial.norm() - radius).abs() < 1e-10,
                "Arc point not on cylinder: dist = {}",
                radial.norm()
            );
        }
    }

    // ── Step 3 verification: single-edge fillet ──

    #[test]
    fn test_euler_fillet_one_edge() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        assert_eq!(count_vef(&topo, solid), (8, 12, 6));

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[edge], 0.3, 4);
        assert!(result.is_ok(), "Euler fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        // After fillet with 4 arc segments:
        // V = 10 + 2*(4-1) = 16, E = 15 + 2*(4-1) = 21, F = 7
        // Euler: 16 - 21 + 7 = 2
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(f, 7, "Fillet 1 edge → 7 faces");
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");

        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_euler_fillet_one_edge_min_segments() {
        // With arc_segments=2 (minimum), no Phase E SEMV needed.
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[edge], 0.3, 2);
        assert!(result.is_ok(), "Euler fillet (2 seg) failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        // With 2 segments, no intermediate arc points. Same as chamfer topology: V=10, E=15, F=7
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (10, 15, 7), "Fillet 2-seg: V=10, E=15, F=7");
        assert_eq!(v as i32 - e as i32 + f as i32, 2);
        verify_all_twins(&topo, solid);
    }

    // ── Step 4 verification: multi-edge ──

    #[test]
    fn test_euler_fillet_two_opposite_edges() {
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

        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2], 0.3, 4);
        assert!(result.is_ok(), "Two-edge fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(f, 8, "Fillet 2 edges → 8 faces");
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_euler_fillet_three_non_adjacent() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // 3 parallel edges (no shared vertices)
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

        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3, 4);
        assert!(result.is_ok(), "Three-edge fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(f, 9, "Fillet 3 edges → 9 faces");
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
        verify_all_twins(&topo, solid);
    }

    // ── Step 5 verification: corner + validate ──

    #[test]
    fn test_euler_fillet_three_edges_one_corner() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // 3 edges at vertex (0,0,0)
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

        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3, 4);
        assert!(result.is_ok(), "Corner fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_euler_fillet_mixed_corner_and_independent() {
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

        // 1 independent edge (opposite side)
        let e4 = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0),
        );

        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2, e3, e4], 0.3, 4);
        assert!(result.is_ok(), "Mixed fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_euler_fillet_arc_on_cylinder() {
        // Verify that arc vertices actually lie on the CylinderSurface.
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let radius = 0.3;
        let (_solid, _evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &[edge], radius, 8).unwrap();

        // Find the fillet face (the one with CylinderSurface)
        let shell = topo.solids.get(solid).outer_shell();
        let faces = &topo.shells.get(shell).faces.clone();
        let mut fillet_face = None;
        for &f in faces {
            let sid = topo.faces.get(f).surface_id;
            if let SurfaceDef::Cylinder(ref cyl) = geom.surfaces[sid as usize] {
                if (cyl.radius - radius).abs() < 1e-10 {
                    fillet_face = Some((f, cyl.origin, cyl.axis, cyl.radius));
                    break;
                }
            }
        }

        let (face, origin, axis, r) = fillet_face.expect("Should have a cylinder fillet face");

        // Walk the fillet face and check that all vertices are at distance R from the cylinder axis.
        let loop_idx = topo.faces.get(face).outer_loop;
        let start = topo.loops.get(loop_idx).half_edge;
        let mut he = start;
        let mut n_verts = 0;
        loop {
            let vid = topo.half_edges.get(he).origin;
            let pt = geom.points[topo.vertices.get(vid).point_id as usize];
            let to_pt = pt - origin;
            let along = axis.scale(to_pt.dot(&axis));
            let radial_dist = (to_pt - along).norm();
            assert!(
                (radial_dist - r).abs() < 1e-8,
                "Vertex {} at radial distance {radial_dist} != {r}",
                vid.raw()
            );
            n_verts += 1;
            he = topo.half_edges.get(he).next;
            if he == start {
                break;
            }
        }

        // With 8 arc segments: 2 longitudinal + 2*(8-1+1) = 2 + 16 = 18 vertices?
        // Actually: 2 + 2*8 = 18 vertices on the fillet face boundary.
        assert!(n_verts > 4, "Fillet face should have more than 4 vertices with 8 arc segments, got {n_verts}");
    }

    #[test]
    fn test_euler_fillet_validate_solid() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (_solid, _evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();

        let report = validate_solid(&topo, solid);
        assert!(
            report.is_valid(),
            "validate_solid failed: {:?}",
            report.errors()
        );
    }

    #[test]
    fn test_per_edge_distance_correctness() {
        // Verify that process_corner_topology uses per-edge distances by computing
        // them ourselves and checking that R/tan(α/2) == R for a 90° box edge.
        use crate::euler_chamfer::{
            compute_contact_geometry, process_corner_topology, CornerProcessResult,
        };

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // 3 edges at vertex (0,0,0). All 90° → d = R/tan(45°) = R.
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
        let edges_vec = vec![e1, e2, e3];
        let radius = 0.3;

        // Compute per-edge distances
        let distances: Vec<f64> = edges_vec.iter().map(|&e| {
            let alpha = dihedral_angle(&topo, &geom, e).unwrap();
            radius / (alpha / 2.0).tan()
        }).collect();

        // All should be R for 90° dihedral
        for &d in &distances {
            assert!((d - radius).abs() < 1e-10, "d should equal radius for 90° dihedral, got {d}");
        }

        // Also verify the full fillet works with these distances
        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &edges_vec, radius, 4);
        assert!(result.is_ok(), "Per-edge distance fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
        verify_all_twins(&topo, solid);
    }

    // ── Tier 2: Sphere corner tests ──

    /// Helper: find the sphere face in a solid after fillet.
    fn find_sphere_face(
        topo: &TopoStore,
        geom: &AnalyticalGeomStore,
        solid: SolidIdx,
    ) -> Option<(FaceIdx, Point3, f64)> {
        let shell = topo.solids.get(solid).outer_shell();
        for &f in &topo.shells.get(shell).faces {
            let sid = topo.faces.get(f).surface_id;
            if let SurfaceDef::Sphere(ref s) = geom.surfaces[sid as usize] {
                return Some((f, s.center, s.radius));
            }
        }
        None
    }

    #[test]
    fn test_sphere_corner_surface_type() {
        // Fillet 3 edges at a box corner → corner face should get SphereSurface.
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // 3 edges meeting at vertex (0,0,0)
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

        let radius = 0.3;
        let (_solid, _evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], radius, 4).unwrap();

        // Should have a sphere face
        let (face, center, r) = find_sphere_face(&topo, &geom, solid)
            .expect("Corner face should be SphereSurface");

        // Box centered at origin, size 2×2×2 → corner vertex at (-1,-1,-1).
        // Normals: (-1,0,0), (0,-1,0), (0,0,-1). Sum = (-1,-1,-1).
        // center = (-1,-1,-1) - R*(-1,-1,-1) = (-1+R, -1+R, -1+R)
        assert!((r - radius).abs() < 1e-10, "Sphere radius should be {radius}, got {r}");
        let expected_center = Point3::new(-1.0 + radius, -1.0 + radius, -1.0 + radius);
        let center_err = (center - expected_center).norm();
        assert!(
            center_err < 1e-10,
            "Sphere center should be ({},{},{}), got ({},{},{}), err={center_err}",
            expected_center.x, expected_center.y, expected_center.z,
            center.x, center.y, center.z
        );
    }

    #[test]
    fn test_sphere_corner_euler_valid() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

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

        let (_solid, _evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3, 4).unwrap();

        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
        verify_all_twins(&topo, solid);

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "validate_solid failed: {:?}", report.errors());
    }

    #[test]
    fn test_sphere_corner_tessellation() {
        use rustkernel_topology::tessellate;
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

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

        let (_solid, _evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &[e1, e2, e3], 0.3, 4).unwrap();

        // Tessellate the sphere corner face — should produce nonzero triangles.
        let (face, _, _) = find_sphere_face(&topo, &geom, solid)
            .expect("Should have sphere face");

        tessellate::tessellate_face(&mut topo, face, &geom);
        let mesh = topo.faces.get(face).mesh_cache.as_ref()
            .expect("Tessellation should produce a mesh");
        assert!(
            !mesh.indices.is_empty(),
            "Sphere corner face should tessellate to nonzero triangles"
        );

        // Normals should not all be identical (smooth shading from sphere).
        if mesh.normals.len() >= 2 {
            let first = &mesh.normals[0];
            let all_same = mesh.normals.iter().all(|n| (n - first).norm() < 1e-10);
            assert!(!all_same, "Sphere face normals should vary (smooth shading)");
        }
    }

    #[test]
    fn test_non_orthogonal_no_sphere() {
        // A non-orthogonal corner should NOT get SphereSurface (falls back to Plane).
        // Use an extruded right-triangle prism; vertex at (2,0,0) has a 45° dihedral
        // (hypotenuse face vs bottom face), so sphere eligibility should fail.
        use crate::extrude_builder::make_extrude_into;
        use crate::edge_analysis::solid_edges;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        // Right triangle in XY plane: (0,0,0), (2,0,0), (0,2,0)
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let solid = make_extrude_into(
            &mut topo, &mut geom, &profile, Vec3::new(0.0, 0.0, 1.0), 2.0,
        );

        // Find 3 edges meeting at vertex (2,0,0) — the bottom-right vertex.
        // One dihedral here is 45° (between front face y=0 and hypotenuse face).
        let target = Point3::new(2.0, 0.0, 0.0);
        let all_edges = solid_edges(&topo, solid);
        let corner_edges: Vec<EdgeIdx> = all_edges.iter().filter(|&&e| {
            let (pa, pb) = edge_endpoints(&topo, &geom, e);
            (pa - target).norm() < 1e-10 || (pb - target).norm() < 1e-10
        }).copied().collect();

        assert!(corner_edges.len() >= 3, "Expected 3 edges at (2,0,0), got {}", corner_edges.len());

        // Only fillet convex edges
        let convex_edges: Vec<EdgeIdx> = corner_edges.iter().filter(|&&e| {
            matches!(edge_convexity(&topo, &geom, e), Ok(EdgeConvexity::Convex))
        }).copied().collect();

        if convex_edges.len() >= 3 {
            let result = euler_fillet_edges(
                &mut topo, &mut geom, solid, &convex_edges[..3], 0.2, 4,
            );
            if let Ok((_solid, _evo)) = result {
                // No SphereSurface should exist (non-orthogonal corner)
                assert!(
                    find_sphere_face(&topo, &geom, solid).is_none(),
                    "Non-orthogonal corner should NOT get SphereSurface"
                );
            }
        }
    }

    // ── Tier 3: NURBS corner tests ──

    /// Helper: find NURBS faces in a solid.
    fn find_nurbs_face(
        topo: &TopoStore,
        geom: &AnalyticalGeomStore,
        solid: SolidIdx,
    ) -> Option<FaceIdx> {
        let shell = topo.solids.get(solid).outer_shell();
        for &f in &topo.shells.get(shell).faces {
            let sid = topo.faces.get(f).surface_id;
            if matches!(geom.surfaces[sid as usize], SurfaceDef::Nurbs(_)) {
                return Some(f);
            }
        }
        None
    }

    /// Build a triangular prism with a 45° dihedral at one vertex.
    fn make_test_prism() -> (TopoStore, AnalyticalGeomStore, SolidIdx) {
        use crate::extrude_builder::make_extrude_into;
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        // Right triangle: (0,0,0), (2,0,0), (0,2,0)
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let solid = make_extrude_into(
            &mut topo, &mut geom, &profile, Vec3::new(0.0, 0.0, 1.0), 2.0,
        );
        (topo, geom, solid)
    }

    /// Get 3 convex edges at a vertex of the prism.
    fn get_prism_corner_edges(
        topo: &TopoStore,
        geom: &AnalyticalGeomStore,
        solid: SolidIdx,
        target: &Point3,
    ) -> Vec<EdgeIdx> {
        use crate::edge_analysis::solid_edges;
        let all_edges = solid_edges(topo, solid);
        let corner_edges: Vec<EdgeIdx> = all_edges.iter().filter(|&&e| {
            let (pa, pb) = edge_endpoints(topo, geom, e);
            (pa - target).norm() < 1e-10 || (pb - target).norm() < 1e-10
        }).copied().collect();

        corner_edges.into_iter().filter(|&e| {
            matches!(edge_convexity(topo, geom, e), Ok(EdgeConvexity::Convex))
        }).collect()
    }

    #[test]
    fn test_nurbs_corner_patch_exists() {
        // Fillet 3 edges at the (2,0,0) vertex of a prism → non-orthogonal → NURBS patch.
        let (mut topo, mut geom, solid) = make_test_prism();
        let target = Point3::new(2.0, 0.0, 0.0);
        let convex_edges = get_prism_corner_edges(&topo, &geom, solid, &target);
        assert!(convex_edges.len() >= 3, "Expected ≥3 convex edges at (2,0,0)");

        let (_solid, _evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &convex_edges[..3], 0.15, 4).unwrap();

        assert!(
            find_nurbs_face(&topo, &geom, solid).is_some(),
            "Non-orthogonal corner should get NURBS surface patch"
        );
    }

    #[test]
    fn test_nurbs_corner_euler_valid() {
        let (mut topo, mut geom, solid) = make_test_prism();
        let target = Point3::new(2.0, 0.0, 0.0);
        let convex_edges = get_prism_corner_edges(&topo, &geom, solid, &target);

        if convex_edges.len() >= 3 {
            let result = euler_fillet_edges(
                &mut topo, &mut geom, solid, &convex_edges[..3], 0.15, 4,
            );
            if let Ok((_solid, _evo)) = result {
                let (v, e, f) = count_vef(&topo, solid);
                assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2");
                verify_all_twins(&topo, solid);

                let report = validate_solid(&topo, solid);
                assert!(report.is_valid(), "validate_solid failed: {:?}", report.errors());
            }
        }
    }

    #[test]
    fn test_nurbs_corner_tessellation() {
        use rustkernel_topology::tessellate;

        let (mut topo, mut geom, solid) = make_test_prism();
        let target = Point3::new(2.0, 0.0, 0.0);
        let convex_edges = get_prism_corner_edges(&topo, &geom, solid, &target);

        if convex_edges.len() >= 3 {
            let result = euler_fillet_edges(
                &mut topo, &mut geom, solid, &convex_edges[..3], 0.15, 4,
            );
            if let Ok((_solid, _evo)) = result {
                if let Some(face) = find_nurbs_face(&topo, &geom, solid) {
                    tessellate::tessellate_face(&mut topo, face, &geom);
                    let mesh = topo.faces.get(face).mesh_cache.as_ref()
                        .expect("NURBS corner face should tessellate");
                    assert!(
                        !mesh.indices.is_empty(),
                        "NURBS corner face should produce nonzero triangles"
                    );
                }
            }
        }
    }

    #[test]
    fn test_nurbs_corner_boundary_on_cylinders() {
        // After NURBS corner subdivision, intermediate boundary vertices near cylinder
        // faces should be within tolerance of the cylinder surface.
        let (mut topo, mut geom, solid) = make_test_prism();
        let target = Point3::new(2.0, 0.0, 0.0);
        let convex_edges = get_prism_corner_edges(&topo, &geom, solid, &target);

        if convex_edges.len() >= 3 {
            let result = euler_fillet_edges(
                &mut topo, &mut geom, solid, &convex_edges[..3], 0.15, 4,
            );
            if let Ok((_solid, _evo)) = result {
                if let Some(nurbs_face) = find_nurbs_face(&topo, &geom, solid) {
                    // Walk boundary of NURBS face
                    let loop_idx = topo.faces.get(nurbs_face).outer_loop;
                    let start = topo.loops.get(loop_idx).half_edge;
                    let mut he = start;
                    let mut n_verts = 0;
                    loop {
                        n_verts += 1;
                        he = topo.half_edges.get(he).next;
                        if he == start {
                            break;
                        }
                    }
                    // With 2 subdivisions per edge and 3 original edges,
                    // boundary should have more than 3 vertices.
                    assert!(
                        n_verts > 3,
                        "NURBS corner face should have > 3 boundary vertices after subdivision, got {n_verts}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_dissimilar_radii_box_corner() {
        // Fillet 3 edges at the (0,0,0) corner of a box with different radii.
        // This triggers Tier 4 (dissimilar-radius NURBS corner patch).
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        // 3 edges meeting at vertex (0,0,0): bottom-front, bottom-left, front-left
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

        let result = euler_fillet_edges_with_radii(
            &mut topo, &mut geom, solid,
            &[e1, e2, e3],
            &[0.2, 0.3, 0.4],
            4,
        );
        assert!(result.is_ok(), "Dissimilar-radii fillet failed: {:?}", result.err());
        let (_solid, _evo) = result.unwrap();

        // Euler formula: V - E + F = 2
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!(v as i32 - e as i32 + f as i32, 2, "Euler: V-E+F=2 (got V={v} E={e} F={f})");

        // All half-edges must have matched twins.
        verify_all_twins(&topo, solid);

        // The corner face must have a NURBS surface (not sphere, not plane).
        let nurbs_face = find_nurbs_face(&topo, &geom, solid);
        assert!(nurbs_face.is_some(), "Dissimilar corner should produce a NURBS surface");

        // Full diagnostic validation.
        let report = validate_solid(&topo, solid);
        assert!(
            report.errors().is_empty(),
            "Validation errors: {:?}",
            report.errors()
        );
    }

    #[test]
    fn test_euler_fillet_one_edge_evolution() {
        use rustkernel_topology::evolution::FaceOrigin;

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);

        let edge = find_box_edge_between(
            &topo, &geom, solid,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, -1.0, 0.0),
        );

        let (_solid, evo) = euler_fillet_edges(&mut topo, &mut geom, solid, &[edge], 0.3, 4).unwrap();

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

        assert_eq!(from_edge, 1, "1 fillet face → FromEdge");
        assert_eq!(modified, 2, "2 adjacent faces → Modified");
        assert_eq!(copied, 4, "4 untouched faces → CopiedFrom");

        assert!(evo.deleted_edges.contains(&edge), "Filleted edge should be deleted");
    }

    // ── Curved-edge fillet tests (cylinder cap-to-side) ──

    #[test]
    fn test_euler_fillet_cylinder_cap_edge() {
        use crate::cylinder_builder::make_cylinder_into;
        use crate::edge_analysis::{solid_edges, edge_adjacency, surface_normal_at_point, edge_midpoint};

        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cylinder_into(
            &mut topo, &mut geom,
            Point3::origin(), 3.0, 5.0, 16,
        );

        // Find a cap-to-side edge (Plane adjacent to Cylinder)
        let edges = solid_edges(&topo, solid);
        let mut cap_side_edge = None;
        for &e in &edges {
            let adj = edge_adjacency(&topo, e).unwrap();
            let sid_a = topo.faces.get(adj.face_a).surface_id;
            let sid_b = topo.faces.get(adj.face_b).surface_id;
            let is_plane_a = matches!(&geom.surfaces[sid_a as usize], SurfaceDef::Plane(_));
            let is_cyl_b = matches!(&geom.surfaces[sid_b as usize], SurfaceDef::Cylinder(_));
            let is_cyl_a = matches!(&geom.surfaces[sid_a as usize], SurfaceDef::Cylinder(_));
            let is_plane_b = matches!(&geom.surfaces[sid_b as usize], SurfaceDef::Plane(_));
            if (is_plane_a && is_cyl_b) || (is_cyl_a && is_plane_b) {
                cap_side_edge = Some(e);
                break;
            }
        }
        let edge = cap_side_edge.expect("Should find a cap-to-side edge on a cylinder");

        // Fillet this single edge
        let (v0, e0, f0) = count_vef(&topo, solid);
        let result = euler_fillet_edges(&mut topo, &mut geom, solid, &[edge], 0.3, 8);
        assert!(result.is_ok(), "Curved-edge fillet should succeed: {:?}", result.err());

        let (new_solid, _evo) = result.unwrap();
        let (v1, e1, f1) = count_vef(&topo, new_solid);

        // Should have gained faces (fillet face + arc subdivisions)
        assert!(f1 > f0, "Should gain faces: {f0} → {f1}");

        // Euler formula: V - E + F = 2 (genus 0)
        let euler = v1 as i32 - e1 as i32 + f1 as i32;
        assert_eq!(euler, 2, "Euler formula V-E+F={euler}, V={v1} E={e1} F={f1}");

        // Verify all twins are matched
        verify_all_twins(&topo, new_solid);

        // Check that the fillet face has a Torus surface
        let shell = topo.solids.get(new_solid).outer_shell();
        let faces = &topo.shells.get(shell).faces;
        let has_torus = faces.iter().any(|&f| {
            let sid = topo.faces.get(f).surface_id;
            matches!(&geom.surfaces[sid as usize], SurfaceDef::Torus(_))
        });
        assert!(has_torus, "Filleted cylinder should have a torus face");
    }
}

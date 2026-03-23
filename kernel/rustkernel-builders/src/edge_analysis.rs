use std::collections::HashSet;

use rustkernel_geom::{AnalyticalGeomStore, SurfaceDef};
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::geom_store::GeomAccess;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

/// Information about the two faces adjacent to an edge.
pub struct EdgeAdjacency {
    pub edge: EdgeIdx,
    pub face_a: FaceIdx,
    pub face_b: FaceIdx,
    pub he_a: HalfEdgeIdx,
    pub he_b: HalfEdgeIdx,
}

/// Classification of an edge's convexity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeConvexity {
    Convex,
    Concave,
    Flat,
}

/// Errors that can occur during edge analysis.
#[derive(Debug)]
pub enum EdgeAnalysisError {
    UnmatchedTwin(EdgeIdx),
    NonPlanarFace(FaceIdx),
    CurvedEdge(EdgeIdx),
}

impl std::fmt::Display for EdgeAnalysisError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EdgeAnalysisError::UnmatchedTwin(e) => write!(f, "Edge {:?} has no matched twin", e),
            EdgeAnalysisError::NonPlanarFace(face) => {
                write!(f, "Face {:?} is not planar", face)
            }
            EdgeAnalysisError::CurvedEdge(e) => write!(f, "Edge {:?} is not straight", e),
        }
    }
}

impl std::error::Error for EdgeAnalysisError {}

/// Find the face that contains a given half-edge by walking loop_ref → face.
pub fn face_of_half_edge(topo: &TopoStore, he: HalfEdgeIdx) -> FaceIdx {
    let loop_idx = topo.half_edges.get(he).loop_ref;
    topo.loops.get(loop_idx).face
}

/// Get the two adjacent faces and half-edges for an edge.
pub fn edge_adjacency(
    topo: &TopoStore,
    edge_idx: EdgeIdx,
) -> Result<EdgeAdjacency, EdgeAnalysisError> {
    let edge = topo.edges.get(edge_idx);
    let he_a = edge.half_edges[0];

    // Use twin link rather than edge.half_edges[1], because the twin-matching
    // loop processes both (a,b) and (b,a) directions and may overwrite half_edges[1].
    let he_b = topo
        .half_edges
        .get(he_a)
        .twin
        .ok_or(EdgeAnalysisError::UnmatchedTwin(edge_idx))?;

    let face_a = face_of_half_edge(topo, he_a);
    let face_b = face_of_half_edge(topo, he_b);

    Ok(EdgeAdjacency {
        edge: edge_idx,
        face_a,
        face_b,
        he_a,
        he_b,
    })
}

/// Get the outward normal of a planar face.
pub fn plane_normal(geom: &AnalyticalGeomStore, topo: &TopoStore, face_idx: FaceIdx) -> rustkernel_math::Vec3 {
    let surface_id = topo.faces.get(face_idx).surface_id;
    match &geom.surfaces[surface_id as usize] {
        SurfaceDef::Plane(p) => p.normal,
        _ => panic!("Expected planar face"),
    }
}

/// Get the origin point of a planar face's surface.
pub fn plane_origin(geom: &AnalyticalGeomStore, topo: &TopoStore, face_idx: FaceIdx) -> rustkernel_math::Point3 {
    let surface_id = topo.faces.get(face_idx).surface_id;
    match &geom.surfaces[surface_id as usize] {
        SurfaceDef::Plane(p) => p.origin,
        _ => panic!("Expected planar face"),
    }
}

/// Compute the centroid of a face's outer loop vertices.
pub fn face_centroid(topo: &TopoStore, geom: &AnalyticalGeomStore, face_idx: FaceIdx) -> rustkernel_math::Point3 {
    let loop_idx = topo.faces.get(face_idx).outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;
    let mut sum = rustkernel_math::Vec3::zeros();
    let mut count = 0usize;
    let mut he = start_he;
    loop {
        let vert = topo.half_edges.get(he).origin;
        let pt = geom.points[topo.vertices.get(vert).point_id as usize];
        sum += pt.coords;
        count += 1;
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }
    rustkernel_math::Point3::from(sum / count as f64)
}

/// Get the midpoint of the edge by averaging its two endpoint positions.
pub fn edge_midpoint(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
) -> rustkernel_math::Point3 {
    let edge = topo.edges.get(edge_idx);
    let he_a = edge.half_edges[0];
    let he_b = topo.half_edges.get(he_a).twin.unwrap();
    let v0 = topo.half_edges.get(he_a).origin;
    let v1 = topo.half_edges.get(he_b).origin;
    let p0 = geom.points[topo.vertices.get(v0).point_id as usize];
    let p1 = geom.points[topo.vertices.get(v1).point_id as usize];
    rustkernel_math::Point3::from((p0.coords + p1.coords) * 0.5)
}

/// Get the two endpoint positions of an edge.
pub fn edge_endpoints(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
) -> (rustkernel_math::Point3, rustkernel_math::Point3) {
    let edge = topo.edges.get(edge_idx);
    let he_a = edge.half_edges[0];
    let he_b = topo.half_edges.get(he_a).twin.unwrap();
    let v0 = topo.half_edges.get(he_a).origin;
    let v1 = topo.half_edges.get(he_b).origin;
    let p0 = geom.points[topo.vertices.get(v0).point_id as usize];
    let p1 = geom.points[topo.vertices.get(v1).point_id as usize];
    (p0, p1)
}

/// Evaluate the surface normal of a face at a given 3D point.
///
/// Uses inverse UV mapping to find the parametric coordinates, then evaluates
/// the surface normal. Works for any surface type (Plane, Cylinder, Sphere, etc.).
pub fn surface_normal_at_point(
    geom: &AnalyticalGeomStore,
    topo: &TopoStore,
    face_idx: FaceIdx,
    point: &Point3,
) -> Vec3 {
    let surface_id = topo.faces.get(face_idx).surface_id;
    let (u, v) = geom.surface_inverse_uv(surface_id, point);
    geom.surface_normal(surface_id, u, v)
}

/// Classify an edge as convex, concave, or flat.
///
/// Works for any surface type — evaluates surface normals at the edge midpoint.
///
/// Uses the centroid-based inward direction: the direction from the edge midpoint
/// toward face_a's centroid represents the "into face_a" direction. If this direction
/// is opposite to face_b's outward normal (dot < 0), the edge is convex.
pub fn edge_convexity(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
) -> Result<EdgeConvexity, EdgeAnalysisError> {
    let adj = edge_adjacency(topo, edge_idx)?;

    let mid = edge_midpoint(topo, geom, edge_idx);
    let centroid_a = face_centroid(topo, geom, adj.face_a);
    let normal_b = surface_normal_at_point(geom, topo, adj.face_b, &mid);

    // inward_a: direction from edge midpoint toward face_a's interior.
    let inward_a = centroid_a - mid;
    let dot = inward_a.dot(&normal_b);

    const TOL: f64 = 1e-10;
    if dot < -TOL {
        // face_a's interior is on the negative-normal side of face_b → convex
        Ok(EdgeConvexity::Convex)
    } else if dot > TOL {
        Ok(EdgeConvexity::Concave)
    } else {
        Ok(EdgeConvexity::Flat)
    }
}

/// Compute the interior dihedral angle between two faces sharing an edge.
///
/// Works for any surface type — evaluates surface normals at the edge midpoint.
/// Returns the angle in radians. For a box, all edges have dihedral = PI/2.
pub fn dihedral_angle(
    topo: &TopoStore,
    geom: &AnalyticalGeomStore,
    edge_idx: EdgeIdx,
) -> Result<f64, EdgeAnalysisError> {
    let adj = edge_adjacency(topo, edge_idx)?;

    let mid = edge_midpoint(topo, geom, edge_idx);
    let n1 = surface_normal_at_point(geom, topo, adj.face_a, &mid);
    let n2 = surface_normal_at_point(geom, topo, adj.face_b, &mid);

    // The outward normals of two faces sharing a convex edge point away from each other.
    // The interior dihedral angle α satisfies: cos(π - α) = n1 · n2
    // So α = π - acos(n1 · n2).
    let dot = n1.dot(&n2).clamp(-1.0, 1.0);
    let angle = std::f64::consts::PI - dot.acos();
    Ok(angle)
}

/// Collect all unique EdgeIdx values belonging to a solid.
pub fn solid_edges(topo: &TopoStore, solid: SolidIdx) -> Vec<EdgeIdx> {
    let shell_idx = topo.solids.get(solid).outer_shell();
    let faces = &topo.shells.get(shell_idx).faces;
    let mut edge_set = HashSet::new();

    for &face_idx in faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            edge_set.insert(topo.half_edges.get(he).edge);
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }

    edge_set.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box_builder::make_box_into;
    use crate::cylinder_builder::make_cylinder_into;
    use crate::sphere_builder::make_sphere_into;
    use rustkernel_math::Point3;

    fn make_test_box() -> (TopoStore, AnalyticalGeomStore, SolidIdx) {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        (topo, geom, solid)
    }

    #[test]
    fn test_box_edge_adjacency() {
        let (topo, _geom, solid) = make_test_box();
        let edges = solid_edges(&topo, solid);
        assert_eq!(edges.len(), 12, "Box should have 12 unique edges");

        for &edge_idx in &edges {
            let adj = edge_adjacency(&topo, edge_idx).unwrap();
            assert_ne!(
                adj.face_a, adj.face_b,
                "Edge should connect two different faces"
            );
        }
    }

    #[test]
    fn test_box_dihedral_angles() {
        let (topo, geom, solid) = make_test_box();
        let edges = solid_edges(&topo, solid);
        let half_pi = std::f64::consts::FRAC_PI_2;

        for &edge_idx in &edges {
            let angle = dihedral_angle(&topo, &geom, edge_idx).unwrap();
            assert!(
                (angle - half_pi).abs() < 1e-10,
                "Box dihedral should be PI/2, got {}",
                angle
            );
        }
    }

    #[test]
    fn test_box_all_convex() {
        let (topo, geom, solid) = make_test_box();
        let edges = solid_edges(&topo, solid);

        for &edge_idx in &edges {
            let conv = edge_convexity(&topo, &geom, edge_idx).unwrap();
            assert_eq!(
                conv,
                EdgeConvexity::Convex,
                "All box edges should be convex"
            );
        }
    }

    #[test]
    fn test_solid_edges_count() {
        let (topo, _geom, solid) = make_test_box();
        let edges = solid_edges(&topo, solid);
        assert_eq!(edges.len(), 12, "Box should have 12 unique edges");
    }

    // ── Cylinder tests (curved edges between non-planar faces) ──

    fn make_test_cylinder() -> (TopoStore, AnalyticalGeomStore, SolidIdx) {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cylinder_into(
            &mut topo, &mut geom,
            Point3::origin(), 2.0, 5.0, 16,
        );
        (topo, geom, solid)
    }

    #[test]
    fn test_cylinder_edges_all_convex() {
        let (topo, geom, solid) = make_test_cylinder();
        let edges = solid_edges(&topo, solid);

        for &edge_idx in &edges {
            let conv = edge_convexity(&topo, &geom, edge_idx).unwrap();
            assert_eq!(
                conv,
                EdgeConvexity::Convex,
                "All cylinder edges should be convex, edge {:?}",
                edge_idx
            );
        }
    }

    #[test]
    fn test_cylinder_dihedral_angles() {
        let (topo, geom, solid) = make_test_cylinder();
        let edges = solid_edges(&topo, solid);

        for &edge_idx in &edges {
            let angle = dihedral_angle(&topo, &geom, edge_idx).unwrap();
            // All cylinder edges have dihedral angle PI/2 (cap-to-side)
            // or close to PI (side-to-side on a 16-segment polygon).
            assert!(
                angle > 0.0 && angle < std::f64::consts::PI + 0.01,
                "Dihedral angle out of range: {} for edge {:?}",
                angle, edge_idx
            );
        }
    }

    #[test]
    fn test_cylinder_cap_side_dihedral_is_right_angle() {
        // A cylinder with many segments should have cap-to-side edges
        // with dihedral angle close to PI/2.
        let (topo, geom, solid) = make_test_cylinder();
        let edges = solid_edges(&topo, solid);

        let mut found_right_angle = false;
        let half_pi = std::f64::consts::FRAC_PI_2;
        for &edge_idx in &edges {
            let angle = dihedral_angle(&topo, &geom, edge_idx).unwrap();
            if (angle - half_pi).abs() < 0.15 {
                // Within ~9° of 90° — this is a cap-to-side edge
                found_right_angle = true;
            }
        }
        assert!(found_right_angle, "Should find cap-to-side edges near PI/2");
    }

    #[test]
    fn test_surface_normal_at_point_plane() {
        let (topo, geom, solid) = make_test_box();
        let edges = solid_edges(&topo, solid);
        let edge = edges[0];
        let adj = edge_adjacency(&topo, edge).unwrap();

        // For a planar face, surface_normal_at_point should match plane_normal
        let mid = edge_midpoint(&topo, &geom, edge);
        let n_general = surface_normal_at_point(&geom, &topo, adj.face_a, &mid);
        let n_plane = plane_normal(&geom, &topo, adj.face_a);
        assert!(
            (n_general - n_plane).norm() < 1e-10,
            "surface_normal_at_point should match plane_normal for planar faces"
        );
    }

    #[test]
    fn test_sphere_all_edges_convex() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_sphere_into(&mut topo, &mut geom, Point3::origin(), 3.0, 8, 6);
        let edges = solid_edges(&topo, solid);

        for &edge_idx in &edges {
            let conv = edge_convexity(&topo, &geom, edge_idx).unwrap();
            assert_eq!(
                conv,
                EdgeConvexity::Convex,
                "All sphere edges should be convex, edge {:?}",
                edge_idx
            );
        }
    }
}

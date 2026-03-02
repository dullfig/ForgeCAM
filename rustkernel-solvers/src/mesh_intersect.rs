use rustkernel_math::tri_tri::tri_tri_intersect;
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::mesh_cache::FaceMesh;
use tracing::debug;

/// Find all intersection segments between two triangle meshes.
///
/// For each triangle in `mesh_a` × each triangle in `mesh_b`, runs triangle-triangle
/// intersection and collects the resulting segments. O(n*m) brute force — adequate
/// for the moderate tessellation densities used in the kernel (typically 16×16).
pub fn mesh_mesh_intersect(mesh_a: &FaceMesh, mesh_b: &FaceMesh) -> Vec<(Point3, Point3)> {
    let mut segments = Vec::new();

    for ti_a in 0..mesh_a.triangle_count() {
        let ia0 = mesh_a.indices[ti_a * 3] as usize;
        let ia1 = mesh_a.indices[ti_a * 3 + 1] as usize;
        let ia2 = mesh_a.indices[ti_a * 3 + 2] as usize;
        let a0 = &mesh_a.positions[ia0];
        let a1 = &mesh_a.positions[ia1];
        let a2 = &mesh_a.positions[ia2];

        for ti_b in 0..mesh_b.triangle_count() {
            let ib0 = mesh_b.indices[ti_b * 3] as usize;
            let ib1 = mesh_b.indices[ti_b * 3 + 1] as usize;
            let ib2 = mesh_b.indices[ti_b * 3 + 2] as usize;
            let b0 = &mesh_b.positions[ib0];
            let b1 = &mesh_b.positions[ib1];
            let b2 = &mesh_b.positions[ib2];

            let result = tri_tri_intersect(a0, a1, a2, b0, b1, b2);
            if let Some((p0, p1)) = result.segment {
                if (p1 - p0).norm() > 1e-14 {
                    segments.push((p0, p1));
                }
            }
        }
    }

    debug!(segment_count = segments.len(), "mesh_mesh_intersect complete");
    segments
}

/// Find all intersection segments between a triangle mesh and a plane.
///
/// For each triangle in the mesh, finds edges that cross the plane and computes
/// the crossing points. Pairs the crossings per-triangle into segments.
/// This is a fast path for Plane×NURBS intersection.
pub fn mesh_plane_intersect(
    mesh: &FaceMesh,
    plane_origin: &Point3,
    plane_normal: &Vec3,
) -> Vec<(Point3, Point3)> {
    let mut segments = Vec::new();
    let n = plane_normal.normalize();
    let d = n.dot(&plane_origin.coords);

    for ti in 0..mesh.triangle_count() {
        let i0 = mesh.indices[ti * 3] as usize;
        let i1 = mesh.indices[ti * 3 + 1] as usize;
        let i2 = mesh.indices[ti * 3 + 2] as usize;
        let v0 = &mesh.positions[i0];
        let v1 = &mesh.positions[i1];
        let v2 = &mesh.positions[i2];

        let d0 = n.dot(&v0.coords) - d;
        let d1 = n.dot(&v1.coords) - d;
        let d2 = n.dot(&v2.coords) - d;

        let mut crossings = Vec::with_capacity(2);

        edge_plane_crossing(v0, v1, d0, d1, &mut crossings);
        edge_plane_crossing(v1, v2, d1, d2, &mut crossings);
        edge_plane_crossing(v2, v0, d2, d0, &mut crossings);

        // Include vertices exactly on the plane
        if d0.abs() < 1e-12 { crossings.push(*v0); }
        if d1.abs() < 1e-12 { crossings.push(*v1); }
        if d2.abs() < 1e-12 { crossings.push(*v2); }

        // Deduplicate near-coincident points
        crossings.dedup_by(|a, b| (*a - *b).norm() < 1e-12);

        if crossings.len() >= 2 {
            let seg_len = (crossings[1] - crossings[0]).norm();
            if seg_len > 1e-14 {
                segments.push((crossings[0], crossings[1]));
            }
        }
    }

    debug!(segment_count = segments.len(), "mesh_plane_intersect complete");
    segments
}

/// If an edge (va, vb) crosses the plane (signed distances d_a, d_b), push the crossing point.
fn edge_plane_crossing(va: &Point3, vb: &Point3, d_a: f64, d_b: f64, crossings: &mut Vec<Point3>) {
    if (d_a > 1e-12 && d_b < -1e-12) || (d_a < -1e-12 && d_b > 1e-12) {
        let t = d_a / (d_a - d_b);
        let crossing = Point3::from(va.coords * (1.0 - t) + vb.coords * t);
        crossings.push(crossing);
    }
}

/// Chain raw intersection segments into ordered polylines.
///
/// Uses a greedy nearest-endpoint algorithm:
/// 1. Start with the first unused segment.
/// 2. Find the segment whose endpoint is closest to the chain's tail.
/// 3. Append (flipping if necessary) and continue.
/// 4. Also try extending the chain's head.
/// 5. Start a new chain when no segment is within tolerance.
///
/// Returns multiple polylines (NURBS SSI can produce disjoint curves).
pub fn chain_segments(segments: Vec<(Point3, Point3)>, tolerance: f64) -> Vec<Vec<Point3>> {
    if segments.is_empty() {
        return Vec::new();
    }

    let mut used = vec![false; segments.len()];
    let mut chains: Vec<Vec<Point3>> = Vec::new();
    let tol_sq = tolerance * tolerance;

    loop {
        // Find first unused segment to start a new chain.
        let start = used.iter().position(|&u| !u);
        let start = match start {
            Some(i) => i,
            None => break,
        };

        used[start] = true;
        let mut chain = vec![segments[start].0, segments[start].1];

        // Greedily extend the chain in both directions.
        let mut changed = true;
        while changed {
            changed = false;

            // Try to extend the tail
            let tail = *chain.last().unwrap();
            let mut best_idx = None;
            let mut best_dist_sq = tol_sq;
            let mut best_flip = false;

            for (i, seg) in segments.iter().enumerate() {
                if used[i] { continue; }
                let d0 = (seg.0 - tail).norm_squared();
                let d1 = (seg.1 - tail).norm_squared();
                if d0 < best_dist_sq {
                    best_dist_sq = d0;
                    best_idx = Some(i);
                    best_flip = false;
                }
                if d1 < best_dist_sq {
                    best_dist_sq = d1;
                    best_idx = Some(i);
                    best_flip = true;
                }
            }

            if let Some(idx) = best_idx {
                used[idx] = true;
                if best_flip {
                    chain.push(segments[idx].0);
                } else {
                    chain.push(segments[idx].1);
                }
                changed = true;
            }

            // Try to extend the head
            let head = chain[0];
            best_idx = None;
            best_dist_sq = tol_sq;
            best_flip = false;

            for (i, seg) in segments.iter().enumerate() {
                if used[i] { continue; }
                let d0 = (seg.0 - head).norm_squared();
                let d1 = (seg.1 - head).norm_squared();
                if d1 < best_dist_sq {
                    best_dist_sq = d1;
                    best_idx = Some(i);
                    best_flip = false;
                }
                if d0 < best_dist_sq {
                    best_dist_sq = d0;
                    best_idx = Some(i);
                    best_flip = true;
                }
            }

            if let Some(idx) = best_idx {
                used[idx] = true;
                if best_flip {
                    chain.insert(0, segments[idx].1);
                } else {
                    chain.insert(0, segments[idx].0);
                }
                changed = true;
            }
        }

        chains.push(chain);
    }

    debug!(chain_count = chains.len(), "chain_segments complete");
    chains
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_quad_mesh(v0: Point3, v1: Point3, v2: Point3, v3: Point3) -> FaceMesh {
        FaceMesh {
            positions: vec![v0, v1, v2, v3],
            normals: vec![Vec3::z(); 4],
            indices: vec![0, 1, 2, 0, 2, 3],
            uvs: vec![[0.0, 0.0]; 4],
        }
    }

    #[test]
    fn test_mesh_mesh_two_crossing_quads() {
        // Horizontal quad in XY plane at z=0
        let mesh_a = make_quad_mesh(
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
        );
        // Vertical quad in XZ plane at y=0
        let mesh_b = make_quad_mesh(
            Point3::new(-1.0, 0.0, -1.0),
            Point3::new(1.0, 0.0, -1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(-1.0, 0.0, 1.0),
        );

        let segments = mesh_mesh_intersect(&mesh_a, &mesh_b);
        assert!(!segments.is_empty(), "crossing quads should produce segments");

        // Chain them
        let chains = chain_segments(segments, 0.1);
        assert!(!chains.is_empty());
    }

    #[test]
    fn test_mesh_plane_intersect_quad() {
        let mesh = make_quad_mesh(
            Point3::new(-1.0, -1.0, -1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(-1.0, 1.0, 1.0),
        );

        let segments = mesh_plane_intersect(
            &mesh,
            &Point3::origin(),
            &Vec3::new(0.0, 0.0, 1.0),
        );
        assert!(!segments.is_empty(), "tilted quad should cross z=0 plane");
    }

    #[test]
    fn test_chain_segments_simple() {
        let segments = vec![
            (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
            (Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)),
            (Point3::new(2.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)),
        ];
        let chains = chain_segments(segments, 0.01);
        assert_eq!(chains.len(), 1, "should produce one chain");
        assert_eq!(chains[0].len(), 4, "chain should have 4 points");
    }

    #[test]
    fn test_chain_segments_disjoint() {
        let segments = vec![
            (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
            (Point3::new(10.0, 0.0, 0.0), Point3::new(11.0, 0.0, 0.0)),
        ];
        let chains = chain_segments(segments, 0.01);
        assert_eq!(chains.len(), 2, "disjoint segments should produce two chains");
    }

    #[test]
    fn test_chain_segments_reversed() {
        let segments = vec![
            (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
            (Point3::new(2.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)), // reversed
        ];
        let chains = chain_segments(segments, 0.01);
        assert_eq!(chains.len(), 1, "reversed segment should still chain");
        assert_eq!(chains[0].len(), 3);
    }
}

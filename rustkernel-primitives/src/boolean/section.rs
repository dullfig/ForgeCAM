use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::face_util::{face_boundary_points, PlaneFrame};
use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
use rustkernel_topology::intersection::{
    IntersectionCurve, IntersectionPipeline, SurfaceSurfaceResult,
};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use crate::boolean::ops::BooleanError;
use crate::geom::{AnalyticalGeomStore, LineSegment};
use rustkernel_math::polygon2d::Polygon2D;

/// Intersect a solid with a plane, returning the intersection as one or more
/// closed wire loops (LoopIdx in the topology store).
pub fn section_solid(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    pipeline: &IntersectionPipeline,
    solid: SolidIdx,
    plane_origin: Point3,
    plane_normal: Vec3,
) -> Result<Vec<LoopIdx>, BooleanError> {
    let shell_idx = topo.solids.get(solid).shell;
    let faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

    // Add the section plane to the geom store once upfront.
    let section_sid = geom.add_surface(crate::geom::Plane {
        origin: plane_origin,
        normal: plane_normal.normalize(),
    });

    // Collect trimmed segments from intersecting each face with the section plane.
    let mut segments: Vec<(Point3, Point3)> = Vec::new();

    for &face_idx in &faces {
        let sid = topo.faces.get(face_idx).surface_id;
        let face_kind = geom.surface_kind(sid);

        let result = match &face_kind {
            SurfaceKind::Plane { .. } => pipeline.solve(geom, sid, section_sid),
            _ => continue,
        };

        let ssi_result = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        match ssi_result {
            SurfaceSurfaceResult::Empty | SurfaceSurfaceResult::Coincident => continue,
            SurfaceSurfaceResult::Curves(curves) => {
                for curve in curves {
                    match &curve {
                        IntersectionCurve::Line(line) => {
                            // Clip line to the face boundary only (one-sided trim).
                            let frame = {
                                let (fo, fn_) = match face_kind {
                                    SurfaceKind::Plane { origin, normal } => (origin, normal),
                                    _ => continue,
                                };
                                PlaneFrame::from_normal(fo, fn_)
                            };

                            let boundary = face_boundary_points(topo, geom, face_idx);
                            let poly = Polygon2D {
                                vertices: boundary.iter().map(|p| frame.project_to_2d(p)).collect(),
                            };
                            let line_origin_2d = frame.project_to_2d(&line.origin);
                            let line_dir_2d = [
                                line.direction.dot(&frame.u_axis),
                                line.direction.dot(&frame.v_axis),
                            ];

                            let dir_len_sq = line_dir_2d[0].powi(2) + line_dir_2d[1].powi(2);
                            if dir_len_sq < 1e-20 {
                                continue;
                            }

                            let hits = poly.intersect_line(line_origin_2d, line_dir_2d);

                            let mut i = 0;
                            while i + 1 < hits.len() {
                                let t0 = hits[i].t_line;
                                let t1 = hits[i + 1].t_line;
                                if (t1 - t0).abs() > 1e-10 {
                                    let p0 = line.origin + t0 * line.direction;
                                    let p1 = line.origin + t1 * line.direction;
                                    segments.push((p0, p1));
                                }
                                i += 2;
                            }
                        }
                    }
                }
            }
        }
    }

    if segments.is_empty() {
        return Ok(Vec::new());
    }

    // Chain segments into closed loops by matching endpoints.
    let loops = chain_segments_into_loops(topo, geom, &segments);
    Ok(loops)
}

/// Chain line segments into closed loops by matching coincident endpoints.
fn chain_segments_into_loops(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    segments: &[(Point3, Point3)],
) -> Vec<LoopIdx> {
    let tol = 1e-8;
    let tol_sq = tol * tol;

    // Build adjacency: for each endpoint, find connected segments.
    let mut used = vec![false; segments.len()];
    let mut loops = Vec::new();

    for start_idx in 0..segments.len() {
        if used[start_idx] {
            continue;
        }

        // Try to build a chain starting from this segment.
        let mut chain: Vec<(Point3, Point3)> = Vec::new();
        chain.push(segments[start_idx]);
        used[start_idx] = true;

        // Walk forward, matching the end of the chain to the start of another segment.
        loop {
            let current_end = chain.last().unwrap().1;
            let mut found = false;

            for i in 0..segments.len() {
                if used[i] {
                    continue;
                }
                // Check if this segment's start matches current end.
                if (segments[i].0 - current_end).norm_squared() < tol_sq {
                    chain.push(segments[i]);
                    used[i] = true;
                    found = true;
                    break;
                }
                // Check reversed.
                if (segments[i].1 - current_end).norm_squared() < tol_sq {
                    chain.push((segments[i].1, segments[i].0));
                    used[i] = true;
                    found = true;
                    break;
                }
            }

            if !found {
                break;
            }

            // Check if the chain is closed.
            let chain_start = chain[0].0;
            let chain_end = chain.last().unwrap().1;
            if (chain_start - chain_end).norm_squared() < tol_sq {
                break;
            }
        }

        // Check if closed.
        let chain_start = chain[0].0;
        let chain_end = chain.last().unwrap().1;
        if (chain_start - chain_end).norm_squared() >= tol_sq {
            // Not closed — skip (shouldn't happen for a valid solid section).
            continue;
        }

        // Build topology for this closed loop.
        // We don't create a Face — just a Loop with vertices and half-edges.
        let dummy_face = topo.faces.alloc(Face {
            outer_loop: Idx::from_raw(0),
            surface_id: 0,
            mesh_cache: None,
            shell: Idx::from_raw(0),
        });
        let loop_idx = topo.loops.alloc(Loop {
            half_edge: Idx::from_raw(0),
            face: dummy_face,
        });
        topo.faces.get_mut(dummy_face).outer_loop = loop_idx;

        // Create vertices for each segment start point.
        let mut verts = Vec::new();
        for &(start, _) in &chain {
            let pid = geom.add_point(start);
            let vert = topo.vertices.alloc(Vertex { point_id: pid });
            verts.push(vert);
        }

        let n = verts.len();
        let mut he_idxs = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            let start_pt = geom.point(topo.vertices.get(verts[i]).point_id);
            let end_pt = geom.point(topo.vertices.get(verts[j]).point_id);
            let curve_id = geom.add_curve(LineSegment {
                start: start_pt,
                end: end_pt,
            });
            let edge_idx = topo.edges.alloc(Edge {
                half_edges: [Idx::from_raw(0), Idx::from_raw(0)],
                curve_id,
            });
            let he_idx = topo.half_edges.alloc(HalfEdge {
                origin: verts[i],
                twin: None,
                next: Idx::from_raw(0),
                edge: edge_idx,
                loop_ref: loop_idx,
            });
            topo.edges.get_mut(edge_idx).half_edges[0] = he_idx;
            he_idxs.push(he_idx);
        }

        for i in 0..n {
            let next = he_idxs[(i + 1) % n];
            topo.half_edges.get_mut(he_idxs[i]).next = next;
        }
        topo.loops.get_mut(loop_idx).half_edge = he_idxs[0];

        loops.push(loop_idx);
    }

    loops
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::Kernel;

    #[test]
    fn test_section_box_through_middle() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 2.0, 2.0);

        // Section at z=0 (through the middle of the box).
        let loops = section_solid(
            &mut k.topo,
            &mut k.geom,
            &k.pipeline,
            solid,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        )
        .unwrap();

        assert_eq!(loops.len(), 1, "Expected 1 section loop through middle of box");

        // Count edges in the loop.
        let loop_idx = loops[0];
        let start_he = k.topo.loops.get(loop_idx).half_edge;
        let mut count = 0;
        let mut he = start_he;
        loop {
            count += 1;
            he = k.topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
        assert_eq!(count, 4, "Section of a box should be a rectangle (4 edges)");
    }

    #[test]
    fn test_section_box_miss() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 2.0, 2.0);

        // Section at z=5 (above the box).
        let loops = section_solid(
            &mut k.topo,
            &mut k.geom,
            &k.pipeline,
            solid,
            Point3::new(0.0, 0.0, 5.0),
            Vec3::new(0.0, 0.0, 1.0),
        )
        .unwrap();

        assert_eq!(loops.len(), 0, "Section plane above box should give no loops");
    }
}

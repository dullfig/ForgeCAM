use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{
    AnalyticalGeomStore, CylinderSurface, LineSegment, Plane, SurfaceDef,
};

/// Build a cylinder with given center (bottom circle center), radius, height, and N lateral segments.
/// Axis is +Z. Returns the SolidIdx.
///
/// Topology: N lateral quads + top polygon cap + bottom polygon cap.
/// V = 2N, E = 3N, F = N+2 → V-E+F = 2.
pub fn make_cylinder_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    center: Point3,
    radius: f64,
    height: f64,
    n_segments: usize,
) -> SolidIdx {
    let _span = info_span!("make_cylinder", radius, height, n_segments).entered();
    assert!(n_segments >= 3);
    let axis = Vec3::new(0.0, 0.0, 1.0);

    // Generate vertices on bottom and top circles.
    let mut bottom_verts = Vec::with_capacity(n_segments);
    let mut top_verts = Vec::with_capacity(n_segments);

    for i in 0..n_segments {
        let angle = (i as f64 / n_segments as f64) * std::f64::consts::TAU;
        let x = center.x + radius * angle.cos();
        let y = center.y + radius * angle.sin();

        let bottom_p = Point3::new(x, y, center.z);
        let top_p = Point3::new(x, y, center.z + height);

        let bp_id = geom.add_point(bottom_p);
        let tp_id = geom.add_point(top_p);

        let bv = topo.vertices.alloc(Vertex { point_id: bp_id });
        let tv = topo.vertices.alloc(Vertex { point_id: tp_id });

        bottom_verts.push(bv);
        top_verts.push(tv);
    }

    // Create solid and shell.
    let solid_idx = topo.solids.alloc(Solid {
        shells: vec![Idx::from_raw(0)],
        genus: 0,
    });
    let shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: solid_idx,
    });
    topo.solids.get_mut(solid_idx).shells = vec![shell_idx];

    // Surfaces.
    let cyl_surface_id = geom.add_surface(SurfaceDef::Cylinder(CylinderSurface {
        origin: center,
        axis,
        radius,
    }));
    let bottom_surface_id = geom.add_plane(Plane {
        origin: center,
        normal: -axis,
    });
    let top_surface_id = geom.add_plane(Plane {
        origin: Point3::new(center.x, center.y, center.z + height),
        normal: axis,
    });

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    // N lateral quad faces. Each quad uses vertices:
    // bottom[i], bottom[(i+1)%N], top[(i+1)%N], top[i]
    // Winding: outward-facing (viewed from outside).
    for i in 0..n_segments {
        let j = (i + 1) % n_segments;

        let verts = [bottom_verts[i], bottom_verts[j], top_verts[j], top_verts[i]];
        build_face_from_vert_idxs(
            topo,
            geom,
            &verts,
            cyl_surface_id,
            shell_idx,
            &mut he_map,
        );
    }

    // Bottom cap: polygon with vertices bottom[N-1], bottom[N-2], ..., bottom[0] (reversed for outward normal -Z).
    {
        let cap_verts: Vec<VertexIdx> = bottom_verts.iter().rev().copied().collect();
        build_face_from_vert_idxs(topo, geom, &cap_verts, bottom_surface_id, shell_idx, &mut he_map);
    }

    // Top cap: polygon with vertices top[0], top[1], ..., top[N-1] (outward normal +Z).
    {
        let cap_verts: Vec<VertexIdx> = top_verts.clone();
        build_face_from_vert_idxs(topo, geom, &cap_verts, top_surface_id, shell_idx, &mut he_map);
    }

    // Twin matching.
    match_twins_from_map(topo, &he_map);

    solid_idx
}

/// Build a face from vertex indices, registering half-edges in the map.
pub fn build_face_from_vert_idxs(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    verts: &[VertexIdx],
    surface_id: u32,
    shell_idx: ShellIdx,
    he_map: &mut HashMap<(u32, u32), HalfEdgeIdx>,
) -> FaceIdx {
    let n = verts.len();

    let face_idx = topo.faces.alloc(Face {
        outer_loop: Idx::from_raw(0),
        surface_id,
        mesh_cache: None,
        shell: shell_idx,
    });
    let loop_idx = topo.loops.alloc(Loop {
        half_edge: Idx::from_raw(0),
        face: face_idx,
    });
    topo.faces.get_mut(face_idx).outer_loop = loop_idx;
    topo.shells.get_mut(shell_idx).faces.push(face_idx);

    let mut he_idxs = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let origin = verts[i];
        let dest = verts[j];

        let start_pt = geom.points[topo.vertices.get(origin).point_id as usize];
        let end_pt = geom.points[topo.vertices.get(dest).point_id as usize];
        let curve_id = geom.add_line_segment(LineSegment {
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

        let key = (
            topo.vertices.get(origin).point_id,
            topo.vertices.get(dest).point_id,
        );
        he_map.insert(key, he_idx);
    }

    for i in 0..n {
        let next = he_idxs[(i + 1) % n];
        topo.half_edges.get_mut(he_idxs[i]).next = next;
    }
    topo.loops.get_mut(loop_idx).half_edge = he_idxs[0];

    face_idx
}

/// Match twins using the half-edge map.
pub fn match_twins_from_map(
    topo: &mut TopoStore,
    he_map: &HashMap<(u32, u32), HalfEdgeIdx>,
) {
    let keys: Vec<(u32, u32)> = he_map.keys().cloned().collect();
    for (a, b) in keys {
        if let (Some(&he_ab), Some(&he_ba)) = (he_map.get(&(a, b)), he_map.get(&(b, a))) {
            // Skip if already matched (each pair is visited twice: as (a,b) and (b,a))
            if topo.half_edges.get(he_ab).twin.is_some() {
                continue;
            }
            topo.half_edges.get_mut(he_ab).twin = Some(he_ba);
            topo.half_edges.get_mut(he_ba).twin = Some(he_ab);

            let shared_edge = topo.half_edges.get(he_ab).edge;
            topo.half_edges.get_mut(he_ba).edge = shared_edge;
            topo.edges.get_mut(shared_edge).half_edges[1] = he_ba;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    const N: usize = 32;

    fn verify_euler_genus(topo: &TopoStore, solid: SolidIdx) {
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
                let he_data = topo.half_edges.get(he);
                verts.insert(he_data.origin.raw());
                edges.insert(he_data.edge.raw());
                he = he_data.next;
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
            "Euler: V({v}) - E({e}) + F({f}) = {} != {expected} (genus={genus})",
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

    #[test]
    fn test_cylinder_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cylinder_into(&mut topo, &mut geom, Point3::origin(), 1.0, 2.0, N);
        verify_euler_genus(&topo, solid);
    }

    #[test]
    fn test_cylinder_all_twins() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cylinder_into(&mut topo, &mut geom, Point3::origin(), 1.0, 2.0, N);
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_cylinder_vertices_on_surface() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let radius = 3.0;
        let height = 5.0;
        let center = Point3::new(1.0, 2.0, 0.0);
        let solid = make_cylinder_into(&mut topo, &mut geom, center, radius, height, N);

        let shell = topo.solids.get(solid).outer_shell();
        let faces = &topo.shells.get(shell).faces;
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert = topo.half_edges.get(he).origin;
                let pid = topo.vertices.get(vert).point_id;
                let p = geom.point(pid);
                let dx = p.x - center.x;
                let dy = p.y - center.y;
                let dist = (dx * dx + dy * dy).sqrt();
                assert!(
                    (dist - radius).abs() < 1e-10,
                    "Vertex not on cylinder: dist={dist}, expected {radius}"
                );
                assert!(
                    p.z >= center.z - 1e-10 && p.z <= center.z + height + 1e-10,
                    "Vertex z={} out of range [{}, {}]",
                    p.z,
                    center.z,
                    center.z + height,
                );
                he = topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }
    }

    #[test]
    fn test_cylinder_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cylinder_into(&mut topo, &mut geom, Point3::origin(), 1.0, 2.0, N);
        let shell = topo.solids.get(solid).outer_shell();
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
    fn test_cylinder_face_counts() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cylinder_into(&mut topo, &mut geom, Point3::origin(), 1.0, 2.0, N);

        let shell = topo.solids.get(solid).outer_shell();
        let n_faces = topo.shells.get(shell).faces.len();
        // N lateral + 1 bottom + 1 top = N + 2
        assert_eq!(n_faces, N + 2);
    }
}

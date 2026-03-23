use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::{info_span, warn};

use rustkernel_geom::{AnalyticalGeomStore, Plane};

use crate::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};

/// Extrude a closed 2D profile along a direction to create a solid.
///
/// `profile` is a closed polygon in 3D (CCW winding when viewed from outside).
/// `direction` is the extrude direction (will be normalized internally).
/// `height` is the extrude distance.
///
/// All surfaces are `Plane`, all curves are `LineSegment` — the existing
/// boolean pipeline works unchanged on extruded solids.
pub fn make_extrude_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    profile: &[Point3],
    direction: Vec3,
    height: f64,
) -> SolidIdx {
    let _span = info_span!("make_extrude", profile_len = profile.len(), height).entered();
    let n = profile.len();
    assert!(n >= 3, "Profile must have at least 3 vertices");
    assert!(height > 0.0, "Height must be positive");

    if direction.norm() < 1e-12 {
        warn!("extrude direction is near-zero length");
    }
    if n < 3 {
        warn!(n, "profile has fewer than 3 points");
    }

    let dir = direction.normalize();
    let offset = dir * height;

    // Create bottom and top vertices.
    let mut bottom_verts = Vec::with_capacity(n);
    let mut top_verts = Vec::with_capacity(n);

    for &pt in profile {
        let bp_id = geom.add_point(pt);
        let tp_id = geom.add_point(pt + offset);
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

    // Surfaces for bottom and top caps.
    let bottom_surface_id = geom.add_plane(Plane {
        origin: profile[0],
        normal: -dir,
    });
    let top_surface_id = geom.add_plane(Plane {
        origin: profile[0] + offset,
        normal: dir,
    });

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    // Side faces: one quad per profile edge.
    // Quad vertices: bottom[i], bottom[(i+1)%n], top[(i+1)%n], top[i]
    // This gives outward-facing winding.
    for i in 0..n {
        let j = (i + 1) % n;

        // Compute side face normal from the edge direction and extrude direction.
        let edge_dir = profile[j] - profile[i];
        let side_normal = edge_dir.cross(&dir).normalize();

        let side_surface_id = geom.add_plane(Plane {
            origin: profile[i],
            normal: side_normal,
        });

        let verts = [bottom_verts[i], bottom_verts[j], top_verts[j], top_verts[i]];
        build_face_from_vert_idxs(topo, geom, &verts, side_surface_id, shell_idx, &mut he_map);
    }

    // Bottom cap: reversed winding for outward normal (-dir).
    {
        let cap_verts: Vec<VertexIdx> = bottom_verts.iter().rev().copied().collect();
        build_face_from_vert_idxs(topo, geom, &cap_verts, bottom_surface_id, shell_idx, &mut he_map);
    }

    // Top cap: forward winding for outward normal (+dir).
    {
        let cap_verts: Vec<VertexIdx> = top_verts.clone();
        build_face_from_vert_idxs(topo, geom, &cap_verts, top_surface_id, shell_idx, &mut he_map);
    }

    // Twin matching.
    match_twins_from_map(topo, &he_map);

    solid_idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    fn verify_euler(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).outer_shell();
        let faces = &topo.shells.get(shell_idx).faces;

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
        assert_eq!(v - e + f, 2, "Euler: V({v}) - E({e}) + F({f}) != 2");
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
    fn test_extrude_rectangle_is_box() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        // Rectangle on XY plane.
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];

        let solid = make_extrude_into(&mut topo, &mut geom, &profile, Vec3::new(0.0, 0.0, 1.0), 3.0);

        verify_euler(&topo, solid);
        verify_all_twins(&topo, solid);

        // 4 side + 1 bottom + 1 top = 6 faces (like a box).
        let shell = topo.solids.get(solid).outer_shell();
        assert_eq!(topo.shells.get(shell).faces.len(), 6);
    }

    #[test]
    fn test_extrude_triangle_prism() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];

        let solid = make_extrude_into(&mut topo, &mut geom, &profile, Vec3::new(0.0, 0.0, 1.0), 2.0);

        verify_euler(&topo, solid);
        verify_all_twins(&topo, solid);

        // 3 side + 1 bottom + 1 top = 5 faces.
        let shell = topo.solids.get(solid).outer_shell();
        assert_eq!(topo.shells.get(shell).faces.len(), 5);
    }

    #[test]
    fn test_extrude_l_shape() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        // L-shaped profile: 6 points.
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];

        let solid = make_extrude_into(&mut topo, &mut geom, &profile, Vec3::new(0.0, 0.0, 1.0), 1.0);

        verify_euler(&topo, solid);
        verify_all_twins(&topo, solid);

        // 6 side + 1 bottom + 1 top = 8 faces.
        let shell = topo.solids.get(solid).outer_shell();
        assert_eq!(topo.shells.get(shell).faces.len(), 8);
    }

    #[test]
    fn test_extrude_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];

        let solid = make_extrude_into(&mut topo, &mut geom, &profile, Vec3::new(0.0, 0.0, 1.0), 1.0);
        let shell = topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut topo, shell, &geom);

        let mut total_tris = 0;
        for &face_idx in &topo.shells.get(shell).faces.clone() {
            let mesh = topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        // 6 quad faces → 12 triangles.
        assert_eq!(total_tris, 12);
    }

    #[test]
    fn test_extrude_non_z_direction() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        // Profile on XZ plane, extrude along Y.
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, 1.0),
        ];

        let solid = make_extrude_into(&mut topo, &mut geom, &profile, Vec3::new(0.0, 1.0, 0.0), 2.0);

        verify_euler(&topo, solid);
        verify_all_twins(&topo, solid);
    }
}

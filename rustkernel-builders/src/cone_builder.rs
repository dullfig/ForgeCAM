use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use rustkernel_geom::{
    AnalyticalGeomStore, ConeSurface, Plane, SurfaceDef,
};

/// Build a cone (or truncated cone/frustum).
///
/// - `center`: center of the bottom circle
/// - `r1`: bottom radius (must be > 0)
/// - `r2`: top radius (0 for pointed cone, > 0 for frustum)
/// - `height`: height along +Z
/// - `n_segments`: number of lateral segments
///
/// Pointed cone (r2 = 0): N triangular faces + 1 bottom cap.
///   V = N+1, E = 2N, F = N+1 → V-E+F = 2 ✓
///
/// Frustum (r2 > 0): N quad faces + 2 polygon caps.
///   V = 2N, E = 3N, F = N+2 → V-E+F = 2 ✓
pub fn make_cone_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    center: Point3,
    r1: f64,
    r2: f64,
    height: f64,
    n_segments: usize,
) -> SolidIdx {
    assert!(n_segments >= 3);
    assert!(r1 > 0.0);
    assert!(r2 >= 0.0);

    let axis = Vec3::new(0.0, 0.0, 1.0);

    // Compute cone half-angle and apex.
    let half_angle = if (r1 - r2).abs() < 1e-15 {
        // Cylinder — shouldn't get here, but handle gracefully.
        std::f64::consts::FRAC_PI_4
    } else if r2 < 1e-15 {
        // Pointed cone: apex at z = height, half_angle = atan(r1 / height).
        (r1 / height).atan()
    } else {
        // Frustum: apex above top. apex_height = height * r1 / (r1 - r2).
        (r1 / (height * r1 / (r1 - r2))).atan()
    };

    let apex = if r2 < 1e-15 {
        // Pointed cone: physical apex at top.
        Point3::new(center.x, center.y, center.z + height)
    } else {
        // Virtual apex above or below.
        let apex_h = height * r1 / (r1 - r2);
        Point3::new(center.x, center.y, center.z + apex_h)
    };

    let cone_surface_id = geom.add_surface(SurfaceDef::Cone(ConeSurface {
        apex,
        axis,
        half_angle,
    }));

    let bottom_surface_id = geom.add_plane(Plane {
        origin: center,
        normal: -axis,
    });

    // Bottom circle vertices.
    let mut bottom_verts = Vec::with_capacity(n_segments);
    for i in 0..n_segments {
        let angle = (i as f64 / n_segments as f64) * std::f64::consts::TAU;
        let p = Point3::new(
            center.x + r1 * angle.cos(),
            center.y + r1 * angle.sin(),
            center.z,
        );
        let pid = geom.add_point(p);
        let vert = topo.vertices.alloc(Vertex { point_id: pid });
        bottom_verts.push(vert);
    }

    let solid_idx = topo.solids.alloc(Solid {
        shell: Idx::from_raw(0),
        genus: 0,
    });
    let shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: solid_idx,
    });
    topo.solids.get_mut(solid_idx).shell = shell_idx;

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    if r2 < 1e-15 {
        // Pointed cone.
        let apex_pt = Point3::new(center.x, center.y, center.z + height);
        let apex_pid = geom.add_point(apex_pt);
        let apex_vert = topo.vertices.alloc(Vertex { point_id: apex_pid });

        // N triangular lateral faces.
        for i in 0..n_segments {
            let j = (i + 1) % n_segments;
            let verts = [bottom_verts[i], bottom_verts[j], apex_vert];
            build_face(topo, geom, &verts, cone_surface_id, shell_idx, &mut he_map);
        }

        // Bottom cap (reversed winding).
        let cap_verts: Vec<VertexIdx> = bottom_verts.iter().rev().copied().collect();
        build_face(topo, geom, &cap_verts, bottom_surface_id, shell_idx, &mut he_map);
    } else {
        // Frustum (truncated cone).
        let top_surface_id = geom.add_plane(Plane {
            origin: Point3::new(center.x, center.y, center.z + height),
            normal: axis,
        });

        let mut top_verts = Vec::with_capacity(n_segments);
        for i in 0..n_segments {
            let angle = (i as f64 / n_segments as f64) * std::f64::consts::TAU;
            let p = Point3::new(
                center.x + r2 * angle.cos(),
                center.y + r2 * angle.sin(),
                center.z + height,
            );
            let pid = geom.add_point(p);
            let vert = topo.vertices.alloc(Vertex { point_id: pid });
            top_verts.push(vert);
        }

        // N quad lateral faces.
        for i in 0..n_segments {
            let j = (i + 1) % n_segments;
            let verts = [bottom_verts[i], bottom_verts[j], top_verts[j], top_verts[i]];
            build_face(topo, geom, &verts, cone_surface_id, shell_idx, &mut he_map);
        }

        // Bottom cap.
        let cap_bottom: Vec<VertexIdx> = bottom_verts.iter().rev().copied().collect();
        build_face(topo, geom, &cap_bottom, bottom_surface_id, shell_idx, &mut he_map);

        // Top cap.
        build_face(topo, geom, &top_verts, top_surface_id, shell_idx, &mut he_map);
    }

    crate::cylinder_builder::match_twins_from_map(topo, &he_map);

    solid_idx
}

fn build_face(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    verts: &[VertexIdx],
    surface_id: u32,
    shell_idx: ShellIdx,
    he_map: &mut HashMap<(u32, u32), HalfEdgeIdx>,
) -> FaceIdx {
    crate::cylinder_builder::build_face_from_vert_idxs(topo, geom, verts, surface_id, shell_idx, he_map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    const N: usize = 32;

    fn verify_euler(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).shell;
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
        let shell_idx = topo.solids.get(solid).shell;
        for &face_idx in &topo.shells.get(shell_idx).faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                assert!(topo.half_edges.get(he).twin.is_some(), "Unmatched twin");
                he = topo.half_edges.get(he).next;
                if he == start_he { break; }
            }
        }
    }

    #[test]
    fn test_pointed_cone_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cone_into(&mut topo, &mut geom, Point3::origin(), 2.0, 0.0, 3.0, N);
        verify_euler(&topo, solid);
        verify_all_twins(&topo, solid);

        let shell = topo.solids.get(solid).shell;
        assert_eq!(topo.shells.get(shell).faces.len(), N + 1); // N triangles + 1 cap
    }

    #[test]
    fn test_frustum_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_cone_into(&mut topo, &mut geom, Point3::origin(), 3.0, 1.0, 4.0, N);
        verify_euler(&topo, solid);
        verify_all_twins(&topo, solid);

        let shell = topo.solids.get(solid).shell;
        assert_eq!(topo.shells.get(shell).faces.len(), N + 2); // N quads + 2 caps
    }
}

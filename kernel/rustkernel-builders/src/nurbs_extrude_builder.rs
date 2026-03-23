use std::collections::HashMap;

use curvo::prelude::NurbsCurve3D;
use curvo::prelude::NurbsSurface3D;
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{AnalyticalGeomStore, FlippableNurbs, Plane, SurfaceDef};

use crate::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};

/// Extrude a closed NURBS curve into a solid along direction * height.
///
/// The profile curve is sampled at `n_segments` evenly-spaced parameter values.
/// Lateral faces reference a single NURBS extrusion surface.
/// Top and bottom caps are planar (fan triangulation).
///
/// The curve must be closed (start ≈ end within tolerance).
pub fn make_nurbs_extrude_solid_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    curve: &NurbsCurve3D<f64>,
    direction: Vec3,
    height: f64,
    n_segments: usize,
) -> SolidIdx {
    let _span = info_span!("make_nurbs_extrude", height, n_segments).entered();
    assert!(n_segments >= 3, "Need at least 3 segments");
    assert!(height > 0.0, "Height must be positive");

    let dir = direction.normalize();
    let offset = dir * height;

    // Sample the profile curve at n_segments evenly-spaced parameter values.
    let (t_min, t_max) = curve.knots_domain();
    let mut profile_pts: Vec<Point3> = Vec::with_capacity(n_segments);
    for i in 0..n_segments {
        let t = t_min + (t_max - t_min) * i as f64 / n_segments as f64;
        profile_pts.push(curve.point_at(t));
    }

    // Create bottom and top vertices.
    let mut bottom_verts = Vec::with_capacity(n_segments);
    let mut top_verts = Vec::with_capacity(n_segments);

    for &pt in &profile_pts {
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

    // Create the NURBS extrusion surface.
    let nurbs_surface = NurbsSurface3D::<f64>::extrude(curve, &offset);
    let nurbs_surface_id = geom.add_surface(SurfaceDef::Nurbs(FlippableNurbs {
        surface: nurbs_surface,
        flipped: false,
    }));

    // Planar cap surfaces.
    let bottom_surface_id = geom.add_plane(Plane {
        origin: profile_pts[0],
        normal: -dir,
    });
    let top_surface_id = geom.add_plane(Plane {
        origin: profile_pts[0] + offset,
        normal: dir,
    });

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    // N lateral quad faces, all referencing the single NURBS surface.
    for i in 0..n_segments {
        let j = (i + 1) % n_segments;
        let verts = [bottom_verts[i], bottom_verts[j], top_verts[j], top_verts[i]];
        build_face_from_vert_idxs(
            topo,
            geom,
            &verts,
            nurbs_surface_id,
            shell_idx,
            &mut he_map,
        );
    }

    // Bottom cap: reversed winding for outward normal (-dir).
    {
        let cap_verts: Vec<VertexIdx> = bottom_verts.iter().rev().copied().collect();
        build_face_from_vert_idxs(
            topo,
            geom,
            &cap_verts,
            bottom_surface_id,
            shell_idx,
            &mut he_map,
        );
    }

    // Top cap: forward winding for outward normal (+dir).
    {
        let cap_verts: Vec<VertexIdx> = top_verts.clone();
        build_face_from_vert_idxs(
            topo,
            geom,
            &cap_verts,
            top_surface_id,
            shell_idx,
            &mut he_map,
        );
    }

    // Twin matching.
    match_twins_from_map(topo, &he_map);

    solid_idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use curvo::prelude::Interpolation;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    /// Create a closed circular NURBS curve (square approximation).
    fn make_closed_curve() -> NurbsCurve3D<f64> {
        let pts = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0), // close the loop
        ];
        NurbsCurve3D::<f64>::interpolate(&pts, 3).unwrap()
    }

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
    fn test_nurbs_extrude_solid_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve();
        let solid = make_nurbs_extrude_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
            16,
        );
        verify_euler(&topo, solid);
    }

    #[test]
    fn test_nurbs_extrude_solid_twins() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve();
        let solid = make_nurbs_extrude_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
            16,
        );
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_nurbs_extrude_solid_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve();
        let solid = make_nurbs_extrude_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
            16,
        );
        let shell = topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut topo, shell, &geom);

        let mut total_tris = 0;
        for &face_idx in &topo.shells.get(shell).faces.clone() {
            let mesh = topo
                .faces
                .get(face_idx)
                .mesh_cache
                .as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
            // Verify normals are non-degenerate
            for n in &mesh.normals {
                let len = n.norm();
                assert!(len > 0.5, "Degenerate normal: {:?} (len={})", n, len);
            }
        }
        assert!(total_tris > 0, "Should have triangles");
    }

    #[test]
    fn test_nurbs_extrude_solid_face_count() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve();
        let n = 12;
        let solid = make_nurbs_extrude_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            n,
        );
        let shell = topo.solids.get(solid).outer_shell();
        // N lateral + bottom cap + top cap = N + 2
        assert_eq!(topo.shells.get(shell).faces.len(), n + 2);
    }
}

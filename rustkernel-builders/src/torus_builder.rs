use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::info_span;

use rustkernel_geom::{AnalyticalGeomStore, SurfaceDef, TorusSurface};

/// Build a torus with given center, major radius R, minor radius r.
/// Axis is +Z. N × M quad grid on two periodic parameters.
///
/// Topology: V = NM, E = 2NM, F = NM → V-E+F = 0 (genus 1, correct).
pub fn make_torus_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    center: Point3,
    major_radius: f64,
    minor_radius: f64,
    n_major: usize,
    n_minor: usize,
) -> SolidIdx {
    let _span = info_span!("make_torus", major_radius, minor_radius, n_major, n_minor).entered();
    assert!(n_major >= 3);
    assert!(n_minor >= 3);

    use std::f64::consts::TAU;

    let axis = Vec3::new(0.0, 0.0, 1.0);

    let torus_surface_id = geom.add_surface(SurfaceDef::Torus(TorusSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }));

    // Vertex grid[i][j]: i indexes major angle, j indexes minor angle.
    let mut grid: Vec<Vec<VertexIdx>> = Vec::with_capacity(n_major);
    for i in 0..n_major {
        let u = (i as f64 / n_major as f64) * TAU;
        let mut row = Vec::with_capacity(n_minor);
        for j in 0..n_minor {
            let v = (j as f64 / n_minor as f64) * TAU;
            let x = center.x + (major_radius + minor_radius * v.cos()) * u.cos();
            let y = center.y + (major_radius + minor_radius * v.cos()) * u.sin();
            let z = center.z + minor_radius * v.sin();
            let pid = geom.add_point(Point3::new(x, y, z));
            let vert = topo.vertices.alloc(Vertex { point_id: pid });
            row.push(vert);
        }
        grid.push(row);
    }

    let solid_idx = topo.solids.alloc(Solid {
        shell: Idx::from_raw(0),
        genus: 1,
    });
    let shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: solid_idx,
    });
    topo.solids.get_mut(solid_idx).shell = shell_idx;

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    // N_major × N_minor quad faces.
    for i in 0..n_major {
        let i_next = (i + 1) % n_major;
        for j in 0..n_minor {
            let j_next = (j + 1) % n_minor;
            let verts = [
                grid[i][j],
                grid[i_next][j],
                grid[i_next][j_next],
                grid[i][j_next],
            ];
            crate::cylinder_builder::build_face_from_vert_idxs(
                topo, geom, &verts, torus_surface_id, shell_idx, &mut he_map,
            );
        }
    }

    crate::cylinder_builder::match_twins_from_map(topo, &he_map);

    solid_idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    const N_MAJOR: usize = 16;
    const N_MINOR: usize = 8;

    fn verify_euler_genus(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).shell;
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
            "Euler: V({v}) - E({e}) + F({f}) = {} != {expected} (genus {genus})",
            v - e + f
        );
    }

    fn verify_all_twins(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).shell;
        for &face_idx in &topo.shells.get(shell_idx).faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                assert!(topo.half_edges.get(he).twin.is_some(), "Unmatched twin at he={}", he.raw());
                he = topo.half_edges.get(he).next;
                if he == start_he { break; }
            }
        }
    }

    #[test]
    fn test_torus_euler_genus1() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_torus_into(&mut topo, &mut geom, Point3::origin(), 5.0, 1.0, N_MAJOR, N_MINOR);
        assert_eq!(topo.solids.get(solid).genus, 1);
        verify_euler_genus(&topo, solid);
    }

    #[test]
    fn test_torus_all_twins() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_torus_into(&mut topo, &mut geom, Point3::origin(), 5.0, 1.0, N_MAJOR, N_MINOR);
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_torus_vertices_on_surface() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let center = Point3::new(1.0, 2.0, 3.0);
        let big_r = 5.0;
        let small_r = 1.5;
        let solid = make_torus_into(&mut topo, &mut geom, center, big_r, small_r, N_MAJOR, N_MINOR);

        let shell = topo.solids.get(solid).shell;
        let mut checked = HashSet::new();
        for &face_idx in &topo.shells.get(shell).faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert = topo.half_edges.get(he).origin;
                if checked.insert(vert.raw()) {
                    let pid = topo.vertices.get(vert).point_id;
                    let p = geom.point(pid);
                    // Distance from axis (Z through center).
                    let dx = p.x - center.x;
                    let dy = p.y - center.y;
                    let dz = p.z - center.z;
                    let rho = (dx * dx + dy * dy).sqrt(); // distance from axis in XY
                    let d_tube = ((rho - big_r).powi(2) + dz.powi(2)).sqrt();
                    assert!(
                        (d_tube - small_r).abs() < 1e-10,
                        "Vertex not on torus: tube_dist={d_tube}, expected {small_r}"
                    );
                }
                he = topo.half_edges.get(he).next;
                if he == start_he { break; }
            }
        }
    }

    #[test]
    fn test_torus_face_count() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_torus_into(&mut topo, &mut geom, Point3::origin(), 5.0, 1.0, N_MAJOR, N_MINOR);
        let shell = topo.solids.get(solid).shell;
        assert_eq!(topo.shells.get(shell).faces.len(), N_MAJOR * N_MINOR);
    }

    #[test]
    fn test_torus_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_torus_into(&mut topo, &mut geom, Point3::origin(), 5.0, 1.0, N_MAJOR, N_MINOR);
        let shell = topo.solids.get(solid).shell;
        tessellate_shell(&mut topo, shell, &geom);

        let mut total_tris = 0;
        for &face_idx in &topo.shells.get(shell).faces.clone() {
            let mesh = topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert!(total_tris > 0);
    }
}

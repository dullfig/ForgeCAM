use std::collections::HashMap;

use rustkernel_math::Point3;
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use rustkernel_geom::{AnalyticalGeomStore, LineSegment, SphereSurface, SurfaceDef};

/// Build a sphere with given center and radius using N longitude × M latitude grid.
///
/// Topology:
/// - 2 polar vertices + N*(M-1) body vertices → V = 2 + N*(M-1)
/// - N triangles at each pole, N*(M-2) quads in body
/// - F = 2*N + N*(M-2) = N*M
/// - E = N*(2M-1)
/// - V - E + F = 2 + N*(M-1) - N*(2M-1) + N*M = 2 ✓
pub fn make_sphere_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    center: Point3,
    radius: f64,
    n_lon: usize,
    n_lat: usize,
) -> SolidIdx {
    assert!(n_lon >= 3);
    assert!(n_lat >= 2);

    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    // Create the sphere surface.
    let sphere_surface_id = geom.add_surface(SurfaceDef::Sphere(SphereSurface {
        center,
        radius,
    }));

    // South pole (v = -PI/2).
    let south_pole = Point3::new(center.x, center.y, center.z - radius);
    let sp_id = geom.add_point(south_pole);
    let south_vert = topo.vertices.alloc(Vertex { point_id: sp_id });

    // North pole (v = PI/2).
    let north_pole = Point3::new(center.x, center.y, center.z + radius);
    let np_id = geom.add_point(north_pole);
    let north_vert = topo.vertices.alloc(Vertex { point_id: np_id });

    // Body vertices: grid[lat][lon], lat from 1 to n_lat-1.
    let mut grid: Vec<Vec<VertexIdx>> = Vec::new();
    for j in 1..n_lat {
        let v_angle = -FRAC_PI_2 + (j as f64 / n_lat as f64) * PI;
        let mut row = Vec::with_capacity(n_lon);
        for i in 0..n_lon {
            let u_angle = (i as f64 / n_lon as f64) * TAU;
            let x = center.x + radius * v_angle.cos() * u_angle.cos();
            let y = center.y + radius * v_angle.cos() * u_angle.sin();
            let z = center.z + radius * v_angle.sin();
            let pid = geom.add_point(Point3::new(x, y, z));
            let vert = topo.vertices.alloc(Vertex { point_id: pid });
            row.push(vert);
        }
        grid.push(row);
    }

    // Create solid and shell.
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

    // South polar triangles: south_pole, grid[0][(i+1)%N], grid[0][i].
    // Winding must be opposite to body quads' shared edge for twin matching.
    for i in 0..n_lon {
        let j = (i + 1) % n_lon;
        let verts = [south_vert, grid[0][j], grid[0][i]];
        build_face(topo, geom, &verts, sphere_surface_id, shell_idx, &mut he_map);
    }

    // Body quads: grid[lat][i], grid[lat][(i+1)%N], grid[lat+1][(i+1)%N], grid[lat+1][i].
    for lat in 0..grid.len() - 1 {
        for i in 0..n_lon {
            let j = (i + 1) % n_lon;
            let verts = [grid[lat][i], grid[lat][j], grid[lat + 1][j], grid[lat + 1][i]];
            build_face(topo, geom, &verts, sphere_surface_id, shell_idx, &mut he_map);
        }
    }

    // North polar triangles: north_pole, grid[last][i], grid[last][(i+1)%N].
    // Winding must produce edge grid[last][i] → grid[last][j], which is the twin of
    // the body quad's edge grid[last][j] → grid[last][i].
    let last = grid.len() - 1;
    for i in 0..n_lon {
        let j = (i + 1) % n_lon;
        let verts = [north_vert, grid[last][i], grid[last][j]];
        build_face(topo, geom, &verts, sphere_surface_id, shell_idx, &mut he_map);
    }

    // Twin matching.
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

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    const N_LON: usize = 16;
    const N_LAT: usize = 8;

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
            "Euler: V({v}) - E({e}) + F({f}) = {} != {expected}",
            v - e + f
        );
    }

    fn verify_all_twins(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).shell;
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
    fn test_sphere_euler() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_sphere_into(&mut topo, &mut geom, Point3::origin(), 1.0, N_LON, N_LAT);
        verify_euler_genus(&topo, solid);
    }

    #[test]
    fn test_sphere_all_twins() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_sphere_into(&mut topo, &mut geom, Point3::origin(), 1.0, N_LON, N_LAT);
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_sphere_vertices_on_surface() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let radius = 5.0;
        let center = Point3::new(1.0, 2.0, 3.0);
        let solid = make_sphere_into(&mut topo, &mut geom, center, radius, N_LON, N_LAT);

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
                    let dist = (p - center).norm();
                    assert!(
                        (dist - radius).abs() < 1e-10,
                        "Vertex not on sphere: dist={dist}, expected {radius}"
                    );
                }
                he = topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }
    }

    #[test]
    fn test_sphere_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let solid = make_sphere_into(&mut topo, &mut geom, Point3::origin(), 1.0, N_LON, N_LAT);
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

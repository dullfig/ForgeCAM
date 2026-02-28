use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::geom_store::GeomAccess;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use rustkernel_geom::{AnalyticalGeomStore, LineSegment, Plane};

/// Build an axis-aligned box with dimensions (dx, dy, dz) centered at origin.
/// Returns the TopoStore, GeomStore, and the SolidIdx of the created box.
pub fn make_box(dx: f64, dy: f64, dz: f64) -> (TopoStore, AnalyticalGeomStore, SolidIdx) {
    let mut topo = TopoStore::new();
    let mut geom = AnalyticalGeomStore::new();
    let solid_idx = make_box_into(&mut topo, &mut geom, Point3::origin(), dx, dy, dz);
    (topo, geom, solid_idx)
}

/// Build an axis-aligned box with dimensions (dx, dy, dz) centered at `center`,
/// into existing topology and geometry stores.
pub fn make_box_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    center: Point3,
    dx: f64,
    dy: f64,
    dz: f64,
) -> SolidIdx {
    let hx = dx / 2.0;
    let hy = dy / 2.0;
    let hz = dz / 2.0;

    // 8 corner points: bottom face CCW from -x-y-z, then top face.
    //   0: (-hx, -hy, -hz)   4: (-hx, -hy, +hz)
    //   1: (+hx, -hy, -hz)   5: (+hx, -hy, +hz)
    //   2: (+hx, +hy, -hz)   6: (+hx, +hy, +hz)
    //   3: (-hx, +hy, -hz)   7: (-hx, +hy, +hz)
    let cx = center.x;
    let cy = center.y;
    let cz = center.z;
    let corners = [
        Point3::new(cx - hx, cy - hy, cz - hz), // 0
        Point3::new(cx + hx, cy - hy, cz - hz), // 1
        Point3::new(cx + hx, cy + hy, cz - hz), // 2
        Point3::new(cx - hx, cy + hy, cz - hz), // 3
        Point3::new(cx - hx, cy - hy, cz + hz), // 4
        Point3::new(cx + hx, cy - hy, cz + hz), // 5
        Point3::new(cx + hx, cy + hy, cz + hz), // 6
        Point3::new(cx - hx, cy + hy, cz + hz), // 7
    ];

    // Register points in geometry store and create topology vertices.
    let mut vert_idxs = Vec::new();
    for &corner in &corners {
        let point_id = geom.add_point(corner);
        let vi = topo.vertices.alloc(Vertex { point_id });
        vert_idxs.push(vi);
    }

    // 6 faces defined as CCW vertex loops when viewed from outside.
    // Each face: [v0, v1, v2, v3] with outward normal.
    let face_defs: [([usize; 4], Vec3); 6] = [
        ([3, 2, 1, 0], Vec3::new(0.0, 0.0, -1.0)), // bottom (-Z)
        ([4, 5, 6, 7], Vec3::new(0.0, 0.0, 1.0)),  // top (+Z)
        ([0, 1, 5, 4], Vec3::new(0.0, -1.0, 0.0)),  // front (-Y)
        ([2, 3, 7, 6], Vec3::new(0.0, 1.0, 0.0)),   // back (+Y)
        ([3, 0, 4, 7], Vec3::new(-1.0, 0.0, 0.0)),  // left (-X)
        ([1, 2, 6, 5], Vec3::new(1.0, 0.0, 0.0)),   // right (+X)
    ];

    // We need to allocate a solid and shell first so faces can reference them.
    // We'll use placeholder indices and fix them up.
    let solid_idx = topo.solids.alloc(Solid {
        shell: Idx::from_raw(0), // placeholder
        genus: 0,
    });
    let shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: solid_idx,
    });
    topo.solids.get_mut(solid_idx).shell = shell_idx;

    // Map from (origin_vert, dest_vert) -> HalfEdgeIdx, used for twin matching.
    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    for (verts, normal) in &face_defs {
        // Plane origin is the center of the face (center offset along normal direction).
        let face_origin = center + (*normal) * match normal {
            n if n.x.abs() > 0.5 => hx,
            n if n.y.abs() > 0.5 => hy,
            _ => hz,
        };
        let surface_id = geom.add_plane(Plane {
            origin: face_origin,
            normal: *normal,
        });

        // Allocate face and loop with placeholder half-edge index.
        let face_idx = topo.faces.alloc(Face {
            outer_loop: Idx::from_raw(0), // placeholder
            surface_id,
            mesh_cache: None,
            shell: shell_idx,
        });
        let loop_idx = topo.loops.alloc(Loop {
            half_edge: Idx::from_raw(0), // placeholder
            face: face_idx,
        });
        topo.faces.get_mut(face_idx).outer_loop = loop_idx;
        topo.shells.get_mut(shell_idx).faces.push(face_idx);

        // Create 4 half-edges for this face, with placeholder next pointers.
        let mut he_idxs = Vec::with_capacity(4);
        for i in 0..4 {
            let origin = vert_idxs[verts[i]];
            let dest = vert_idxs[verts[(i + 1) % 4]];

            // Create the curve geometry for this edge.
            let start = geom.point(topo.vertices.get(origin).point_id);
            let end = geom.point(topo.vertices.get(dest).point_id);
            let curve_id = geom.add_line_segment(LineSegment { start, end });

            // Allocate edge (may be shared with twin — for M1 we allocate per half-edge pair later).
            let edge_idx = topo.edges.alloc(Edge {
                half_edges: [Idx::from_raw(0), Idx::from_raw(0)], // placeholder
                curve_id,
            });

            let he_idx = topo.half_edges.alloc(HalfEdge {
                origin,
                twin: None,
                next: Idx::from_raw(0), // placeholder
                edge: edge_idx,
                loop_ref: loop_idx,
            });

            topo.edges.get_mut(edge_idx).half_edges[0] = he_idx;
            he_idxs.push(he_idx);

            // Register for twin matching: key is (origin_raw, dest_raw).
            let key = (
                topo.vertices.get(origin).point_id,
                topo.vertices.get(dest).point_id,
            );
            he_map.insert(key, he_idx);
        }

        // Wire up next pointers in a cycle.
        for i in 0..4 {
            let next = he_idxs[(i + 1) % 4];
            topo.half_edges.get_mut(he_idxs[i]).next = next;
        }

        // Set the loop's entry half-edge.
        topo.loops.get_mut(loop_idx).half_edge = he_idxs[0];
    }

    // Twin matching: for each half-edge (a->b), find the half-edge (b->a).
    let keys: Vec<(u32, u32)> = he_map.keys().cloned().collect();
    for (a, b) in keys {
        if let (Some(&he_ab), Some(&he_ba)) = (he_map.get(&(a, b)), he_map.get(&(b, a))) {
            topo.half_edges.get_mut(he_ab).twin = Some(he_ba);
            topo.half_edges.get_mut(he_ba).twin = Some(he_ab);

            // Share the same edge entity: point he_ba's edge to he_ab's edge.
            let shared_edge = topo.half_edges.get(he_ab).edge;
            topo.half_edges.get_mut(he_ba).edge = shared_edge;
            topo.edges.get_mut(shared_edge).half_edges[1] = he_ba;
        }
    }

    solid_idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::tessellate::tessellate_shell;

    #[test]
    fn test_euler_formula() {
        let (topo, _, _) = make_box(1.0, 1.0, 1.0);

        let v = topo.vertices.len();
        let f = topo.faces.len();

        // Count unique edges: each edge entity shared by twin pair.
        // We allocated one edge per half-edge, but twin matching reuses edges.
        // Count edges that are actually referenced.
        let mut edge_set = std::collections::HashSet::new();
        for (_, he) in topo.half_edges.iter() {
            edge_set.insert(he.edge.raw());
        }
        let e = edge_set.len();

        // Euler formula for a closed polyhedron: V - E + F = 2
        assert_eq!(v, 8, "Expected 8 vertices");
        assert_eq!(e, 12, "Expected 12 edges");
        assert_eq!(f, 6, "Expected 6 faces");
        assert_eq!(
            v as i32 - e as i32 + f as i32,
            2,
            "Euler formula V - E + F = 2 failed"
        );
    }

    #[test]
    fn test_all_twins_matched() {
        let (topo, _, _) = make_box(1.0, 1.0, 1.0);

        for (idx, he) in topo.half_edges.iter() {
            assert!(
                he.twin.is_some(),
                "Half-edge {:?} has no twin (box should be closed)",
                idx.raw()
            );
        }
    }

    #[test]
    fn test_tessellation_triangle_count() {
        let (mut topo, geom, solid_idx) = make_box(2.0, 3.0, 4.0);
        let shell_idx = topo.solids.get(solid_idx).shell;

        tessellate_shell(&mut topo, shell_idx, &geom);

        let mut total_triangles = 0;
        for (_, face) in topo.faces.iter() {
            let mesh = face.mesh_cache.as_ref().expect("Face should be tessellated");
            total_triangles += mesh.triangle_count();
        }

        // 6 quad faces, each fan-triangulated into 2 triangles = 12
        assert_eq!(total_triangles, 12, "Expected 12 triangles for a box");
    }
}

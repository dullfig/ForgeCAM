use crate::geom_store::GeomAccess;
use crate::mesh_cache::FaceMesh;
use crate::store::TopoStore;
use crate::topo::*;

/// Tessellate a single face using fan triangulation.
/// Walks the face's outer loop collecting vertices, then fans from v0.
pub fn tessellate_face(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    // Collect vertex positions and normals by walking the loop.
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();

    let mut current = start_he;
    loop {
        let he = topo.half_edges.get(current);
        let vert = topo.vertices.get(he.origin);
        let pos = geom.point(vert.point_id);
        let normal = geom.surface_normal(surface_id, 0.0, 0.0);

        positions.push(pos);
        normals.push(normal);
        uvs.push([0.0, 0.0]); // Placeholder UVs for M1 planar faces

        let next = he.next;
        if next == start_he {
            break;
        }
        current = next;
    }

    // Fan triangulation from vertex 0.
    let n = positions.len() as u32;
    let mut indices = Vec::new();
    for i in 1..n - 1 {
        indices.push(0);
        indices.push(i);
        indices.push(i + 1);
    }

    let mesh = FaceMesh {
        positions,
        normals,
        indices,
        uvs,
    };

    topo.faces.get_mut(face_idx).mesh_cache = Some(mesh);
}

/// Tessellate all faces of a shell.
pub fn tessellate_shell(topo: &mut TopoStore, shell_idx: ShellIdx, geom: &dyn GeomAccess) {
    let face_indices: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for face_idx in face_indices {
        tessellate_face(topo, face_idx, geom);
    }
}

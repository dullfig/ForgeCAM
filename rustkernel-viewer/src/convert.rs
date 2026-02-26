use rustkernel_math::{point3_to_f32, vec3_to_f32};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::ShellIdx;

/// Merged mesh data ready for GPU upload.
pub struct MergedMesh {
    pub positions: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
}

/// Merge all per-face FaceMesh data from a shell into a single mesh.
pub fn merge_shell_meshes(topo: &TopoStore, shell_idx: ShellIdx) -> MergedMesh {
    let shell = topo.shells.get(shell_idx);
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    for &face_idx in &shell.faces {
        let face = topo.faces.get(face_idx);
        if let Some(ref mesh) = face.mesh_cache {
            let base = positions.len() as u32;
            for p in &mesh.positions {
                positions.push(point3_to_f32(p));
            }
            for n in &mesh.normals {
                normals.push(vec3_to_f32(n));
            }
            for &idx in &mesh.indices {
                indices.push(base + idx);
            }
        }
    }

    MergedMesh {
        positions,
        normals,
        indices,
    }
}

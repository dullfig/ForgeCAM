use crate::geom_store::{GeomAccess, SurfaceKind};
use crate::mesh_cache::FaceMesh;
use crate::store::TopoStore;
use crate::topo::*;

/// Options controlling tessellation quality for curved surfaces.
pub struct TessellationOptions {
    /// Minimum number of segments per full circle (2π).
    pub min_segments: usize,
    /// Angular tolerance in radians — maximum angle per segment.
    pub angular_tolerance: f64,
    /// Chordal tolerance — maximum chord-to-arc deviation.
    pub chordal_tolerance: f64,
}

impl Default for TessellationOptions {
    fn default() -> Self {
        Self {
            min_segments: 16,
            angular_tolerance: std::f64::consts::PI / 16.0, // ~11.25°
            chordal_tolerance: 0.01,
        }
    }
}

/// Tessellate a single face, dispatching based on surface kind.
pub fn tessellate_face(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
    tessellate_face_with_options(topo, face_idx, geom, &TessellationOptions::default())
}

/// Tessellate a single face with custom options.
pub fn tessellate_face_with_options(
    topo: &mut TopoStore,
    face_idx: FaceIdx,
    geom: &dyn GeomAccess,
    _opts: &TessellationOptions,
) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let kind = geom.surface_kind(surface_id);

    match kind {
        SurfaceKind::Plane { .. } => tessellate_planar_face(topo, face_idx, geom),
        SurfaceKind::Nurbs => {
            // Check if face has a real boundary loop (solid face) vs standalone surface
            let face = topo.faces.get(face_idx);
            let loop_idx = face.outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let he = topo.half_edges.get(start_he);
            // If the face has actual boundary vertices (more than 1 edge), tessellate via boundary
            if he.next != start_he {
                tessellate_curved_face(topo, face_idx, geom);
            } else if let Some(mesh) = geom.tessellate_surface(surface_id, 16, 16) {
                topo.faces.get_mut(face_idx).mesh_cache = Some(mesh);
            }
        }
        _ => tessellate_curved_face(topo, face_idx, geom),
    }
}

/// Tessellate a planar face using fan triangulation (original Phase 1 method).
fn tessellate_planar_face(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

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
        uvs.push([0.0, 0.0]);

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

/// Tessellate a curved face using its boundary vertices.
/// Each boundary vertex gets a smooth normal from the analytical surface.
/// Uses fan triangulation like planar faces, but with per-vertex normals
/// computed from the analytical surface evaluation.
fn tessellate_curved_face(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();

    let mut current = start_he;
    loop {
        let he = topo.half_edges.get(current);
        let vert = topo.vertices.get(he.origin);
        let pos = geom.point(vert.point_id);

        // Compute the surface normal at this vertex's position using inverse mapping.
        let (u, v) = geom.surface_inverse_uv(surface_id, &pos);
        let normal = geom.surface_normal(surface_id, u, v);

        positions.push(pos);
        normals.push(normal);
        uvs.push([u, v]);

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

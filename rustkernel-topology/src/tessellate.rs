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
    let kind = geom.surface_kind(surface_id);

    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();

    let mut current = start_he;
    loop {
        let he = topo.half_edges.get(current);
        let vert = topo.vertices.get(he.origin);
        let pos = geom.point(vert.point_id);

        // Compute the surface normal at this vertex's position using inverse mapping.
        let (u, v) = inverse_map(&kind, &pos);
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

/// Closed-form inverse mapping from 3D point to (u, v) parametric coordinates for each quadric.
fn inverse_map(kind: &SurfaceKind, p: &rustkernel_math::Point3) -> (f64, f64) {
    use rustkernel_math::orthonormal_basis;

    match kind {
        SurfaceKind::Plane { .. } => (0.0, 0.0),
        SurfaceKind::Cylinder { origin, axis, .. } => {
            let a = axis.normalize();
            let (ref_x, ref_y) = orthonormal_basis(&a);
            let d = p - origin;
            let v = d.dot(&a);
            let u = d.dot(&ref_y).atan2(d.dot(&ref_x));
            (u, v)
        }
        SurfaceKind::Sphere { center, .. } => {
            let d = p - center;
            let axis = rustkernel_math::Vec3::new(0.0, 0.0, 1.0);
            let (ref_x, ref_y) = orthonormal_basis(&axis);
            let r = d.norm();
            if r < 1e-15 {
                return (0.0, 0.0);
            }
            let v = (d.dot(&axis) / r).asin();
            let u = d.dot(&ref_y).atan2(d.dot(&ref_x));
            (u, v)
        }
        SurfaceKind::Cone { apex, axis, .. } => {
            let a = axis.normalize();
            let (ref_x, ref_y) = orthonormal_basis(&a);
            let d = p - apex;
            let v = d.dot(&a);
            let u = d.dot(&ref_y).atan2(d.dot(&ref_x));
            (u, v)
        }
        SurfaceKind::Torus { center, axis, major_radius, .. } => {
            let a = axis.normalize();
            let (ref_x, ref_y) = orthonormal_basis(&a);
            let d = p - center;
            let d_proj_x = d.dot(&ref_x);
            let d_proj_y = d.dot(&ref_y);
            let u = d_proj_y.atan2(d_proj_x);
            let tube_center = *center + *major_radius * (u.cos() * ref_x + u.sin() * ref_y);
            let to_p = p - tube_center;
            let radial = u.cos() * ref_x + u.sin() * ref_y;
            let v = to_p.dot(&a).atan2(to_p.dot(&radial));
            (u, v)
        }
        SurfaceKind::Unknown => (0.0, 0.0),
    }
}

/// Tessellate all faces of a shell.
pub fn tessellate_shell(topo: &mut TopoStore, shell_idx: ShellIdx, geom: &dyn GeomAccess) {
    let face_indices: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for face_idx in face_indices {
        tessellate_face(topo, face_idx, geom);
    }
}

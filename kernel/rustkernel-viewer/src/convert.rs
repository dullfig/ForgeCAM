use std::collections::{HashMap, HashSet};

use rustkernel_math::point3_to_f32;
use rustkernel_topology::geom_store::GeomAccess;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::{FaceIdx, ShellIdx};

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
                normals.push([n.x as f32, n.y as f32, n.z as f32]);
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

/// Compute the outward polygon normal of a face from its vertex loop (Newell's method).
fn face_polygon_normal(topo: &TopoStore, geom: &dyn GeomAccess, face_idx: FaceIdx) -> [f64; 3] {
    let loop_idx = topo.faces.get(face_idx).outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;
    let mut verts = Vec::new();
    let mut cur = start_he;
    loop {
        let he = topo.half_edges.get(cur);
        let p = geom.point(topo.vertices.get(he.origin).point_id);
        verts.push([p.x, p.y, p.z]);
        cur = he.next;
        if cur == start_he {
            break;
        }
    }
    // Newell's method
    let mut nx = 0.0_f64;
    let mut ny = 0.0_f64;
    let mut nz = 0.0_f64;
    let n = verts.len();
    for i in 0..n {
        let j = (i + 1) % n;
        nx += (verts[i][1] - verts[j][1]) * (verts[i][2] + verts[j][2]);
        ny += (verts[i][2] - verts[j][2]) * (verts[i][0] + verts[j][0]);
        nz += (verts[i][0] - verts[j][0]) * (verts[i][1] + verts[j][1]);
    }
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len > 1e-15 {
        [nx / len, ny / len, nz / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Extract unique B-Rep edge line segments from a shell, filtering out smooth edges.
/// Only edges where adjacent faces meet at an angle above the crease threshold are shown.
pub fn extract_shell_edges(
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    shell_idx: ShellIdx,
) -> Vec<([f32; 3], [f32; 3])> {
    let shell = topo.shells.get(shell_idx);
    let crease_cos = (20.0_f64).to_radians().cos(); // ~0.94 — edges sharper than 20° shown

    // Precompute face normals and surface IDs.
    let mut face_normals: HashMap<u32, [f64; 3]> = HashMap::new();
    let mut face_surface: HashMap<u32, u32> = HashMap::new();
    for &face_idx in &shell.faces {
        face_normals.insert(face_idx.raw(), face_polygon_normal(topo, geom, face_idx));
        face_surface.insert(face_idx.raw(), topo.faces.get(face_idx).surface_id);
    }

    // Build map: half-edge → face.
    let mut he_to_face: HashMap<u32, u32> = HashMap::new();
    for &face_idx in &shell.faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut cur = start_he;
        loop {
            he_to_face.insert(cur.raw(), face_idx.raw());
            cur = topo.half_edges.get(cur).next;
            if cur == start_he {
                break;
            }
        }
    }

    let mut seen = HashSet::new();
    let mut lines = Vec::new();

    for &face_idx in &shell.faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut cur = start_he;
        loop {
            let he = topo.half_edges.get(cur);
            if seen.insert(he.edge.raw()) {
                if let Some(twin_idx) = he.twin {
                    let face_a = face_idx.raw();
                    let face_b = he_to_face.get(&twin_idx.raw()).copied();

                    let is_sharp = match face_b {
                        Some(fb) if fb != face_a => {
                            // Same underlying surface → always smooth (e.g. cylinder segments)
                            let surf_a = face_surface.get(&face_a).unwrap();
                            let surf_b = face_surface.get(&fb).unwrap();
                            if surf_a == surf_b {
                                false
                            } else {
                                // Different surfaces → check dihedral angle
                                let na = face_normals.get(&face_a).unwrap();
                                let nb = face_normals.get(&fb).unwrap();
                                let dot = na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2];
                                dot < crease_cos
                            }
                        }
                        _ => true, // boundary or same face → always show
                    };

                    if is_sharp {
                        let p0 = geom.point(topo.vertices.get(he.origin).point_id);
                        let twin = topo.half_edges.get(twin_idx);
                        let p1 = geom.point(topo.vertices.get(twin.origin).point_id);
                        lines.push((point3_to_f32(&p0), point3_to_f32(&p1)));
                    }
                }
            }
            cur = he.next;
            if cur == start_he {
                break;
            }
        }
    }
    lines
}

/// Build a thin triangle strip mesh for rendering edges as visible lines.
/// Creates two perpendicular strips per edge so edges are visible from all angles.
pub fn build_edge_mesh(edges: &[([f32; 3], [f32; 3])]) -> MergedMesh {
    let width: f32 = 0.003;
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    for &(p0, p1) in edges {
        let dx = p1[0] - p0[0];
        let dy = p1[1] - p0[1];
        let dz = p1[2] - p0[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < 1e-8 {
            continue;
        }
        let dir = [dx / len, dy / len, dz / len];

        // Two perpendicular offset directions
        let up = if dir[2].abs() < 0.9 {
            [0.0, 0.0, 1.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let perp1 = normalize3(cross3(dir, up));
        let perp2 = normalize3(cross3(dir, perp1));

        for perp in &[perp1, perp2] {
            let base = positions.len() as u32;
            let off = [perp[0] * width, perp[1] * width, perp[2] * width];
            let n = normalize3(cross3(dir, *perp));

            positions.push([p0[0] - off[0], p0[1] - off[1], p0[2] - off[2]]);
            positions.push([p0[0] + off[0], p0[1] + off[1], p0[2] + off[2]]);
            positions.push([p1[0] + off[0], p1[1] + off[1], p1[2] + off[2]]);
            positions.push([p1[0] - off[0], p1[1] - off[1], p1[2] - off[2]]);

            for _ in 0..4 {
                normals.push(n);
            }

            indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        }
    }

    MergedMesh {
        positions,
        normals,
        indices,
    }
}

/// A candidate mesh edge that may be a silhouette depending on the view.
/// Stores the two endpoints and the normals of the two adjacent triangles.
pub struct SilhouetteCandidate {
    pub p0: [f32; 3],
    pub p1: [f32; 3],
    pub n0: [f32; 3],
    pub n1: [f32; 3],
}

/// Extract silhouette candidate edges from per-face tessellation meshes.
/// Only includes edges that are internal to a single face's tessellation (shared by
/// two triangles within the same face). Boundary edges between B-Rep faces are excluded —
/// those are already handled by the sharp-edge system.
pub fn extract_silhouette_candidates(topo: &TopoStore, shell_idx: ShellIdx) -> Vec<SilhouetteCandidate> {
    let shell = topo.shells.get(shell_idx);
    let mut candidates = Vec::new();

    for &face_idx in &shell.faces {
        let face = topo.faces.get(face_idx);
        let mesh = match face.mesh_cache {
            Some(ref m) => m,
            None => continue,
        };

        let positions: Vec<[f32; 3]> = mesh.positions.iter().map(|p| point3_to_f32(p)).collect();
        let tri_count = mesh.indices.len() / 3;

        // Compute per-triangle normals for this face.
        let mut tri_normals: Vec<[f32; 3]> = Vec::with_capacity(tri_count);
        for t in 0..tri_count {
            let i0 = mesh.indices[t * 3] as usize;
            let i1 = mesh.indices[t * 3 + 1] as usize;
            let i2 = mesh.indices[t * 3 + 2] as usize;
            let a = positions[i0];
            let b = positions[i1];
            let c = positions[i2];
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            tri_normals.push(normalize3(cross3(ab, ac)));
        }

        // Build edge adjacency within this face.
        let mut edge_tris: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for t in 0..tri_count {
            let verts = [
                mesh.indices[t * 3],
                mesh.indices[t * 3 + 1],
                mesh.indices[t * 3 + 2],
            ];
            for &(a, b) in &[(verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                edge_tris.entry(key).or_default().push(t);
            }
        }

        // Only internal edges (shared by exactly 2 triangles within this face).
        // Boundary edges (1 triangle) are the face outline — handled by B-Rep edges.
        for (&(va, vb), tris) in &edge_tris {
            if tris.len() == 2 {
                candidates.push(SilhouetteCandidate {
                    p0: positions[va as usize],
                    p1: positions[vb as usize],
                    n0: tri_normals[tris[0]],
                    n1: tri_normals[tris[1]],
                });
            }
        }
    }
    candidates
}

/// Compute silhouette edges for the current view.
/// Returns line segments where adjacent triangles face opposite sides of the eye.
pub fn compute_silhouette_lines(
    candidates: &[SilhouetteCandidate],
    eye_pos: [f32; 3],
) -> Vec<([f32; 3], [f32; 3])> {
    let mut lines = Vec::new();
    for c in candidates {
        let mid = [
            (c.p0[0] + c.p1[0]) * 0.5,
            (c.p0[1] + c.p1[1]) * 0.5,
            (c.p0[2] + c.p1[2]) * 0.5,
        ];
        let view = [
            eye_pos[0] - mid[0],
            eye_pos[1] - mid[1],
            eye_pos[2] - mid[2],
        ];
        let d0 = dot3(c.n0, view);
        let d1 = dot3(c.n1, view);
        if d0 * d1 <= 0.0 {
            lines.push((c.p0, c.p1));
        }
    }
    lines
}

fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize3(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-10 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

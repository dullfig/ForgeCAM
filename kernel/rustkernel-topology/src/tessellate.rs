use crate::geom_store::{GeomAccess, SurfaceKind};
use crate::mesh_cache::FaceMesh;
use crate::store::TopoStore;
use crate::topo::*;
use tracing::warn;

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

/// Tessellate a single face with default options.
pub fn tessellate_face(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
    tessellate_face_with_options(topo, face_idx, geom, &TessellationOptions::default())
}

/// Tessellate a single face, dispatching based on surface kind.
pub fn tessellate_face_with_options(
    topo: &mut TopoStore,
    face_idx: FaceIdx,
    geom: &dyn GeomAccess,
    opts: &TessellationOptions,
) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let kind = geom.surface_kind(surface_id);

    match kind {
        SurfaceKind::Plane { .. } => tessellate_planar_face(topo, face_idx, geom),
        SurfaceKind::Cylinder { .. }
        | SurfaceKind::Sphere { .. }
        | SurfaceKind::Cone { .. }
        | SurfaceKind::Torus { .. } => {
            tessellate_analytical_face(topo, face_idx, geom, opts);
        }
        SurfaceKind::Nurbs | SurfaceKind::Unknown => {
            let face = topo.faces.get(face_idx);
            let loop_idx = face.outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let he = topo.half_edges.get(start_he);
            if he.next != start_he {
                tessellate_boundary_fan(topo, face_idx, geom);
            } else if let Some(mesh) = geom.tessellate_surface(surface_id, 16, 16) {
                topo.faces.get_mut(face_idx).mesh_cache = Some(mesh);
            }
        }
    }
}

/// Tessellate all faces of a shell with default options.
pub fn tessellate_shell(topo: &mut TopoStore, shell_idx: ShellIdx, geom: &dyn GeomAccess) {
    tessellate_shell_with_options(topo, shell_idx, geom, &TessellationOptions::default())
}

/// Tessellate all faces of a shell with custom options.
pub fn tessellate_shell_with_options(
    topo: &mut TopoStore,
    shell_idx: ShellIdx,
    geom: &dyn GeomAccess,
    opts: &TessellationOptions,
) {
    let face_indices: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();
    for face_idx in face_indices {
        tessellate_face_with_options(topo, face_idx, geom, opts);
    }
}

// ---------------------------------------------------------------------------
// Planar faces — ear-clipping triangulation (handles concave polygons)
// ---------------------------------------------------------------------------

fn tessellate_planar_face(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    let normal = geom.surface_normal(surface_id, 0.0, 0.0);

    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();

    let mut current = start_he;
    loop {
        let he = topo.half_edges.get(current);
        let vert = topo.vertices.get(he.origin);
        let pos = geom.point(vert.point_id);

        positions.push(pos);
        normals.push(normal);
        uvs.push([0.0, 0.0]);

        let next = he.next;
        if next == start_he {
            break;
        }
        current = next;
    }

    let indices = ear_clip_triangulate(&positions, &normal);

    topo.faces.get_mut(face_idx).mesh_cache = Some(FaceMesh {
        positions,
        normals,
        indices,
        uvs,
    });
}

/// Ear-clipping triangulation for a planar polygon.
/// Projects 3D points onto 2D using the face normal, then clips ears.
/// Handles concave polygons correctly.
fn ear_clip_triangulate(positions: &[rustkernel_math::Point3], normal: &rustkernel_math::Vec3) -> Vec<u32> {
    let n = positions.len();
    if n < 3 {
        return Vec::new();
    }
    if n == 3 {
        return vec![0, 1, 2];
    }

    // Build an orthonormal basis on the face plane for 2D projection
    let (ax_u, ax_v) = ortho_basis(normal);

    // Project to 2D
    let pts2d: Vec<[f64; 2]> = positions
        .iter()
        .map(|p| {
            let d = p - positions[0];
            [d.dot(&ax_u), d.dot(&ax_v)]
        })
        .collect();

    // Ear-clipping on index list
    let mut remaining: Vec<usize> = (0..n).collect();
    let mut indices = Vec::with_capacity((n - 2) * 3);
    let mut safety = n * n; // prevent infinite loop on degenerate input

    while remaining.len() > 3 && safety > 0 {
        safety -= 1;
        let len = remaining.len();
        let mut clipped = false;

        for i in 0..len {
            let prev = remaining[(i + len - 1) % len];
            let curr = remaining[i];
            let next = remaining[(i + 1) % len];

            // Check if this ear has correct winding (convex vertex in polygon)
            if !is_convex_2d(&pts2d[prev], &pts2d[curr], &pts2d[next]) {
                continue;
            }

            // Check no other vertex is inside this triangle
            let mut contains_other = false;
            for j in 0..len {
                let idx = remaining[j];
                if idx == prev || idx == curr || idx == next {
                    continue;
                }
                if point_in_triangle_2d(&pts2d[idx], &pts2d[prev], &pts2d[curr], &pts2d[next]) {
                    contains_other = true;
                    break;
                }
            }

            if !contains_other {
                indices.push(prev as u32);
                indices.push(curr as u32);
                indices.push(next as u32);
                remaining.remove(i);
                clipped = true;
                break;
            }
        }

        if !clipped {
            // Degenerate polygon — fall back to fan from first remaining vertex
            warn!("ear-clip failed to find ear, falling back to fan");
            for i in 1..remaining.len() - 1 {
                indices.push(remaining[0] as u32);
                indices.push(remaining[i] as u32);
                indices.push(remaining[i + 1] as u32);
            }
            break;
        }
    }

    // Final triangle
    if remaining.len() == 3 {
        indices.push(remaining[0] as u32);
        indices.push(remaining[1] as u32);
        indices.push(remaining[2] as u32);
    }

    indices
}

/// Check if vertex B is convex (left turn from A→B→C) in 2D.
fn is_convex_2d(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2]) -> bool {
    let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    cross > 0.0
}

/// Check if point P is inside triangle (A, B, C) in 2D using barycentric coordinates.
fn point_in_triangle_2d(p: &[f64; 2], a: &[f64; 2], b: &[f64; 2], c: &[f64; 2]) -> bool {
    let d = |x: &[f64; 2], y: &[f64; 2], z: &[f64; 2]| -> f64 {
        (x[0] - z[0]) * (y[1] - z[1]) - (y[0] - z[0]) * (x[1] - z[1])
    };

    let d1 = d(p, a, b);
    let d2 = d(p, b, c);
    let d3 = d(p, c, a);

    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

    !(has_neg && has_pos)
}

/// Build an orthonormal basis (u, v) perpendicular to the given normal.
fn ortho_basis(normal: &rustkernel_math::Vec3) -> (rustkernel_math::Vec3, rustkernel_math::Vec3) {
    use rustkernel_math::Vec3;
    let n = normal.normalize();
    // Pick a reference vector not parallel to n
    let ref_vec = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = n.cross(&ref_vec).normalize();
    let v = n.cross(&u);
    (u, v)
}

// ---------------------------------------------------------------------------
// Analytical curved faces — structured (u,v) grid tessellation
// ---------------------------------------------------------------------------

/// Tessellate an analytical curved face using a structured (u,v) grid.
/// Computes grid density from surface curvature and the chordal/angular tolerance.
fn tessellate_analytical_face(
    topo: &mut TopoStore,
    face_idx: FaceIdx,
    geom: &dyn GeomAccess,
    opts: &TessellationOptions,
) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let kind = geom.surface_kind(surface_id);

    // Collect boundary UV values
    let boundary_uvs = collect_boundary_uvs(topo, face_idx, geom);
    if boundary_uvs.len() < 3 {
        warn!(
            face = face_idx.raw(),
            "degenerate boundary (< 3 vertices), skipping tessellation"
        );
        return;
    }

    // Non-rectangular faces (e.g., sphere pole triangles) get boundary-conforming
    // tessellation instead of a rectangular UV grid.  A UV grid on a triangular face
    // extends outside the boundary, creating overlapping geometry with neighbors.
    if boundary_uvs.len() != 4 {
        tessellate_nonrect_analytical(topo, face_idx, geom, opts, &kind);
        return;
    }

    // Compute UV extent with angular wrapping detection
    let (u_min, u_max, v_min, v_max) = compute_uv_bounds(&boundary_uvs, &kind);

    let du = u_max - u_min;
    let dv = v_max - v_min;

    // Degenerate face in parameter space — fall back to boundary fan
    if du < 1e-12 || dv < 1e-12 {
        tessellate_boundary_fan(topo, face_idx, geom);
        return;
    }

    // Compute grid density from curvature
    let (n_u, n_v) = compute_grid_density(&kind, du, dv, opts);

    // Generate grid vertices
    let total_verts = (n_u + 1) * (n_v + 1);
    let mut positions = Vec::with_capacity(total_verts);
    let mut normals = Vec::with_capacity(total_verts);
    let mut uvs = Vec::with_capacity(total_verts);

    for j in 0..=n_v {
        let v = v_min + dv * (j as f64 / n_v as f64);
        for i in 0..=n_u {
            let u = u_min + du * (i as f64 / n_u as f64);
            let pos = geom.surface_eval(surface_id, u, v);
            let nml = geom.surface_normal(surface_id, u, v);
            positions.push(pos);
            normals.push(nml);
            uvs.push([u, v]);
        }
    }

    // Generate triangle indices (two triangles per grid cell)
    let mut indices = Vec::with_capacity(n_u * n_v * 6);
    let w = (n_u + 1) as u32;
    for j in 0..n_v as u32 {
        for i in 0..n_u as u32 {
            let bl = j * w + i;
            let br = bl + 1;
            let tl = bl + w;
            let tr = tl + 1;

            indices.push(bl);
            indices.push(br);
            indices.push(tl);

            indices.push(br);
            indices.push(tr);
            indices.push(tl);
        }
    }

    topo.faces.get_mut(face_idx).mesh_cache = Some(FaceMesh {
        positions,
        normals,
        indices,
        uvs,
    });
}

/// Collect (u,v) parameter coordinates for all boundary vertices of a face.
fn collect_boundary_uvs(
    topo: &TopoStore,
    face_idx: FaceIdx,
    geom: &dyn GeomAccess,
) -> Vec<(f64, f64)> {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    let mut uvs = Vec::new();
    let mut current = start_he;
    loop {
        let he = topo.half_edges.get(current);
        let vert = topo.vertices.get(he.origin);
        let pos = geom.point(vert.point_id);
        let (u, v) = geom.surface_inverse_uv(surface_id, &pos);
        uvs.push((u, v));
        current = he.next;
        if current == start_he {
            break;
        }
    }
    uvs
}

/// Compute the (u,v) bounding box of boundary vertices, handling angular wrapping.
///
/// For surfaces with periodic angular parameters (cylinder, sphere, cone, torus),
/// detects when boundary vertices straddle the 0/2π discontinuity and shifts values
/// to produce a continuous range.
fn compute_uv_bounds(boundary_uvs: &[(f64, f64)], kind: &SurfaceKind) -> (f64, f64, f64, f64) {
    use std::f64::consts::{PI, TAU};

    let us: Vec<f64> = boundary_uvs.iter().map(|(u, _)| *u).collect();
    let vs: Vec<f64> = boundary_uvs.iter().map(|(_, v)| *v).collect();

    let (mut u_min, mut u_max) = min_max(&us);
    let (mut v_min, mut v_max) = min_max(&vs);

    // Detect angular wrapping for u parameter (all analytical surfaces have periodic u)
    let u_periodic = matches!(
        kind,
        SurfaceKind::Cylinder { .. }
            | SurfaceKind::Sphere { .. }
            | SurfaceKind::Cone { .. }
            | SurfaceKind::Torus { .. }
    );

    if u_periodic && u_max - u_min > PI {
        // Boundary straddles the 0/2π discontinuity. Shift small values up by 2π.
        let shifted: Vec<f64> = us.iter().map(|&u| if u < PI { u + TAU } else { u }).collect();
        let (smin, smax) = min_max(&shifted);
        u_min = smin;
        u_max = smax;
    }

    // Torus minor angle (v) is also periodic
    if matches!(kind, SurfaceKind::Torus { .. }) && v_max - v_min > PI {
        let shifted: Vec<f64> = vs.iter().map(|&v| if v < PI { v + TAU } else { v }).collect();
        let (smin, smax) = min_max(&shifted);
        v_min = smin;
        v_max = smax;
    }

    (u_min, u_max, v_min, v_max)
}

/// Compute the number of grid subdivisions in u and v for an analytical surface,
/// based on the surface curvature and tessellation tolerances.
fn compute_grid_density(
    kind: &SurfaceKind,
    du: f64,
    dv: f64,
    opts: &TessellationOptions,
) -> (usize, usize) {
    match kind {
        SurfaceKind::Cylinder { radius, .. } => {
            // u = angular (curvature = 1/r), v = linear (no curvature)
            let n_u = angular_segments(radius.abs(), du, opts);
            (n_u, 1)
        }
        SurfaceKind::Sphere { radius, .. } => {
            // Both u (longitude) and v (latitude) have curvature = 1/r
            let r = radius.abs();
            let n_u = angular_segments(r, du, opts);
            let n_v = angular_segments(r, dv, opts);
            (n_u, n_v)
        }
        SurfaceKind::Cone { half_angle, .. } => {
            // u = angular, curvature depends on distance from apex (v).
            // Use the maximum radius in the face's v range for segment count.
            // Cone: radius_at_v = |v * tan(half_angle)|
            // We don't have v_min/v_max here, but du was computed from the full span,
            // so we approximate with a moderate radius.
            let r = half_angle.abs().tan(); // radius at v=1
            let n_u = if r > 1e-12 {
                angular_segments(r, du, opts)
            } else {
                opts.min_segments.max(1)
            };
            (n_u, 1) // v is linear on cone
        }
        SurfaceKind::Torus {
            major_radius,
            minor_radius,
            ..
        } => {
            // u = major angle (curvature ~ 1/(R+r)), v = minor angle (curvature ~ 1/r)
            let n_u =
                angular_segments(major_radius.abs() + minor_radius.abs(), du, opts);
            let n_v = angular_segments(minor_radius.abs(), dv, opts);
            (n_u, n_v)
        }
        _ => {
            // Fallback for unknown surface types
            (opts.min_segments.max(1), opts.min_segments.max(1))
        }
    }
}

/// Compute the number of angular segments for a circular arc of given radius and
/// angular span, such that chord-to-arc deviation stays within tolerance.
///
/// Uses both chordal tolerance (geometry-based) and angular tolerance (display-based),
/// taking the stricter of the two.
fn angular_segments(radius: f64, angle_span: f64, opts: &TessellationOptions) -> usize {
    if radius < 1e-12 || angle_span.abs() < 1e-12 {
        return 1;
    }

    // Chordal tolerance: chord deviation = r * (1 - cos(θ/2))
    // Solve for max θ: θ_max = 2 * acos(1 - tol/r)
    let max_angle_chordal = if opts.chordal_tolerance >= radius {
        std::f64::consts::PI // tolerance larger than radius — clamp
    } else {
        2.0 * (1.0 - opts.chordal_tolerance / radius).acos()
    };

    // Use stricter of chordal and angular tolerance
    let max_angle = max_angle_chordal.min(opts.angular_tolerance);

    let n = (angle_span.abs() / max_angle).ceil() as usize;

    // Scale min_segments proportionally to the angular span
    let min_for_span =
        ((angle_span.abs() / std::f64::consts::TAU) * opts.min_segments as f64).ceil() as usize;

    n.max(min_for_span).max(1)
}

// ---------------------------------------------------------------------------
// Fallback: boundary fan tessellation (for NURBS or degenerate faces)
// ---------------------------------------------------------------------------

/// Fan tessellation from boundary vertices with per-vertex normals.
/// Used as fallback for NURBS faces or faces with degenerate parameter domains.
fn tessellate_boundary_fan(topo: &mut TopoStore, face_idx: FaceIdx, geom: &dyn GeomAccess) {
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
        let (u, v) = geom.surface_inverse_uv(surface_id, &pos);
        let normal = geom.surface_normal(surface_id, u, v);

        if normal.norm() < 0.5 {
            warn!(
                face = face_idx.raw(),
                surface = surface_id,
                "degenerate normal during tessellation"
            );
        }

        positions.push(pos);
        normals.push(normal);
        uvs.push([u, v]);

        let next = he.next;
        if next == start_he {
            break;
        }
        current = next;
    }

    let n = positions.len() as u32;
    let mut indices = Vec::new();
    for i in 1..n - 1 {
        indices.push(0);
        indices.push(i);
        indices.push(i + 1);
    }

    topo.faces.get_mut(face_idx).mesh_cache = Some(FaceMesh {
        positions,
        normals,
        indices,
        uvs,
    });
}

// ---------------------------------------------------------------------------
// Non-rectangular analytical face tessellation
// ---------------------------------------------------------------------------

/// Tessellate a non-rectangular analytical face (e.g., sphere pole triangle).
/// Subdivides boundary edges in 3D with surface projection, then fan-triangulates
/// from a projected centroid.  Avoids the UV-grid approach which produces geometry
/// outside the face boundary when the face is not a quad.
fn tessellate_nonrect_analytical(
    topo: &mut TopoStore,
    face_idx: FaceIdx,
    geom: &dyn GeomAccess,
    opts: &TessellationOptions,
    kind: &SurfaceKind,
) {
    let face = topo.faces.get(face_idx);
    let surface_id = face.surface_id;
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    // 1. Collect boundary 3D positions
    let mut boundary: Vec<rustkernel_math::Point3> = Vec::new();
    let mut current = start_he;
    loop {
        let he = topo.half_edges.get(current);
        let vert = topo.vertices.get(he.origin);
        boundary.push(geom.point(vert.point_id));
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if boundary.len() < 3 {
        return;
    }

    // 2. Build subdivided boundary ring — interpolate in 3D, project onto surface
    let curvature_r = surface_curvature_radius(kind);
    let n_boundary = boundary.len();
    let mut ring: Vec<rustkernel_math::Point3> = Vec::new();

    for i in 0..n_boundary {
        let j = (i + 1) % n_boundary;
        let p0 = boundary[i];
        let p1 = boundary[j];

        // Compute subdivision count from chord length and surface curvature
        let chord = (p1 - p0).norm();
        let n_sub = if curvature_r > 1e-12 && chord > 1e-12 {
            let half = (chord / (2.0 * curvature_r)).min(1.0);
            let angle = 2.0 * half.asin();
            angular_segments(curvature_r, angle, opts)
        } else {
            1
        };

        ring.push(p0);

        for k in 1..n_sub {
            let t = k as f64 / n_sub as f64;
            let approx = rustkernel_math::Point3::new(
                p0.x * (1.0 - t) + p1.x * t,
                p0.y * (1.0 - t) + p1.y * t,
                p0.z * (1.0 - t) + p1.z * t,
            );
            // Project onto surface via inverse-UV → eval
            let (u, v) = geom.surface_inverse_uv(surface_id, &approx);
            let on_surface = geom.surface_eval(surface_id, u, v);
            ring.push(on_surface);
        }
    }

    // 3. Compute centroid of ring, project onto surface
    let n_ring = ring.len();
    let cx = ring.iter().map(|p| p.x).sum::<f64>() / n_ring as f64;
    let cy = ring.iter().map(|p| p.y).sum::<f64>() / n_ring as f64;
    let cz = ring.iter().map(|p| p.z).sum::<f64>() / n_ring as f64;
    let approx_center = rustkernel_math::Point3::new(cx, cy, cz);
    let (cu, cv) = geom.surface_inverse_uv(surface_id, &approx_center);
    let center_pos = geom.surface_eval(surface_id, cu, cv);
    let center_nml = geom.surface_normal(surface_id, cu, cv);

    // 4. Build mesh: center (vertex 0) + ring vertices
    let mut positions = Vec::with_capacity(n_ring + 1);
    let mut normals_out = Vec::with_capacity(n_ring + 1);
    let mut uvs = Vec::with_capacity(n_ring + 1);

    positions.push(center_pos);
    normals_out.push(center_nml);
    uvs.push([cu, cv]);

    for p in &ring {
        let (u, v) = geom.surface_inverse_uv(surface_id, p);
        let nml = geom.surface_normal(surface_id, u, v);
        positions.push(*p);
        normals_out.push(nml);
        uvs.push([u, v]);
    }

    // 5. Fan triangulate from center to ring
    let n_ring_u32 = n_ring as u32;
    let mut indices = Vec::with_capacity(n_ring * 3);
    for i in 0..n_ring_u32 {
        let next = (i + 1) % n_ring_u32;
        indices.push(0);
        indices.push(i + 1);
        indices.push(next + 1);
    }

    topo.faces.get_mut(face_idx).mesh_cache = Some(FaceMesh {
        positions,
        normals: normals_out,
        indices,
        uvs,
    });
}

/// Extract a representative curvature radius from a surface kind.
fn surface_curvature_radius(kind: &SurfaceKind) -> f64 {
    match kind {
        SurfaceKind::Sphere { radius, .. } => radius.abs(),
        SurfaceKind::Cylinder { radius, .. } => radius.abs(),
        SurfaceKind::Cone { half_angle, .. } => half_angle.abs().tan().max(0.01),
        SurfaceKind::Torus { minor_radius, .. } => minor_radius.abs(),
        _ => 1.0,
    }
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

fn min_max(vals: &[f64]) -> (f64, f64) {
    let min = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (min, max)
}

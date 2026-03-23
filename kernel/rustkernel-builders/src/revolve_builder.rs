use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use tracing::{info_span, warn};

use rustkernel_geom::{
    AnalyticalGeomStore, ConeSurface, CylinderSurface, Plane, SurfaceDef,
};

use crate::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};

const EPS: f64 = 1e-10;

// ── Helpers ──

/// Rotate `point` around an axis (defined by `axis_origin` and unit `axis_dir`)
/// by `angle` radians using Rodrigues' rotation formula.
fn rotate_point_around_axis(
    point: Point3,
    axis_origin: Point3,
    axis_dir: Vec3,
    angle: f64,
) -> Point3 {
    let v = point - axis_origin;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let rotated = v * cos_a + axis_dir.cross(&v) * sin_a + axis_dir * axis_dir.dot(&v) * (1.0 - cos_a);
    axis_origin + rotated
}

/// Info about a profile vertex relative to the revolve axis.
struct ProfileVertexInfo {
    point: Point3,
    radius: f64,
    axis_t: f64,
    is_pole: bool,
}

/// Analyze each profile point: compute distance from axis and projection along axis.
fn analyze_profile(
    profile: &[Point3],
    axis_origin: Point3,
    axis_dir: Vec3,
) -> Vec<ProfileVertexInfo> {
    profile
        .iter()
        .map(|&p| {
            let v = p - axis_origin;
            let axis_t = v.dot(&axis_dir);
            let proj = axis_origin + axis_t * axis_dir;
            let radius = (p - proj).norm();
            ProfileVertexInfo {
                point: p,
                radius,
                axis_t,
                is_pole: radius < EPS,
            }
        })
        .collect()
}

/// Classify the swept surface for a profile edge from vertex k to vertex k_next.
fn classify_edge_surface(
    info_k: &ProfileVertexInfo,
    info_kn: &ProfileVertexInfo,
    axis_origin: Point3,
    axis_dir: Vec3,
) -> Option<SurfaceDef> {
    let both_on_axis = info_k.is_pole && info_kn.is_pole;
    if both_on_axis {
        return None; // degenerate — skip
    }

    let same_t = (info_k.axis_t - info_kn.axis_t).abs() < EPS;
    let same_r = (info_k.radius - info_kn.radius).abs() < EPS;

    if same_t {
        // Edge perpendicular to axis → Plane at that height
        let plane_origin = axis_origin + info_k.axis_t * axis_dir;
        // Normal points away from solid interior (outward) — we use +axis if
        // this is the "bottom" or -axis for the "top", but for a swept annular
        // face the normal is just axis_dir (we'll fix orientation via winding).
        // Use the axis direction; the face winding handles outward orientation.
        Some(SurfaceDef::Plane(Plane {
            origin: plane_origin,
            normal: axis_dir,
        }))
    } else if same_r {
        // Edge parallel to axis → Cylinder
        let r = info_k.radius;
        // Cylinder origin: the point on the axis at the lower axis_t
        let t_min = info_k.axis_t.min(info_kn.axis_t);
        let cyl_origin = axis_origin + t_min * axis_dir;
        Some(SurfaceDef::Cylinder(CylinderSurface {
            origin: cyl_origin,
            axis: axis_dir,
            radius: r,
        }))
    } else {
        // General case or one pole → Cone
        // Compute virtual apex by extrapolating the meridional line to r=0.
        let t_k = info_k.axis_t;
        let r_k = info_k.radius;
        let t_kn = info_kn.axis_t;
        let r_kn = info_kn.radius;

        let apex_t = t_k - r_k * (t_kn - t_k) / (r_kn - r_k);
        let apex_point = axis_origin + apex_t * axis_dir;

        // half_angle = atan2(r, |t - apex_t|) using whichever endpoint is off-axis
        let (r_ref, t_ref) = if !info_k.is_pole {
            (r_k, t_k)
        } else {
            (r_kn, t_kn)
        };
        let half_angle = r_ref.atan2((t_ref - apex_t).abs());

        Some(SurfaceDef::Cone(ConeSurface {
            apex: apex_point,
            axis: axis_dir,
            half_angle,
        }))
    }
}

// ── Core algorithm ──

/// Revolve a closed 3D profile polygon around an axis to create a solid.
///
/// - `profile`: closed polygon in 3D (at least 3 vertices)
/// - `axis_origin`: a point on the revolve axis
/// - `axis_dir`: axis direction (will be normalized)
/// - `angle`: sweep angle in radians (use TAU for full revolve)
/// - `n_segments`: angular discretization (e.g. 32)
///
/// When revolving line segments, the swept surfaces are always Plane, Cylinder,
/// or Cone — all already in the geometry store. No new surface types needed.
pub fn make_revolve_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    profile: &[Point3],
    axis_origin: Point3,
    axis_dir: Vec3,
    angle: f64,
    n_segments: usize,
) -> SolidIdx {
    let _span = info_span!("make_revolve", profile_len = profile.len(), angle, n_segments).entered();
    let n_profile = profile.len();
    assert!(n_profile >= 3, "Profile must have at least 3 vertices");
    assert!(n_segments >= 3, "Need at least 3 angular segments");
    assert!(angle > 0.0, "Angle must be positive");

    if axis_dir.norm() < 1e-12 {
        warn!("revolve axis direction is near-zero length");
    }
    if n_profile < 3 {
        warn!(n_profile, "profile has fewer than 3 points");
    }

    let axis_dir = axis_dir.normalize();
    let full_revolve = (angle - std::f64::consts::TAU).abs() < 1e-6;

    // Number of distinct angular positions
    let n_angles = if full_revolve { n_segments } else { n_segments + 1 };

    // Phase A: Profile analysis
    let info = analyze_profile(profile, axis_origin, axis_dir);

    // Phase B: Vertex generation
    // grid[k][j] where k = profile vertex index, j = angular step
    // Pole rows have length 1 (single shared vertex).
    let mut grid: Vec<Vec<VertexIdx>> = Vec::with_capacity(n_profile);
    let mut has_any_pole = false;

    for k in 0..n_profile {
        if info[k].is_pole {
            has_any_pole = true;
            // Single vertex for pole
            let pid = geom.add_point(info[k].point);
            let v = topo.vertices.alloc(Vertex { point_id: pid });
            grid.push(vec![v]);
        } else {
            let mut row = Vec::with_capacity(n_angles);
            for j in 0..n_angles {
                let theta = angle * (j as f64) / (n_segments as f64);
                let rotated = rotate_point_around_axis(
                    info[k].point,
                    axis_origin,
                    axis_dir,
                    theta,
                );
                let pid = geom.add_point(rotated);
                let v = topo.vertices.alloc(Vertex { point_id: pid });
                row.push(v);
            }
            grid.push(row);
        }
    }

    // Create solid and shell
    let solid_idx = topo.solids.alloc(Solid {
        shells: vec![Idx::from_raw(0)],
        genus: 0, // will fix up at end
    });
    let shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: solid_idx,
    });
    topo.solids.get_mut(solid_idx).shells = vec![shell_idx];

    let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

    // Phase C + D: For each profile edge × each angular step, build faces
    for k in 0..n_profile {
        let kn = (k + 1) % n_profile;

        // Classify surface for this profile edge
        let surface_def = classify_edge_surface(&info[k], &info[kn], axis_origin, axis_dir);
        let surface_def = match surface_def {
            Some(s) => s,
            None => continue, // degenerate on-axis edge — skip
        };
        let surface_id = geom.add_surface(surface_def);

        let k_pole = info[k].is_pole;
        let kn_pole = info[kn].is_pole;

        for j in 0..n_segments {
            let j_next = if full_revolve {
                (j + 1) % n_segments
            } else {
                j + 1
            };

            // Get vertex indices, handling pole rows
            let v_k_j = if k_pole { grid[k][0] } else { grid[k][j] };
            let v_k_jn = if k_pole { grid[k][0] } else { grid[k][j_next] };
            let v_kn_j = if kn_pole { grid[kn][0] } else { grid[kn][j] };
            let v_kn_jn = if kn_pole { grid[kn][0] } else { grid[kn][j_next] };

            if k_pole {
                // Triangle fan: pole_k, kn[j_next], kn[j]
                let verts = [v_k_j, v_kn_jn, v_kn_j];
                build_face_from_vert_idxs(
                    topo, geom, &verts, surface_id, shell_idx, &mut he_map,
                );
            } else if kn_pole {
                // Triangle fan: k[j], k[j_next], pole_kn
                let verts = [v_k_j, v_k_jn, v_kn_j];
                build_face_from_vert_idxs(
                    topo, geom, &verts, surface_id, shell_idx, &mut he_map,
                );
            } else {
                // Quad: k[j], k[j_next], kn[j_next], kn[j]
                let verts = [v_k_j, v_k_jn, v_kn_jn, v_kn_j];
                build_face_from_vert_idxs(
                    topo, geom, &verts, surface_id, shell_idx, &mut he_map,
                );
            }
        }
    }

    // Phase E: Caps (partial revolve only)
    //
    // Side face quad [k[j], k[j+1], k+1[j+1], k+1[j]] has edge k+1[0]→k[0] at j=0.
    // For twin matching, the start cap needs edge k[0]→k+1[0], so forward winding.
    // Similarly the end cap needs reversed winding to twin with the last angular step.
    if !full_revolve {
        // Start cap at j=0: forward winding
        {
            let cap_verts: Vec<VertexIdx> = (0..n_profile)
                .map(|k| grid[k][0])
                .collect();
            let cap_normal = compute_profile_plane_normal(profile, axis_origin, axis_dir);
            let start_surface_id = geom.add_plane(Plane {
                origin: geom.points[topo.vertices.get(cap_verts[0]).point_id as usize],
                normal: -cap_normal,
            });
            build_face_from_vert_idxs(
                topo, geom, &cap_verts, start_surface_id, shell_idx, &mut he_map,
            );
        }

        // End cap at j=n_segments: reversed winding
        {
            let cap_verts: Vec<VertexIdx> = (0..n_profile)
                .rev()
                .map(|k| {
                    if info[k].is_pole {
                        grid[k][0]
                    } else {
                        grid[k][n_segments] // last angular position
                    }
                })
                .collect();
            let cap_normal = compute_profile_plane_normal(profile, axis_origin, axis_dir);
            let rotated_normal = rotate_vec_around_axis(cap_normal, axis_dir, angle);
            let end_surface_id = geom.add_plane(Plane {
                origin: geom.points[topo.vertices.get(cap_verts[0]).point_id as usize],
                normal: rotated_normal,
            });
            build_face_from_vert_idxs(
                topo, geom, &cap_verts, end_surface_id, shell_idx, &mut he_map,
            );
        }
    }

    // Twin matching
    match_twins_from_map(topo, &he_map);

    // Phase F: Genus
    let genus = if full_revolve && !has_any_pole {
        1 // torus-like: no poles, full revolve
    } else {
        0 // sphere-like or partial
    };
    topo.solids.get_mut(solid_idx).genus = genus;

    solid_idx
}

/// Compute a normal for the profile plane (the plane containing the profile polygon).
/// This is used for partial revolve cap faces.
fn compute_profile_plane_normal(
    profile: &[Point3],
    axis_origin: Point3,
    axis_dir: Vec3,
) -> Vec3 {
    // The profile lies in a meridional half-plane. The cap normal is perpendicular
    // to the axis and to the radial direction. We compute it as the cross product
    // of two profile edges, but a simpler approach: the cap face is a meridional
    // slice, so its normal is perpendicular to the axis and lies in the tangential
    // direction at angle=0.
    //
    // Find the first non-pole profile point and compute its radial direction.
    let mut radial = Vec3::zeros();
    for p in profile {
        let v = p - axis_origin;
        let along = v.dot(&axis_dir) * axis_dir;
        let r = v - along;
        if r.norm() > EPS {
            radial = r.normalize();
            break;
        }
    }
    // Cap normal = axis × radial (tangential direction)
    axis_dir.cross(&radial).normalize()
}

/// Rotate a vector around an axis by an angle (Rodrigues formula for vectors).
fn rotate_vec_around_axis(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    v * cos_a + axis.cross(&v) * sin_a + axis * axis.dot(&v) * (1.0 - cos_a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    fn verify_euler_genus(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).outer_shell();
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
            "Euler: V({v}) - E({e}) + F({f}) = {} != {expected} (genus={genus})",
            v - e + f
        );
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

    /// Rectangle profile in the XZ meridional plane, offset from axis.
    /// When revolved fully, produces a torus-like solid (genus=1).
    fn rect_profile_off_axis() -> Vec<Point3> {
        // Rectangle: x in [2,3], z in [0,1] — fully off the Y=0 axis (axis = +Z)
        vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
        ]
    }

    /// Right triangle with two vertices on the Z axis (poles).
    /// Full revolve → sphere-like (genus=0).
    fn triangle_two_poles() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),  // on axis (pole)
            Point3::new(1.0, 0.0, 0.5),  // off axis
            Point3::new(0.0, 0.0, 1.0),  // on axis (pole)
        ]
    }

    /// Half-rectangle: one edge on axis + one off-axis vertex.
    /// One on-axis edge, produces cone+cylinder combo, genus=0.
    fn half_rect_one_axis_edge() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),  // on axis
            Point3::new(1.0, 0.0, 0.0),  // off axis
            Point3::new(1.0, 0.0, 1.0),  // off axis
            Point3::new(0.0, 0.0, 1.0),  // on axis
        ]
    }

    const N: usize = 16;

    // ── Test 1: Full revolve of rectangle (no poles) → torus-like, genus=1 ──
    #[test]
    fn test_full_revolve_rect_torus() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = rect_profile_off_axis();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            TAU, N,
        );

        assert_eq!(topo.solids.get(solid).genus, 1);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);
    }

    // ── Test 2: Full revolve of triangle (2 poles) → genus=0 ──
    #[test]
    fn test_full_revolve_triangle_two_poles() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = triangle_two_poles();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            TAU, N,
        );

        assert_eq!(topo.solids.get(solid).genus, 0);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);
    }

    // ── Test 3: Full revolve of half-rect (1 on-axis edge, 2 poles) → genus=0 ──
    #[test]
    fn test_full_revolve_half_rect() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = half_rect_one_axis_edge();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            TAU, N,
        );

        assert_eq!(topo.solids.get(solid).genus, 0);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);
    }

    // ── Test 4: Partial revolve 90° of rectangle → genus=0, has caps ──
    #[test]
    fn test_partial_revolve_90_rect() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = rect_profile_off_axis();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2, N,
        );

        assert_eq!(topo.solids.get(solid).genus, 0);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);

        // Should have 4 edge strips * N quads + 2 caps = 4*N + 2 faces
        let shell = topo.solids.get(solid).outer_shell();
        let n_faces = topo.shells.get(shell).faces.len();
        assert_eq!(n_faces, 4 * N + 2);
    }

    // ── Test 5: Partial revolve 180° → same validation ──
    #[test]
    fn test_partial_revolve_180_rect() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = rect_profile_off_axis();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            PI, N,
        );

        assert_eq!(topo.solids.get(solid).genus, 0);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);
    }

    // ── Test 6: Tessellate a revolved solid → valid mesh, no panics ──
    #[test]
    fn test_revolve_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = rect_profile_off_axis();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            TAU, N,
        );

        let shell = topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut topo, shell, &geom);

        let mut total_tris = 0;
        for &face_idx in &topo.shells.get(shell).faces.clone() {
            let mesh = topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert!(total_tris > 0, "Revolved solid should have triangles");
    }

    // ── Test 7: Edge ∥ axis → Cylinder surface ──
    #[test]
    fn test_surface_classify_parallel_is_cylinder() {
        let axis_origin = Point3::origin();
        let axis_dir = Vec3::new(0.0, 0.0, 1.0);

        let a = ProfileVertexInfo { point: Point3::new(2.0, 0.0, 0.0), radius: 2.0, axis_t: 0.0, is_pole: false };
        let b = ProfileVertexInfo { point: Point3::new(2.0, 0.0, 3.0), radius: 2.0, axis_t: 3.0, is_pole: false };

        let surface = classify_edge_surface(&a, &b, axis_origin, axis_dir).unwrap();
        assert!(matches!(surface, SurfaceDef::Cylinder(_)));
        if let SurfaceDef::Cylinder(c) = &surface {
            assert!((c.radius - 2.0).abs() < EPS);
        }
    }

    // ── Test 8: Edge ⊥ axis → Plane surface ──
    #[test]
    fn test_surface_classify_perp_is_plane() {
        let axis_origin = Point3::origin();
        let axis_dir = Vec3::new(0.0, 0.0, 1.0);

        let a = ProfileVertexInfo { point: Point3::new(1.0, 0.0, 2.0), radius: 1.0, axis_t: 2.0, is_pole: false };
        let b = ProfileVertexInfo { point: Point3::new(3.0, 0.0, 2.0), radius: 3.0, axis_t: 2.0, is_pole: false };

        let surface = classify_edge_surface(&a, &b, axis_origin, axis_dir).unwrap();
        assert!(matches!(surface, SurfaceDef::Plane(_)));
    }

    // ── Test 9: Angled edge → Cone surface ──
    #[test]
    fn test_surface_classify_angled_is_cone() {
        let axis_origin = Point3::origin();
        let axis_dir = Vec3::new(0.0, 0.0, 1.0);

        let a = ProfileVertexInfo { point: Point3::new(2.0, 0.0, 0.0), radius: 2.0, axis_t: 0.0, is_pole: false };
        let b = ProfileVertexInfo { point: Point3::new(1.0, 0.0, 1.0), radius: 1.0, axis_t: 1.0, is_pole: false };

        let surface = classify_edge_surface(&a, &b, axis_origin, axis_dir).unwrap();
        assert!(matches!(surface, SurfaceDef::Cone(_)));
    }

    // ── Test 10: Vertices at correct radius from axis ──
    #[test]
    fn test_revolve_vertices_on_cylinder_radius() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let profile = rect_profile_off_axis();
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            TAU, N,
        );

        let shell = topo.solids.get(solid).outer_shell();
        for &face_idx in &topo.shells.get(shell).faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert = topo.half_edges.get(he).origin;
                let pid = topo.vertices.get(vert).point_id;
                let p = geom.points[pid as usize];
                let r = (p.x * p.x + p.y * p.y).sqrt();
                // All vertices should be at radius 2 or 3 (the rect profile)
                assert!(
                    (r - 2.0).abs() < 0.01 || (r - 3.0).abs() < 0.01,
                    "Vertex at unexpected radius {r}"
                );
                he = topo.half_edges.get(he).next;
                if he == start_he { break; }
            }
        }
    }

    // ── Test 11: Full revolve of L-profile → genus=1 ──
    #[test]
    fn test_full_revolve_l_profile() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        // L-shaped profile, all off-axis
        let profile = vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 1.0),
            Point3::new(3.0, 0.0, 1.0),
            Point3::new(3.0, 0.0, 2.0),
            Point3::new(2.0, 0.0, 2.0),
        ];
        let solid = make_revolve_into(
            &mut topo, &mut geom, &profile,
            Point3::origin(), Vec3::new(0.0, 0.0, 1.0),
            TAU, N,
        );

        assert_eq!(topo.solids.get(solid).genus, 1);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);

        // 6 edges * N faces = 6*N faces (no caps for full revolve)
        let shell = topo.solids.get(solid).outer_shell();
        assert_eq!(topo.shells.get(shell).faces.len(), 6 * N);
    }

    // ── Test 12: Degenerate on-axis edge skipped gracefully ──
    #[test]
    fn test_degenerate_on_axis_edge_skipped() {
        let axis_origin = Point3::origin();
        let axis_dir = Vec3::new(0.0, 0.0, 1.0);

        let a = ProfileVertexInfo { point: Point3::new(0.0, 0.0, 0.0), radius: 0.0, axis_t: 0.0, is_pole: true };
        let b = ProfileVertexInfo { point: Point3::new(0.0, 0.0, 1.0), radius: 0.0, axis_t: 1.0, is_pole: true };

        assert!(classify_edge_surface(&a, &b, axis_origin, axis_dir).is_none());
    }
}

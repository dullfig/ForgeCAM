use std::collections::HashMap;

use curvo::prelude::NurbsCurve3D;
use curvo::prelude::NurbsSurface3D;
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use rustkernel_geom::{AnalyticalGeomStore, Plane, SurfaceDef};

use crate::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};

const EPS: f64 = 1e-10;

/// Rotate a point around an axis using Rodrigues' rotation formula.
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

/// Rotate a vector around an axis (Rodrigues formula for vectors).
fn rotate_vec_around_axis(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    v * cos_a + axis.cross(&v) * sin_a + axis * axis.dot(&v) * (1.0 - cos_a)
}

/// Revolve a NURBS curve around an axis to create a solid.
///
/// The profile curve is sampled at `n_profile_segments` evenly-spaced parameter values.
/// The revolution surface is a NURBS surface created via `curvo::try_revolve`.
/// Vertices on the axis (poles) are shared across angular steps.
///
/// - Full revolve (angle ≈ TAU): no caps, genus=1 if no poles
/// - Partial revolve: planar start/end caps, genus=0
pub fn make_nurbs_revolve_solid_into(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    curve: &NurbsCurve3D<f64>,
    axis_origin: Point3,
    axis_dir: Vec3,
    angle: f64,
    n_profile_segments: usize,
    n_angular_segments: usize,
) -> SolidIdx {
    assert!(n_profile_segments >= 3, "Need at least 3 profile segments");
    assert!(n_angular_segments >= 3, "Need at least 3 angular segments");
    assert!(angle > 0.0, "Angle must be positive");

    let axis_dir = axis_dir.normalize();
    let full_revolve = (angle - std::f64::consts::TAU).abs() < 1e-6;

    // Sample profile curve.
    let (t_min, t_max) = curve.knots_domain();
    let mut profile_pts: Vec<Point3> = Vec::with_capacity(n_profile_segments);
    for i in 0..n_profile_segments {
        let t = t_min + (t_max - t_min) * i as f64 / n_profile_segments as f64;
        profile_pts.push(curve.point_at(t));
    }

    // Analyze profile points: compute radius from axis.
    let profile_info: Vec<(Point3, f64, bool)> = profile_pts
        .iter()
        .map(|&p| {
            let v = p - axis_origin;
            let along = v.dot(&axis_dir);
            let proj = axis_origin + along * axis_dir;
            let radius = (p - proj).norm();
            let is_pole = radius < EPS;
            (p, radius, is_pole)
        })
        .collect();

    let has_any_pole = profile_info.iter().any(|info| info.2);

    // Number of distinct angular positions
    let n_angles = if full_revolve { n_angular_segments } else { n_angular_segments + 1 };

    // Create NURBS revolution surface.
    let nurbs_surface = NurbsSurface3D::<f64>::try_revolve(curve, &axis_origin, &axis_dir, angle)
        .expect("NURBS revolve failed");
    let nurbs_surface_id = geom.add_surface(SurfaceDef::Nurbs(nurbs_surface));

    // Vertex generation: grid[k][j] where k = profile index, j = angular step.
    // Pole rows have length 1.
    let mut grid: Vec<Vec<VertexIdx>> = Vec::with_capacity(n_profile_segments);
    for k in 0..n_profile_segments {
        let (pt, _, is_pole) = profile_info[k];
        if is_pole {
            let pid = geom.add_point(pt);
            let v = topo.vertices.alloc(Vertex { point_id: pid });
            grid.push(vec![v]);
        } else {
            let mut row = Vec::with_capacity(n_angles);
            for j in 0..n_angles {
                let theta = angle * (j as f64) / (n_angular_segments as f64);
                let rotated = rotate_point_around_axis(pt, axis_origin, axis_dir, theta);
                let pid = geom.add_point(rotated);
                let v = topo.vertices.alloc(Vertex { point_id: pid });
                row.push(v);
            }
            grid.push(row);
        }
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

    // Build lateral faces.
    for k in 0..n_profile_segments {
        let kn = (k + 1) % n_profile_segments;
        let k_pole = profile_info[k].2;
        let kn_pole = profile_info[kn].2;

        // Skip degenerate on-axis edges (both endpoints are poles).
        if k_pole && kn_pole {
            continue;
        }

        for j in 0..n_angular_segments {
            let j_next = if full_revolve {
                (j + 1) % n_angular_segments
            } else {
                j + 1
            };

            let v_k_j = if k_pole { grid[k][0] } else { grid[k][j] };
            let v_k_jn = if k_pole { grid[k][0] } else { grid[k][j_next] };
            let v_kn_j = if kn_pole { grid[kn][0] } else { grid[kn][j] };
            let v_kn_jn = if kn_pole { grid[kn][0] } else { grid[kn][j_next] };

            if k_pole {
                // Triangle: pole, kn[j_next], kn[j]
                let verts = [v_k_j, v_kn_jn, v_kn_j];
                build_face_from_vert_idxs(
                    topo, geom, &verts, nurbs_surface_id, shell_idx, &mut he_map,
                );
            } else if kn_pole {
                // Triangle: k[j], k[j_next], pole
                let verts = [v_k_j, v_k_jn, v_kn_j];
                build_face_from_vert_idxs(
                    topo, geom, &verts, nurbs_surface_id, shell_idx, &mut he_map,
                );
            } else {
                // Quad: k[j], k[j_next], kn[j_next], kn[j]
                let verts = [v_k_j, v_k_jn, v_kn_jn, v_kn_j];
                build_face_from_vert_idxs(
                    topo, geom, &verts, nurbs_surface_id, shell_idx, &mut he_map,
                );
            }
        }
    }

    // Caps (partial revolve only).
    if !full_revolve {
        // Start cap at j=0: forward winding.
        {
            let cap_verts: Vec<VertexIdx> = (0..n_profile_segments)
                .map(|k| grid[k][0])
                .collect();
            let cap_normal = compute_profile_plane_normal(&profile_pts, axis_origin, axis_dir);
            let start_surface_id = geom.add_plane(Plane {
                origin: geom.points[topo.vertices.get(cap_verts[0]).point_id as usize],
                normal: -cap_normal,
            });
            build_face_from_vert_idxs(
                topo, geom, &cap_verts, start_surface_id, shell_idx, &mut he_map,
            );
        }

        // End cap at j=n_angular_segments: reversed winding.
        {
            let cap_verts: Vec<VertexIdx> = (0..n_profile_segments)
                .rev()
                .map(|k| {
                    if profile_info[k].2 {
                        grid[k][0]
                    } else {
                        grid[k][n_angular_segments]
                    }
                })
                .collect();
            let cap_normal = compute_profile_plane_normal(&profile_pts, axis_origin, axis_dir);
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

    // Twin matching.
    match_twins_from_map(topo, &he_map);

    // Genus.
    let genus = if full_revolve && !has_any_pole {
        1
    } else {
        0
    };
    topo.solids.get_mut(solid_idx).genus = genus;

    solid_idx
}

/// Compute a normal for the profile plane (meridional half-plane).
fn compute_profile_plane_normal(
    profile: &[Point3],
    axis_origin: Point3,
    axis_dir: Vec3,
) -> Vec3 {
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
    axis_dir.cross(&radial).normalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use curvo::prelude::Interpolation;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;
    use std::f64::consts::{FRAC_PI_2, TAU};

    /// Closed curve off the Z axis — full revolve gives genus=1.
    fn make_closed_curve_off_axis() -> NurbsCurve3D<f64> {
        let pts = vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.5),
            Point3::new(3.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 0.0), // close the loop
        ];
        NurbsCurve3D::<f64>::interpolate(&pts, 3).unwrap()
    }

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
            "Euler: V({v}) - E({e}) + F({f}) = {} != {expected} (genus={genus})",
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
    fn test_nurbs_revolve_solid_full() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve_off_axis();
        let solid = make_nurbs_revolve_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            TAU,
            16,
            16,
        );

        assert_eq!(topo.solids.get(solid).genus, 1);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_nurbs_revolve_solid_partial() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve_off_axis();
        let solid = make_nurbs_revolve_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
            16,
            16,
        );

        assert_eq!(topo.solids.get(solid).genus, 0);
        verify_euler_genus(&topo, solid);
        verify_all_twins(&topo, solid);
    }

    #[test]
    fn test_nurbs_revolve_solid_tessellation() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let curve = make_closed_curve_off_axis();
        let solid = make_nurbs_revolve_solid_into(
            &mut topo,
            &mut geom,
            &curve,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            TAU,
            16,
            16,
        );

        let shell = topo.solids.get(solid).shell;
        tessellate_shell(&mut topo, shell, &geom);

        let mut total_tris = 0;
        for &face_idx in &topo.shells.get(shell).faces.clone() {
            let mesh = topo
                .faces
                .get(face_idx)
                .mesh_cache
                .as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
            for n in &mesh.normals {
                let len = n.norm();
                assert!(len > 0.5, "Degenerate normal: {:?} (len={})", n, len);
            }
        }
        assert!(total_tris > 0, "Revolved solid should have triangles");
    }
}

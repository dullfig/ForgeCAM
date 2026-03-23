use rustkernel_math::{Mat4, Point3, Vec3};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use rustkernel_geom::AnalyticalGeomStore;
use tracing::info_span;

use crate::Kernel;

/// Apply a Mat4 transform to all geometry of an already-copied solid.
/// Walks faces→loop→half-edges, deduplicating by point_id/curve_id/surface_id.
pub(crate) fn transform_solid_geometry(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    m: &Mat4,
    flip_normals: bool,
) {
    use std::collections::HashSet;

    let faces: Vec<FaceIdx> = topo.solid_faces(solid);

    let mut visited_points: HashSet<u32> = HashSet::new();
    let mut visited_curves: HashSet<u32> = HashSet::new();
    let mut visited_surfaces: HashSet<u32> = HashSet::new();

    for face_idx in faces {
        let surface_id = topo.faces.get(face_idx).surface_id;
        if visited_surfaces.insert(surface_id) {
            geom.surfaces[surface_id as usize].apply_transform(m);
            if flip_normals {
                geom.surfaces[surface_id as usize].flip_normal();
            }
        }

        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let vert_idx = topo.half_edges.get(he).origin;
            let point_id = topo.vertices.get(vert_idx).point_id;
            if visited_points.insert(point_id) {
                geom.points[point_id as usize] = m.transform_point(&geom.points[point_id as usize]);
            }

            let edge_idx = topo.half_edges.get(he).edge;
            let curve_id = topo.edges.get(edge_idx).curve_id;
            if visited_curves.insert(curve_id) {
                geom.curves[curve_id as usize].apply_transform(m);
            }

            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }
}

impl Kernel {
    /// Deep-copy a solid and translate all its points by `delta`.
    pub fn translate(&mut self, solid: SolidIdx, delta: Vec3) -> SolidIdx {
        use std::collections::HashSet;

        let (new_solid, evo) = self.copy_solid_with_evolution(solid);

        let faces: Vec<FaceIdx> = self.topo.solid_faces(new_solid);

        let mut visited_points: HashSet<u32> = HashSet::new();
        let mut visited_curves: HashSet<u32> = HashSet::new();
        let mut visited_surfaces: HashSet<u32> = HashSet::new();

        for face_idx in faces {
            let surface_id = self.topo.faces.get(face_idx).surface_id;
            if visited_surfaces.insert(surface_id) {
                self.geom.surfaces[surface_id as usize].translate(&delta);
            }

            let loop_idx = self.topo.faces.get(face_idx).outer_loop;
            let start_he = self.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert_idx = self.topo.half_edges.get(he).origin;
                let point_id = self.topo.vertices.get(vert_idx).point_id;
                if visited_points.insert(point_id) {
                    self.geom.points[point_id as usize] += delta;
                }

                let edge_idx = self.topo.half_edges.get(he).edge;
                let curve_id = self.topo.edges.get(edge_idx).curve_id;
                if visited_curves.insert(curve_id) {
                    self.geom.curves[curve_id as usize].translate(&delta);
                }

                he = self.topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }

        self.last_evolution = Some(evo);
        new_solid
    }

    /// Deep-copy a solid and rotate it around an axis.
    pub fn rotate(
        &mut self,
        solid: SolidIdx,
        axis_origin: Point3,
        axis_dir: Vec3,
        angle: f64,
    ) -> SolidIdx {
        let _span = info_span!("kernel.rotate", solid = solid.raw(), angle).entered();
        let (new_solid, evo) = self.copy_solid_with_evolution(solid);

        let unit_axis = rustkernel_math::nalgebra::Unit::new_normalize(axis_dir);
        let rot = rustkernel_math::nalgebra::Rotation3::from_axis_angle(&unit_axis, angle);
        let t = Mat4::new_translation(&axis_origin.coords);
        let t_inv = Mat4::new_translation(&(-axis_origin.coords));
        let m = t * rot.to_homogeneous() * t_inv;

        transform_solid_geometry(&mut self.topo, &mut self.geom, new_solid, &m, false);
        self.last_evolution = Some(evo);
        new_solid
    }

    /// Deep-copy a solid and scale it uniformly around a center point.
    pub fn scale(&mut self, solid: SolidIdx, center: Point3, factor: f64) -> SolidIdx {
        let _span = info_span!("kernel.scale", solid = solid.raw(), factor).entered();
        let (new_solid, evo) = self.copy_solid_with_evolution(solid);

        let t = Mat4::new_translation(&center.coords);
        let t_inv = Mat4::new_translation(&(-center.coords));
        let s = Mat4::new_nonuniform_scaling(&Vec3::new(factor, factor, factor));
        let m = t * s * t_inv;

        transform_solid_geometry(
            &mut self.topo,
            &mut self.geom,
            new_solid,
            &m,
            factor < 0.0,
        );
        self.last_evolution = Some(evo);
        new_solid
    }

    /// Deep-copy a solid and mirror it across a plane.
    pub fn mirror(
        &mut self,
        solid: SolidIdx,
        plane_origin: Point3,
        plane_normal: Vec3,
    ) -> SolidIdx {
        let _span = info_span!("kernel.mirror", solid = solid.raw()).entered();
        let (new_solid, evo) = self.copy_solid_with_evolution(solid);

        let n = plane_normal.normalize();
        // Householder reflection: I - 2*n*nᵀ
        let reflect_3x3 =
            rustkernel_math::nalgebra::Matrix3::identity() - 2.0 * n * n.transpose();
        let mut reflect = Mat4::identity();
        reflect
            .fixed_view_mut::<3, 3>(0, 0)
            .copy_from(&reflect_3x3);
        let t = Mat4::new_translation(&plane_origin.coords);
        let t_inv = Mat4::new_translation(&(-plane_origin.coords));
        let m = t * reflect * t_inv;

        transform_solid_geometry(&mut self.topo, &mut self.geom, new_solid, &m, true);
        self.last_evolution = Some(evo);
        new_solid
    }

    /// Deep-copy a solid and apply a general 4×4 transform.
    /// Assumes rigid + uniform scale for correct analytical geometry.
    pub fn transform(&mut self, solid: SolidIdx, m: Mat4) -> SolidIdx {
        let _span = info_span!("kernel.transform", solid = solid.raw()).entered();
        let (new_solid, evo) = self.copy_solid_with_evolution(solid);

        let flip = m.fixed_view::<3, 3>(0, 0).determinant() < 0.0;
        transform_solid_geometry(&mut self.topo, &mut self.geom, new_solid, &m, flip);
        self.last_evolution = Some(evo);
        new_solid
    }
}

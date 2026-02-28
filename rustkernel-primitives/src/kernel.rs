use std::collections::HashMap;

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::intersection::IntersectionPipeline;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;

use crate::box_builder::make_box_into;
use crate::geom::{AnalyticalGeomStore, LineSegment, Plane};
use crate::solvers::default_pipeline;

/// Central kernel holding all topology, geometry, and the intersection pipeline.
pub struct Kernel {
    pub(crate) topo: TopoStore,
    pub(crate) geom: AnalyticalGeomStore,
    pub(crate) pipeline: IntersectionPipeline,
}

impl Kernel {
    pub fn new() -> Self {
        Self {
            topo: TopoStore::new(),
            geom: AnalyticalGeomStore::new(),
            pipeline: default_pipeline(),
        }
    }

    /// Build an axis-aligned box centered at the origin.
    pub fn make_box(&mut self, dx: f64, dy: f64, dz: f64) -> SolidIdx {
        make_box_into(&mut self.topo, &mut self.geom, Point3::origin(), dx, dy, dz)
    }

    /// Build an axis-aligned box centered at `center`.
    pub fn make_box_at(&mut self, center: [f64; 3], dx: f64, dy: f64, dz: f64) -> SolidIdx {
        let c = Point3::new(center[0], center[1], center[2]);
        make_box_into(&mut self.topo, &mut self.geom, c, dx, dy, dz)
    }

    /// Deep-copy a solid into fresh arena slots. Returns a new solid that shares
    /// no indices with the original.
    pub fn copy_solid(&mut self, solid: SolidIdx) -> SolidIdx {
        let mut remap = IndexRemap::new();
        copy_solid_impl(&mut self.topo, &mut self.geom, solid, &mut remap)
    }

    /// Deep-copy a solid and translate all its points by `delta`.
    pub fn translate(&mut self, solid: SolidIdx, delta: Vec3) -> SolidIdx {
        use std::collections::HashSet;

        let new_solid = self.copy_solid(solid);

        let shell_idx = self.topo.solids.get(new_solid).shell;
        let faces: Vec<FaceIdx> = self.topo.shells.get(shell_idx).faces.clone();

        let mut visited_points: HashSet<u32> = HashSet::new();
        let mut visited_curves: HashSet<u32> = HashSet::new();
        let mut visited_surfaces: HashSet<u32> = HashSet::new();

        for face_idx in faces {
            let surface_id = self.topo.faces.get(face_idx).surface_id;
            if visited_surfaces.insert(surface_id) {
                self.geom.surfaces[surface_id as usize].origin += delta;
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
                    self.geom.curves[curve_id as usize].start += delta;
                    self.geom.curves[curve_id as usize].end += delta;
                }

                he = self.topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }

        new_solid
    }

    // --- Boolean operations ---

    /// Union of two solids.
    pub fn fuse(
        &mut self,
        a: SolidIdx,
        b: SolidIdx,
    ) -> Result<SolidIdx, crate::boolean::ops::BooleanError> {
        crate::boolean::ops::boolean_op(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            a,
            b,
            crate::boolean::face_selector::BooleanOp::Fuse,
        )
    }

    /// Subtraction: A minus B.
    pub fn cut(
        &mut self,
        a: SolidIdx,
        b: SolidIdx,
    ) -> Result<SolidIdx, crate::boolean::ops::BooleanError> {
        crate::boolean::ops::boolean_op(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            a,
            b,
            crate::boolean::face_selector::BooleanOp::Cut,
        )
    }

    /// Intersection of two solids.
    pub fn common(
        &mut self,
        a: SolidIdx,
        b: SolidIdx,
    ) -> Result<SolidIdx, crate::boolean::ops::BooleanError> {
        crate::boolean::ops::boolean_op(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            a,
            b,
            crate::boolean::face_selector::BooleanOp::Common,
        )
    }

    /// Union of multiple solids (sequential fold).
    pub fn fuse_many(
        &mut self,
        solids: &[SolidIdx],
    ) -> Result<SolidIdx, crate::boolean::ops::BooleanError> {
        if solids.is_empty() {
            return Err(crate::boolean::ops::BooleanError::DegenerateInput(
                "empty list".into(),
            ));
        }
        let mut result = solids[0];
        for &s in &solids[1..] {
            result = self.fuse(result, s)?;
        }
        Ok(result)
    }

    /// Intersect a solid with a plane, returning section wire loops.
    pub fn section(
        &mut self,
        solid: SolidIdx,
        plane_origin: [f64; 3],
        plane_normal: [f64; 3],
    ) -> Result<Vec<LoopIdx>, crate::boolean::ops::BooleanError> {
        crate::boolean::section::section_solid(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            solid,
            Point3::new(plane_origin[0], plane_origin[1], plane_origin[2]),
            Vec3::new(plane_normal[0], plane_normal[1], plane_normal[2]),
        )
    }

    /// Access the topology store (read-only).
    pub fn topo(&self) -> &TopoStore {
        &self.topo
    }

    /// Access the geometry store (read-only).
    pub fn geom(&self) -> &AnalyticalGeomStore {
        &self.geom
    }
}

impl Default for Kernel {
    fn default() -> Self {
        Self::new()
    }
}

/// Tracks old-to-new index mappings during a deep copy.
struct IndexRemap {
    vertices: HashMap<VertexIdx, VertexIdx>,
    half_edges: HashMap<HalfEdgeIdx, HalfEdgeIdx>,
    edges: HashMap<EdgeIdx, EdgeIdx>,
    loops: HashMap<LoopIdx, LoopIdx>,
    faces: HashMap<FaceIdx, FaceIdx>,
    shells: HashMap<ShellIdx, ShellIdx>,
    points: HashMap<u32, u32>,
    curves: HashMap<u32, u32>,
    surfaces: HashMap<u32, u32>,
}

impl IndexRemap {
    fn new() -> Self {
        Self {
            vertices: HashMap::new(),
            half_edges: HashMap::new(),
            edges: HashMap::new(),
            loops: HashMap::new(),
            faces: HashMap::new(),
            shells: HashMap::new(),
            points: HashMap::new(),
            curves: HashMap::new(),
            surfaces: HashMap::new(),
        }
    }
}

/// Deep-copy a solid's full topology + geometry into fresh arena slots.
fn copy_solid_impl(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    remap: &mut IndexRemap,
) -> SolidIdx {
    let shell_idx = topo.solids.get(solid).shell;

    // Collect all face indices first (to avoid borrowing issues).
    let face_idxs: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

    // Pass 1: Copy all vertices + geometry (points) reachable from faces.
    for &face_idx in &face_idxs {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let vert_idx = topo.half_edges.get(he).origin;
            if !remap.vertices.contains_key(&vert_idx) {
                let old_vert = topo.vertices.get(vert_idx);
                let new_point_id = copy_point(geom, old_vert.point_id, &mut remap.points);
                let new_vert_idx = topo.vertices.alloc(Vertex {
                    point_id: new_point_id,
                });
                remap.vertices.insert(vert_idx, new_vert_idx);
            }
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }

    // Pass 2: Copy edges + geometry (curves).
    for &face_idx in &face_idxs {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let edge_idx = topo.half_edges.get(he).edge;
            if !remap.edges.contains_key(&edge_idx) {
                let old_edge = topo.edges.get(edge_idx);
                let new_curve_id = copy_curve(geom, old_edge.curve_id, &mut remap.curves);
                let new_edge_idx = topo.edges.alloc(Edge {
                    half_edges: [Idx::from_raw(0), Idx::from_raw(0)], // fix up later
                    curve_id: new_curve_id,
                });
                remap.edges.insert(edge_idx, new_edge_idx);
            }
            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }

    // Pass 3: Copy faces, loops, and half-edges with remapped indices.
    // Allocate the new solid and shell first.
    let new_solid_idx = topo.solids.alloc(Solid {
        shell: Idx::from_raw(0),
    });
    let new_shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: new_solid_idx,
    });
    topo.solids.get_mut(new_solid_idx).shell = new_shell_idx;
    remap.shells.insert(shell_idx, new_shell_idx);

    for &face_idx in &face_idxs {
        let old_face = topo.faces.get(face_idx);
        let old_surface_id = old_face.surface_id;
        let old_loop_idx = old_face.outer_loop;

        let new_surface_id = copy_surface(geom, old_surface_id, &mut remap.surfaces);

        // Allocate new face and loop with placeholders.
        let new_face_idx = topo.faces.alloc(Face {
            outer_loop: Idx::from_raw(0),
            surface_id: new_surface_id,
            mesh_cache: None,
            shell: new_shell_idx,
        });
        remap.faces.insert(face_idx, new_face_idx);

        let new_loop_idx = topo.loops.alloc(Loop {
            half_edge: Idx::from_raw(0),
            face: new_face_idx,
        });
        remap.loops.insert(old_loop_idx, new_loop_idx);
        topo.faces.get_mut(new_face_idx).outer_loop = new_loop_idx;
        topo.shells.get_mut(new_shell_idx).faces.push(new_face_idx);

        // Copy half-edges for this face's loop.
        let start_he = topo.loops.get(old_loop_idx).half_edge;
        let mut he = start_he;
        let mut new_he_idxs = Vec::new();
        loop {
            let old_origin = topo.half_edges.get(he).origin;
            let old_edge = topo.half_edges.get(he).edge;
            let old_next = topo.half_edges.get(he).next;
            let new_origin = remap.vertices[&old_origin];
            let new_edge = remap.edges[&old_edge];

            let new_he_idx = topo.half_edges.alloc(HalfEdge {
                origin: new_origin,
                twin: None, // fix up below
                next: Idx::from_raw(0), // fix up below
                edge: new_edge,
                loop_ref: new_loop_idx,
            });
            remap.half_edges.insert(he, new_he_idx);
            new_he_idxs.push(new_he_idx);

            he = old_next;
            if he == start_he {
                break;
            }
        }

        // Wire up next pointers.
        let n = new_he_idxs.len();
        for i in 0..n {
            let next = new_he_idxs[(i + 1) % n];
            topo.half_edges.get_mut(new_he_idxs[i]).next = next;
        }
        topo.loops.get_mut(new_loop_idx).half_edge = new_he_idxs[0];

        // Fix up edge half_edges references.
        let mut he = start_he;
        loop {
            let old_he = topo.half_edges.get(he);
            let old_edge_idx = old_he.edge;
            let old_edge = topo.edges.get(old_edge_idx);
            let new_he_idx = remap.half_edges[&he];
            let new_edge_idx = remap.edges[&old_edge_idx];

            if old_edge.half_edges[0] == he {
                topo.edges.get_mut(new_edge_idx).half_edges[0] = new_he_idx;
            } else if old_edge.half_edges[1] == he {
                topo.edges.get_mut(new_edge_idx).half_edges[1] = new_he_idx;
            }

            let next = old_he.next;
            he = next;
            if he == start_he {
                break;
            }
        }
    }

    // Pass 4: Fix up twin pointers.
    // For each old half-edge that had a twin, set the new twin if both were copied.
    let he_pairs: Vec<(HalfEdgeIdx, HalfEdgeIdx)> = remap.half_edges.iter()
        .map(|(&old, &new)| (old, new))
        .collect();
    for (old_he, new_he) in he_pairs {
        if let Some(old_twin) = topo.half_edges.get(old_he).twin {
            if let Some(&new_twin) = remap.half_edges.get(&old_twin) {
                topo.half_edges.get_mut(new_he).twin = Some(new_twin);
            }
        }
    }

    new_solid_idx
}

fn copy_point(geom: &mut AnalyticalGeomStore, old_id: u32, map: &mut HashMap<u32, u32>) -> u32 {
    if let Some(&new_id) = map.get(&old_id) {
        return new_id;
    }
    let p = geom.points[old_id as usize];
    let new_id = geom.add_point(p);
    map.insert(old_id, new_id);
    new_id
}

fn copy_curve(geom: &mut AnalyticalGeomStore, old_id: u32, map: &mut HashMap<u32, u32>) -> u32 {
    if let Some(&new_id) = map.get(&old_id) {
        return new_id;
    }
    let seg = &geom.curves[old_id as usize];
    let new_seg = LineSegment {
        start: seg.start,
        end: seg.end,
    };
    let new_id = geom.add_curve(new_seg);
    map.insert(old_id, new_id);
    new_id
}

fn copy_surface(
    geom: &mut AnalyticalGeomStore,
    old_id: u32,
    map: &mut HashMap<u32, u32>,
) -> u32 {
    if let Some(&new_id) = map.get(&old_id) {
        return new_id;
    }
    let plane = &geom.surfaces[old_id as usize];
    let new_plane = Plane {
        origin: plane.origin,
        normal: plane.normal,
    };
    let new_id = geom.add_surface(new_plane);
    map.insert(old_id, new_id);
    new_id
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;

    fn verify_euler(topo: &TopoStore, solid: SolidIdx) {
        let shell_idx = topo.solids.get(solid).shell;
        let faces = &topo.shells.get(shell_idx).faces;

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
        assert_eq!(v - e + f, 2, "Euler: V({v}) - E({e}) + F({f}) != 2");
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
    fn test_kernel_make_box() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_box(2.0, 3.0, 4.0);
        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);
    }

    #[test]
    fn test_kernel_two_boxes_distinct_indices() {
        let mut kernel = Kernel::new();
        let a = kernel.make_box(1.0, 1.0, 1.0);
        let b = kernel.make_box_at([5.0, 0.0, 0.0], 1.0, 1.0, 1.0);

        assert_ne!(a.raw(), b.raw());

        // Both should independently satisfy Euler.
        verify_euler(&kernel.topo, a);
        verify_euler(&kernel.topo, b);

        // Shells should be different.
        let shell_a = kernel.topo.solids.get(a).shell;
        let shell_b = kernel.topo.solids.get(b).shell;
        assert_ne!(shell_a.raw(), shell_b.raw());

        // Face sets should be disjoint.
        let faces_a: HashSet<u32> = kernel.topo.shells.get(shell_a).faces.iter().map(|f| f.raw()).collect();
        let faces_b: HashSet<u32> = kernel.topo.shells.get(shell_b).faces.iter().map(|f| f.raw()).collect();
        assert!(faces_a.is_disjoint(&faces_b));
    }

    #[test]
    fn test_copy_solid_euler() {
        let mut kernel = Kernel::new();
        let original = kernel.make_box(2.0, 2.0, 2.0);
        let copy = kernel.copy_solid(original);

        assert_ne!(original.raw(), copy.raw());
        verify_euler(&kernel.topo, original);
        verify_euler(&kernel.topo, copy);
        verify_all_twins(&kernel.topo, original);
        verify_all_twins(&kernel.topo, copy);
    }

    #[test]
    fn test_copy_solid_independence() {
        let mut kernel = Kernel::new();
        let original = kernel.make_box(1.0, 1.0, 1.0);

        // Record a vertex position from original.
        let orig_shell = kernel.topo.solids.get(original).shell;
        let orig_face = kernel.topo.shells.get(orig_shell).faces[0];
        let orig_loop = kernel.topo.faces.get(orig_face).outer_loop;
        let orig_he = kernel.topo.loops.get(orig_loop).half_edge;
        let orig_vert = kernel.topo.half_edges.get(orig_he).origin;
        let orig_point_id = kernel.topo.vertices.get(orig_vert).point_id;
        let orig_pos = kernel.geom.points[orig_point_id as usize];

        let copy = kernel.copy_solid(original);

        // Mutate copy's first vertex.
        let copy_shell = kernel.topo.solids.get(copy).shell;
        let copy_face = kernel.topo.shells.get(copy_shell).faces[0];
        let copy_loop = kernel.topo.faces.get(copy_face).outer_loop;
        let copy_he = kernel.topo.loops.get(copy_loop).half_edge;
        let copy_vert = kernel.topo.half_edges.get(copy_he).origin;
        let copy_point_id = kernel.topo.vertices.get(copy_vert).point_id;
        kernel.geom.points[copy_point_id as usize] = Point3::new(99.0, 99.0, 99.0);

        // Original should be unchanged.
        let orig_pos_after = kernel.geom.points[orig_point_id as usize];
        assert_eq!(orig_pos, orig_pos_after, "Original vertex was mutated by copy modification");
    }

    #[test]
    fn test_translate() {
        let mut kernel = Kernel::new();
        let original = kernel.make_box(2.0, 2.0, 2.0);
        let delta = Vec3::new(10.0, 0.0, 0.0);
        let translated = kernel.translate(original, delta);

        verify_euler(&kernel.topo, translated);
        verify_all_twins(&kernel.topo, translated);

        // Collect all vertex positions of the translated solid.
        let shell = kernel.topo.solids.get(translated).shell;
        let faces = kernel.topo.shells.get(shell).faces.clone();
        let mut positions = HashSet::new();
        for face_idx in &faces {
            let loop_idx = kernel.topo.faces.get(*face_idx).outer_loop;
            let start_he = kernel.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert = kernel.topo.half_edges.get(he).origin;
                let pid = kernel.topo.vertices.get(vert).point_id;
                let p = kernel.geom.points[pid as usize];
                // Quantize to avoid floating point issues in HashSet.
                let key = (
                    (p.x * 1000.0).round() as i64,
                    (p.y * 1000.0).round() as i64,
                    (p.z * 1000.0).round() as i64,
                );
                positions.insert(key);
                he = kernel.topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }

        // Original box centered at origin with dx=dy=dz=2: corners at ±1.
        // After translate by (10,0,0): corners at x ∈ {9, 11}, y ∈ {-1, 1}, z ∈ {-1, 1}.
        assert_eq!(positions.len(), 8);
        for &(x, y, z) in &positions {
            assert!(x == 9000 || x == 11000, "Unexpected x: {x}");
            assert!(y == -1000 || y == 1000, "Unexpected y: {y}");
            assert!(z == -1000 || z == 1000, "Unexpected z: {z}");
        }

        // Original should still be at origin.
        let orig_shell = kernel.topo.solids.get(original).shell;
        let orig_face = kernel.topo.shells.get(orig_shell).faces[0];
        let orig_loop = kernel.topo.faces.get(orig_face).outer_loop;
        let orig_he = kernel.topo.loops.get(orig_loop).half_edge;
        let orig_vert = kernel.topo.half_edges.get(orig_he).origin;
        let orig_pid = kernel.topo.vertices.get(orig_vert).point_id;
        let p = kernel.geom.points[orig_pid as usize];
        assert!(p.x.abs() <= 1.0 + 1e-10, "Original was mutated: x={}", p.x);
    }

    #[test]
    fn test_tessellate_after_copy() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_box(2.0, 2.0, 2.0);
        let copy = kernel.copy_solid(solid);

        let shell = kernel.topo.solids.get(copy).shell;
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert_eq!(total_tris, 12, "Expected 12 triangles for copied box");
    }
}

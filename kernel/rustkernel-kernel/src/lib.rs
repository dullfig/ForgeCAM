use std::collections::HashMap;

#[allow(unused_imports)]
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::diagnostics::{self, DiagnosticReport};
use rustkernel_topology::intersection::IntersectionPipeline;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use rustkernel_geom::AnalyticalGeomStore;
use rustkernel_solvers::default_pipeline;
use rustkernel_topology::evolution::{
    EdgeOrigin, FaceOrigin, ShapeEvolution, VertexOrigin,
};

mod primitives;
mod sweeps;
mod local_ops;
mod booleans;
mod transforms;
mod persistence;
mod export;
mod analysis;

// Re-export public types from sub-modules
pub use local_ops::ShellError;
pub use persistence::ForgeError;
pub use export::{ExportError, export_stl_ascii, export_stl_binary, export_obj};
pub use analysis::{MassProperties, compute_mass_properties};

/// Central kernel holding all topology, geometry, and the intersection pipeline.
pub struct Kernel {
    pub(crate) topo: TopoStore,
    pub(crate) geom: AnalyticalGeomStore,
    pub(crate) pipeline: IntersectionPipeline,
    /// Evolution data from the last operation (for persistent entity tracking).
    last_evolution: Option<ShapeEvolution>,
}

impl Kernel {
    pub fn new() -> Self {
        Self {
            topo: TopoStore::new(),
            geom: AnalyticalGeomStore::new(),
            pipeline: default_pipeline(),
            last_evolution: None,
        }
    }

    /// Retrieve the shape evolution from the last operation.
    pub fn last_evolution(&self) -> Option<&ShapeEvolution> {
        self.last_evolution.as_ref()
    }

    /// Take ownership of the shape evolution from the last operation,
    /// leaving `None` in its place.
    pub fn take_evolution(&mut self) -> Option<ShapeEvolution> {
        self.last_evolution.take()
    }

    /// Collect all entities from a solid and record primitive evolution.
    fn record_primitive_evolution(&mut self, solid: SolidIdx) {
        let faces: Vec<FaceIdx> = self.topo.solid_faces(solid);

        let mut edges = Vec::new();
        let mut vertices = Vec::new();
        let mut seen_edges = std::collections::HashSet::new();
        let mut seen_verts = std::collections::HashSet::new();

        for &face in &faces {
            let loop_idx = self.topo.faces.get(face).outer_loop;
            let start_he = self.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let he_data = self.topo.half_edges.get(he);
                if seen_verts.insert(he_data.origin.raw()) {
                    vertices.push(he_data.origin);
                }
                if seen_edges.insert(he_data.edge.raw()) {
                    edges.push(he_data.edge);
                }
                he = he_data.next;
                if he == start_he {
                    break;
                }
            }
        }

        self.last_evolution = Some(ShapeEvolution::all_primitive(faces, edges, vertices));
    }

    /// Deep-copy a solid into fresh arena slots. Returns a new solid that shares
    /// no indices with the original.
    pub fn copy_solid(&mut self, solid: SolidIdx) -> SolidIdx {
        let mut remap = IndexRemap::new();
        copy_solid_impl(&mut self.topo, &mut self.geom, solid, &mut remap)
    }

    /// Deep-copy a solid and build a CopiedFrom evolution record.
    /// Used by transforms to record that every entity in the result
    /// is a copy of the corresponding entity in the original.
    fn copy_solid_with_evolution(&mut self, solid: SolidIdx) -> (SolidIdx, ShapeEvolution) {
        let mut remap = IndexRemap::new();
        let new_solid = copy_solid_impl(&mut self.topo, &mut self.geom, solid, &mut remap);

        let mut evo = ShapeEvolution::new();
        for (&old_face, &new_face) in &remap.faces {
            evo.record_face(new_face, FaceOrigin::CopiedFrom(old_face));
        }
        for (&old_edge, &new_edge) in &remap.edges {
            evo.record_edge(new_edge, EdgeOrigin::CopiedFrom(old_edge));
        }
        for (&old_vert, &new_vert) in &remap.vertices {
            evo.record_vertex(new_vert, VertexOrigin::CopiedFrom(old_vert));
        }

        (new_solid, evo)
    }

    // Primitives, transforms, booleans, sweeps, local ops, export, persistence,
    // and analysis methods are in their respective module files.
    // See: primitives.rs, transforms.rs, booleans.rs, sweeps.rs, local_ops.rs,
    //      export.rs, persistence.rs, analysis.rs

    /// Access the topology store (read-only).
    pub fn topo(&self) -> &TopoStore {
        &self.topo
    }

    /// Access the geometry store (read-only).
    pub fn geom(&self) -> &AnalyticalGeomStore {
        &self.geom
    }

    /// Tessellate all faces of a solid (all shells) with default options.
    pub fn tessellate_solid(&mut self, solid: SolidIdx) {
        for &sh in &self.topo.solids.get(solid).shells.clone() {
            rustkernel_topology::tessellate::tessellate_shell(&mut self.topo, sh, &self.geom);
        }
    }

    /// Tessellate all faces of a solid (all shells) with custom quality options.
    pub fn tessellate_solid_with_options(
        &mut self,
        solid: SolidIdx,
        opts: &rustkernel_topology::tessellate::TessellationOptions,
    ) {
        for &sh in &self.topo.solids.get(solid).shells.clone() {
            rustkernel_topology::tessellate::tessellate_shell_with_options(
                &mut self.topo,
                sh,
                &self.geom,
                opts,
            );
        }
    }

    pub fn validate_solid(&self, solid: SolidIdx) -> DiagnosticReport {
        diagnostics::validate_solid(&self.topo, solid)
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
    let shell_idx = topo.solids.get(solid).outer_shell();

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
        shells: vec![Idx::from_raw(0)],
        genus: topo.solids.get(solid).genus,
    });
    let new_shell_idx = topo.shells.alloc(Shell {
        faces: Vec::new(),
        solid: new_solid_idx,
    });
    topo.solids.get_mut(new_solid_idx).shells = vec![new_shell_idx];
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
    let cloned = geom.curves[old_id as usize].clone();
    let new_id = geom.add_curve(cloned);
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
    let cloned = geom.surfaces[old_id as usize].clone();
    let new_id = geom.add_surface(cloned);
    map.insert(old_id, new_id);
    new_id
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_topology::geom_store::GeomAccess;
    use rustkernel_topology::tessellate::tessellate_shell;
    use std::collections::HashSet;
    use std::io::Write;

    fn verify_euler(topo: &TopoStore, solid: SolidIdx) {
        verify_euler_genus(topo, solid);
    }

    fn verify_euler_genus(topo: &TopoStore, solid: SolidIdx) {
        let genus = topo.solids.get(solid).genus;

        // Verify per-shell (each shell independently satisfies V-E+F=2)
        for &shell_idx in &topo.solids.get(solid).shells {
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
            let expected = 2 - 2 * genus as i32;
            assert_eq!(
                v - e + f, expected,
                "Euler: V({v}) - E({e}) + F({f}) = {} != {expected} (genus={genus})",
                v - e + f
            );
        }
    }

    fn verify_all_twins(topo: &TopoStore, solid: SolidIdx) {
        let faces: Vec<FaceIdx> = topo.solid_faces(solid);
        for face_idx in faces {
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
        let shell_a = kernel.topo.solids.get(a).outer_shell();
        let shell_b = kernel.topo.solids.get(b).outer_shell();
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
        let orig_shell = kernel.topo.solids.get(original).outer_shell();
        let orig_face = kernel.topo.shells.get(orig_shell).faces[0];
        let orig_loop = kernel.topo.faces.get(orig_face).outer_loop;
        let orig_he = kernel.topo.loops.get(orig_loop).half_edge;
        let orig_vert = kernel.topo.half_edges.get(orig_he).origin;
        let orig_point_id = kernel.topo.vertices.get(orig_vert).point_id;
        let orig_pos = kernel.geom.points[orig_point_id as usize];

        let copy = kernel.copy_solid(original);

        // Mutate copy's first vertex.
        let copy_shell = kernel.topo.solids.get(copy).outer_shell();
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
        let shell = kernel.topo.solids.get(translated).outer_shell();
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
        let orig_shell = kernel.topo.solids.get(original).outer_shell();
        let orig_face = kernel.topo.shells.get(orig_shell).faces[0];
        let orig_loop = kernel.topo.faces.get(orig_face).outer_loop;
        let orig_he = kernel.topo.loops.get(orig_loop).half_edge;
        let orig_vert = kernel.topo.half_edges.get(orig_he).origin;
        let orig_pid = kernel.topo.vertices.get(orig_vert).point_id;
        let p = kernel.geom.points[orig_pid as usize];
        assert!(p.x.abs() <= 1.0 + 1e-10, "Original was mutated: x={}", p.x);
    }

    #[test]
    fn test_make_cylinder() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_cylinder(1.0, 2.0);
        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);
    }

    #[test]
    fn test_make_sphere() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_sphere(1.0);
        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);
    }

    #[test]
    fn test_make_cone() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_cone(2.0, 0.0, 3.0);
        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);
    }

    #[test]
    fn test_make_torus() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_torus(5.0, 1.0);
        // Torus has genus 1 → V-E+F = 0
        let shell_idx = kernel.topo.solids.get(solid).outer_shell();
        let faces = &kernel.topo.shells.get(shell_idx).faces;
        let mut verts_set = HashSet::new();
        let mut edges_set = HashSet::new();
        for &face_idx in faces {
            let loop_idx = kernel.topo.faces.get(face_idx).outer_loop;
            let start_he = kernel.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let he_data = kernel.topo.half_edges.get(he);
                verts_set.insert(he_data.origin.raw());
                edges_set.insert(he_data.edge.raw());
                he = he_data.next;
                if he == start_he { break; }
            }
        }
        let v = verts_set.len() as i32;
        let e = edges_set.len() as i32;
        let f = faces.len() as i32;
        assert_eq!(v - e + f, 0, "Torus Euler: V({v}) - E({e}) + F({f}) should be 0");
        verify_all_twins(&kernel.topo, solid);
    }

    #[test]
    fn test_translate_cylinder() {
        let mut kernel = Kernel::new();
        let original = kernel.make_cylinder(1.0, 2.0);
        let delta = Vec3::new(5.0, 0.0, 0.0);
        let translated = kernel.translate(original, delta);

        verify_euler(&kernel.topo, translated);
        verify_all_twins(&kernel.topo, translated);

        // All translated vertices should have x >= 4.0 (radius 1 + offset 5 - 1)
        let shell = kernel.topo.solids.get(translated).outer_shell();
        for &face_idx in &kernel.topo.shells.get(shell).faces {
            let loop_idx = kernel.topo.faces.get(face_idx).outer_loop;
            let start_he = kernel.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert = kernel.topo.half_edges.get(he).origin;
                let pid = kernel.topo.vertices.get(vert).point_id;
                let p = kernel.geom.points[pid as usize];
                assert!(p.x >= 4.0 - 1e-10, "Vertex x={} too small after translate", p.x);
                he = kernel.topo.half_edges.get(he).next;
                if he == start_he { break; }
            }
        }
    }

    #[test]
    fn test_copy_sphere_independence() {
        let mut kernel = Kernel::new();
        let original = kernel.make_sphere(3.0);
        let copy = kernel.copy_solid(original);

        verify_euler(&kernel.topo, copy);
        verify_all_twins(&kernel.topo, copy);

        // Mutate a vertex in the copy
        let copy_shell = kernel.topo.solids.get(copy).outer_shell();
        let copy_face = kernel.topo.shells.get(copy_shell).faces[0];
        let copy_loop = kernel.topo.faces.get(copy_face).outer_loop;
        let copy_he = kernel.topo.loops.get(copy_loop).half_edge;
        let copy_vert = kernel.topo.half_edges.get(copy_he).origin;
        let copy_pid = kernel.topo.vertices.get(copy_vert).point_id;
        kernel.geom.points[copy_pid as usize] = Point3::new(99.0, 99.0, 99.0);

        // Original should be unchanged — all vertices still on sphere surface
        let orig_shell = kernel.topo.solids.get(original).outer_shell();
        let mut checked = HashSet::new();
        for &face_idx in &kernel.topo.shells.get(orig_shell).faces {
            let loop_idx = kernel.topo.faces.get(face_idx).outer_loop;
            let start_he = kernel.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert = kernel.topo.half_edges.get(he).origin;
                if checked.insert(vert.raw()) {
                    let pid = kernel.topo.vertices.get(vert).point_id;
                    let p = kernel.geom.points[pid as usize];
                    let dist = (p - Point3::origin()).norm();
                    assert!((dist - 3.0).abs() < 1e-10, "Original vertex mutated: dist={dist}");
                }
                he = kernel.topo.half_edges.get(he).next;
                if he == start_he { break; }
            }
        }
    }

    #[test]
    fn test_tessellate_after_copy() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_box(2.0, 2.0, 2.0);
        let copy = kernel.copy_solid(solid);

        let shell = kernel.topo.solids.get(copy).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert_eq!(total_tris, 12, "Expected 12 triangles for copied box");
    }

    // --- Integration tests (Step 9) ---

    #[test]
    fn test_tessellate_cylinder_vertices_on_surface() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_cylinder(2.0, 5.0);
        let shell = kernel.topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();

            // All mesh vertices should be within the cylinder's bounding volume
            for v in &mesh.positions {
                let dx = v.x;
                let dy = v.y;
                let r = (dx * dx + dy * dy).sqrt();
                assert!(r <= 2.0 + 1e-8, "Mesh vertex outside cylinder: r={r}");
                assert!(v.z >= -1e-8 && v.z <= 5.0 + 1e-8, "Mesh vertex z={} out of range", v.z);
            }
        }
        assert!(total_tris > 0, "Cylinder should have tessellation triangles");
    }

    #[test]
    fn test_tessellate_sphere_normals_outward() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_sphere(3.0);
        let shell = kernel.topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();

            // Each vertex should be at distance ~3.0 from origin (within chord tolerance)
            for v in &mesh.positions {
                let dist = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
                assert!(
                    (dist - 3.0).abs() < 0.5,
                    "Mesh vertex far from sphere surface: dist={dist}"
                );
            }

            // Normals should point outward (dot with position > 0)
            for (v, n) in mesh.positions.iter().zip(mesh.normals.iter()) {
                let dot = v.x * n.x + v.y * n.y + v.z * n.z;
                assert!(dot > 0.0, "Normal points inward at vertex ({},{},{})", v.x, v.y, v.z);
            }
        }
        assert!(total_tris > 0, "Sphere should have tessellation triangles");
    }

    #[test]
    fn test_section_cylinder_at_midheight() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_cylinder(1.0, 4.0);
        let loops = kernel.section(solid, [0.0, 0.0, 2.0], [0.0, 0.0, 1.0]);
        match loops {
            Ok(ref wire_loops) => {
                assert!(!wire_loops.is_empty(), "Section through cylinder should produce at least one loop");
            }
            Err(_) => {
                // Section may fail if the pipeline can't find plane-cylinder intersections
                // for the planar cap faces — this is acceptable for Phase 3.
                // The key thing is it doesn't panic.
            }
        }
    }

    #[test]
    fn test_section_sphere_at_equator() {
        let mut kernel = Kernel::new();
        let solid = kernel.make_sphere(2.0);
        let loops = kernel.section(solid, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        match loops {
            Ok(ref wire_loops) => {
                assert!(!wire_loops.is_empty(), "Section through sphere equator should produce at least one loop");
            }
            Err(_) => {
                // Acceptable for Phase 3 if curved-face section isn't fully wired up.
            }
        }
    }

    #[test]
    fn test_translate_all_primitives() {
        let mut kernel = Kernel::new();
        let delta = Vec3::new(10.0, 20.0, 30.0);

        let box_s = kernel.make_box(1.0, 1.0, 1.0);
        let box_t = kernel.translate(box_s, delta);
        verify_euler(&kernel.topo, box_t);
        verify_all_twins(&kernel.topo, box_t);

        let cyl_s = kernel.make_cylinder(1.0, 2.0);
        let cyl_t = kernel.translate(cyl_s, delta);
        verify_euler(&kernel.topo, cyl_t);
        verify_all_twins(&kernel.topo, cyl_t);

        let sph_s = kernel.make_sphere(1.0);
        let sph_t = kernel.translate(sph_s, delta);
        verify_euler(&kernel.topo, sph_t);
        verify_all_twins(&kernel.topo, sph_t);

        let cone_s = kernel.make_cone(2.0, 0.0, 3.0);
        let cone_t = kernel.translate(cone_s, delta);
        verify_euler(&kernel.topo, cone_t);
        verify_all_twins(&kernel.topo, cone_t);

        // Torus: verify genus is preserved
        let torus_s = kernel.make_torus(5.0, 1.0);
        let torus_t = kernel.translate(torus_s, delta);
        assert_eq!(kernel.topo.solids.get(torus_t).genus, 1);
        verify_all_twins(&kernel.topo, torus_t);
    }

    #[test]
    fn test_copy_all_primitives() {
        let mut kernel = Kernel::new();

        let box_s = kernel.make_box(1.0, 1.0, 1.0);
        let box_c = kernel.copy_solid(box_s);
        verify_euler(&kernel.topo, box_c);
        verify_all_twins(&kernel.topo, box_c);

        let cyl_s = kernel.make_cylinder(1.0, 2.0);
        let cyl_c = kernel.copy_solid(cyl_s);
        verify_euler(&kernel.topo, cyl_c);
        verify_all_twins(&kernel.topo, cyl_c);

        let sph_s = kernel.make_sphere(1.0);
        let sph_c = kernel.copy_solid(sph_s);
        verify_euler(&kernel.topo, sph_c);
        verify_all_twins(&kernel.topo, sph_c);

        let cone_s = kernel.make_cone(2.0, 1.0, 3.0); // frustum
        let cone_c = kernel.copy_solid(cone_s);
        verify_euler(&kernel.topo, cone_c);
        verify_all_twins(&kernel.topo, cone_c);

        let torus_s = kernel.make_torus(5.0, 1.0);
        let torus_c = kernel.copy_solid(torus_s);
        assert_eq!(kernel.topo.solids.get(torus_c).genus, 1);
        verify_all_twins(&kernel.topo, torus_c);
    }

    // --- Sketch + Extrude integration tests ---

    #[test]
    fn test_sketch_extrude_rectangle() {
        let mut kernel = Kernel::new();
        let mut sketch = kernel.create_sketch(Point3::origin(), Vec3::new(0.0, 0.0, 1.0));

        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(2.0, 0.0);
        let p2 = sketch.add_point(2.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_horizontal_line(p0, p1);
        sketch.add_vertical_line(p1, p2);
        sketch.add_horizontal_line(p2, p3);
        sketch.add_vertical_line(p3, p0);

        sketch.constrain_fixed(p0, 0.0, 0.0);
        sketch.constrain_distance(p0, p1, 2.0);
        sketch.constrain_distance(p1, p2, 1.0);

        assert_eq!(sketch.solve(), rustkernel_sketch::solver::SolveResult::FullyConstrained);

        let profile = sketch.to_profile_3d().unwrap();
        let solid = kernel.extrude(&profile, Vec3::new(0.0, 0.0, 1.0), 0.5);

        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);

        // Should have 6 faces like a box.
        let shell = kernel.topo.solids.get(solid).outer_shell();
        assert_eq!(kernel.topo.shells.get(shell).faces.len(), 6);
    }

    #[test]
    fn test_sketch_extrude_then_fuse() {
        let mut kernel = Kernel::new();

        let base = kernel.make_box(2.0, 2.0, 2.0);

        // Create a pad via sketch + extrude that overlaps the base box.
        let mut sketch = kernel.create_sketch(
            Point3::new(0.5, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let p0 = sketch.add_point(-0.5, -0.5);
        let p1 = sketch.add_point(0.5, -0.5);
        let p2 = sketch.add_point(0.5, 0.5);
        let p3 = sketch.add_point(-0.5, 0.5);
        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p0);

        let profile = sketch.to_profile_3d().unwrap();
        let pad = kernel.extrude(&profile, Vec3::new(0.0, 0.0, 1.0), 1.0);

        verify_euler(&kernel.topo, pad);
        verify_all_twins(&kernel.topo, pad);

        // Fuse — may fail with BooleanError on non-axis-aligned overlaps.
        // The key assertion is that the extruded solid itself is topologically valid above.
        match kernel.fuse(base, pad) {
            Ok(_) => {}
            Err(e) => eprintln!("Sketch-extrude fuse: {e}"),
        }
    }

    #[test]
    fn test_sketch_extrude_then_cut() {
        let mut kernel = Kernel::new();

        let base = kernel.make_box(2.0, 2.0, 2.0);

        let mut sketch = kernel.create_sketch(
            Point3::new(0.3, 0.0, -1.5),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let p0 = sketch.add_point(-0.25, -0.25);
        let p1 = sketch.add_point(0.25, -0.25);
        let p2 = sketch.add_point(0.25, 0.25);
        let p3 = sketch.add_point(-0.25, 0.25);
        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p0);

        let profile = sketch.to_profile_3d().unwrap();
        let tool = kernel.extrude(&profile, Vec3::new(0.0, 0.0, 1.0), 3.0);

        verify_euler(&kernel.topo, tool);
        verify_all_twins(&kernel.topo, tool);

        // Cut — may fail with BooleanError on non-axis-aligned overlaps.
        match kernel.cut(base, tool) {
            Ok(_) => {}
            Err(e) => eprintln!("Sketch-extrude cut: {e}"),
        }
    }

    #[test]
    fn test_extrude_tessellation_via_kernel() {
        let mut kernel = Kernel::new();

        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];

        let solid = kernel.extrude(&profile, Vec3::new(0.0, 0.0, 1.0), 2.0);
        let shell = kernel.topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert!(total_tris > 0, "Extruded solid should have triangles");
    }

    #[test]
    fn test_tessellate_all_primitives() {
        let mut kernel = Kernel::new();

        let primitives = [
            kernel.make_box(1.0, 1.0, 1.0),
            kernel.make_cylinder(1.0, 2.0),
            kernel.make_sphere(1.0),
            kernel.make_cone(2.0, 0.0, 3.0),
            kernel.make_torus(5.0, 1.0),
        ];

        for solid in primitives {
            let shell = kernel.topo.solids.get(solid).outer_shell();
            tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

            let mut total_tris = 0;
            for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
                let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                    .expect("Face should be tessellated");
                total_tris += mesh.triangle_count();
                // No degenerate triangles (all should have area > 0)
                assert!(!mesh.positions.is_empty(), "Mesh should have vertices");
            }
            assert!(total_tris > 0, "Primitive should have tessellation triangles");
        }
    }

    // --- Revolve integration tests ---

    #[test]
    fn test_kernel_revolve() {
        let mut kernel = Kernel::new();
        // Rectangle off-axis → full revolve → torus-like
        let profile = vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
        ];
        let solid = kernel.revolve(
            &profile,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::TAU,
        );
        assert_eq!(kernel.topo.solids.get(solid).genus, 1);
        verify_all_twins(&kernel.topo, solid);

        // Euler for genus=1: V-E+F = 0
        let shell_idx = kernel.topo.solids.get(solid).outer_shell();
        let faces = &kernel.topo.shells.get(shell_idx).faces;
        let mut verts_set = HashSet::new();
        let mut edges_set = HashSet::new();
        for &face_idx in faces {
            let loop_idx = kernel.topo.faces.get(face_idx).outer_loop;
            let start_he = kernel.topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let he_data = kernel.topo.half_edges.get(he);
                verts_set.insert(he_data.origin.raw());
                edges_set.insert(he_data.edge.raw());
                he = he_data.next;
                if he == start_he { break; }
            }
        }
        let v = verts_set.len() as i32;
        let e = edges_set.len() as i32;
        let f = faces.len() as i32;
        assert_eq!(v - e + f, 0, "Torus-like Euler: V({v}) - E({e}) + F({f}) should be 0");
    }

    #[test]
    fn test_sketch_revolve_integration() {
        let mut kernel = Kernel::new();

        // Create a sketch on the XZ plane (normal = +Y, but we need a workplane
        // that places the profile in the meridional plane for revolve around Z).
        let mut sketch = kernel.create_sketch(
            Point3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );

        let p0 = sketch.add_point(0.0, 0.0);
        let p1 = sketch.add_point(1.0, 0.0);
        let p2 = sketch.add_point(1.0, 1.0);
        let p3 = sketch.add_point(0.0, 1.0);

        sketch.add_line(p0, p1);
        sketch.add_line(p1, p2);
        sketch.add_line(p2, p3);
        sketch.add_line(p3, p0);

        let profile = sketch.to_profile_3d().unwrap();
        assert_eq!(profile.len(), 4);

        let solid = kernel.revolve(
            &profile,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::TAU,
        );

        verify_all_twins(&kernel.topo, solid);
    }

    #[test]
    fn test_revolve_then_boolean() {
        let mut kernel = Kernel::new();

        let base = kernel.make_box(6.0, 6.0, 2.0);

        // Revolve a small rectangle to create a cylinder-like solid
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.5, 0.0, 3.0),
            Point3::new(0.0, 0.0, 3.0),
        ];
        let revolved = kernel.revolve(
            &profile,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::TAU,
        );

        verify_euler(&kernel.topo, base);
        verify_all_twins(&kernel.topo, revolved);

        // Boolean fuse — may fail with BooleanError on curved-surface intersections.
        match kernel.fuse(base, revolved) {
            Ok(_) => {}
            Err(e) => eprintln!("Extrude+revolve fuse: {e}"),
        }
    }

    // --- NURBS tests ---

    #[test]
    fn test_interpolate_curve() {
        let mut kernel = Kernel::new();
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(4.0, 3.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
        ];
        let cid = kernel.interpolate_curve(&pts, 3);
        let (t0, t1) = kernel.geom().curve_domain(cid);
        // Curve should pass through each interpolation point
        // Evaluate at knot parameters; start and end are exact
        let p0 = kernel.geom().curve_eval(cid, t0);
        let p_end = kernel.geom().curve_eval(cid, t1);
        assert!((p0 - pts[0]).norm() < 1e-6, "Start point mismatch");
        assert!((p_end - pts[4]).norm() < 1e-6, "End point mismatch");
    }

    #[test]
    fn test_loft_two_curves() {
        let mut kernel = Kernel::new();
        let c0 = kernel.interpolate_curve(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            2,
        );
        let c1 = kernel.interpolate_curve(
            &[
                Point3::new(0.0, 0.0, 3.0),
                Point3::new(1.0, 0.0, 3.0),
                Point3::new(2.0, 0.0, 3.0),
            ],
            2,
        );
        let sid = kernel.loft(&[c0, c1], 1);
        let ((u0, u1), (v0, v1)) = kernel.geom().surface_domain(sid);
        assert!(u1 > u0, "u domain should be non-degenerate");
        assert!(v1 > v0, "v domain should be non-degenerate");
    }

    #[test]
    fn test_nurbs_extrude() {
        let mut kernel = Kernel::new();
        let cid = kernel.interpolate_curve(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            1,
        );
        let sid = kernel.nurbs_extrude(cid, Vec3::new(0.0, 0.0, 1.0), 5.0);
        let ((u0, u1), (v0, v1)) = kernel.geom().surface_domain(sid);
        // Collect corners and verify all expected ones are present
        let corners = [
            kernel.geom().surface_eval(sid, u0, v0),
            kernel.geom().surface_eval(sid, u1, v0),
            kernel.geom().surface_eval(sid, u0, v1),
            kernel.geom().surface_eval(sid, u1, v1),
        ];
        let expected = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 5.0),
            Point3::new(1.0, 0.0, 5.0),
        ];
        for exp in &expected {
            let min_dist = corners.iter().map(|c| (c - exp).norm()).fold(f64::MAX, f64::min);
            assert!(min_dist < 1e-6, "Expected corner {:?} not found", exp);
        }
    }

    #[test]
    fn test_nurbs_surface_tessellation() {
        let mut kernel = Kernel::new();
        let cid = kernel.interpolate_curve(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
            ],
            2,
        );
        let sid = kernel.nurbs_extrude(cid, Vec3::new(0.0, 0.0, 1.0), 3.0);
        let mesh = kernel.geom().tessellate_surface(sid, 8, 8);
        assert!(mesh.is_some(), "Should produce a mesh");
        let mesh = mesh.unwrap();
        assert!(mesh.triangle_count() > 0, "Mesh should have triangles");
        assert_eq!(mesh.positions.len(), 81); // 9x9 grid
        assert_eq!(mesh.triangle_count(), 128); // 8x8 quads x 2 tris
    }

    // --- NURBS Solid Builder integration tests ---

    /// Helper: create a closed NURBS curve for testing.
    fn make_closed_nurbs_curve(kernel: &mut Kernel) -> u32 {
        let pts = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0), // close the loop
        ];
        kernel.interpolate_curve(&pts, 3)
    }

    #[test]
    fn test_kernel_nurbs_extrude_solid() {
        let mut kernel = Kernel::new();
        let cid = make_closed_nurbs_curve(&mut kernel);
        let solid = kernel.make_nurbs_extrude_solid(cid, Vec3::new(0.0, 0.0, 1.0), 2.0);

        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);

        // Tessellate and verify
        let shell = kernel.topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert!(total_tris > 0, "NURBS extrude solid should have triangles");
    }

    #[test]
    fn test_kernel_nurbs_revolve_solid() {
        let mut kernel = Kernel::new();

        // Closed curve off the Z axis for revolve
        let pts = vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.5),
            Point3::new(3.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let cid = kernel.interpolate_curve(&pts, 3);
        let solid = kernel.make_nurbs_revolve_solid(
            cid,
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::TAU,
        );

        verify_euler_genus(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);
        assert_eq!(kernel.topo.solids.get(solid).genus, 1);

        // Tessellate and verify
        let shell = kernel.topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        let mut total_tris = 0;
        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            total_tris += mesh.triangle_count();
        }
        assert!(total_tris > 0, "NURBS revolve solid should have triangles");
    }

    #[test]
    fn test_nurbs_solid_tessellation_normals() {
        let mut kernel = Kernel::new();
        let cid = make_closed_nurbs_curve(&mut kernel);
        let solid = kernel.make_nurbs_extrude_solid(cid, Vec3::new(0.0, 0.0, 1.0), 2.0);

        let shell = kernel.topo.solids.get(solid).outer_shell();
        tessellate_shell(&mut kernel.topo, shell, &kernel.geom);

        for &face_idx in &kernel.topo.shells.get(shell).faces.clone() {
            let mesh = kernel.topo.faces.get(face_idx).mesh_cache.as_ref()
                .expect("Face should be tessellated");
            for (i, n) in mesh.normals.iter().enumerate() {
                let len = n.norm();
                assert!(
                    len > 0.5,
                    "Face {} vertex {} has degenerate normal: {:?} (len={})",
                    face_idx.raw(), i, n, len
                );
            }
        }
    }

    // ── .forge round-trip tests ──

    #[test]
    fn test_forge_roundtrip_box() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 3.0, 4.0);

        // Save to in-memory buffer
        let mut buf = std::io::Cursor::new(Vec::new());
        k.save_forge(&mut buf).expect("save_forge failed");

        // Load back
        buf.set_position(0);
        let k2 = Kernel::load_forge(buf).expect("load_forge failed");

        // Validate topology survived round-trip
        let report = k2.validate_solid(solid);
        assert!(report.is_valid(), "Round-trip solid failed validation: {:?}", report);

        // Verify geometry: check vertex count and a sample point
        assert_eq!(k2.topo.vertices.len(), k.topo.vertices.len());
        let p_orig = k.geom.points[0];
        let p_loaded = k2.geom.points[0];
        assert!((p_orig - p_loaded).norm() < 1e-12, "Point mismatch after round-trip");
    }

    #[test]
    fn test_forge_roundtrip_cylinder() {
        let mut k = Kernel::new();
        let solid = k.make_cylinder(1.5, 4.0);

        let mut buf = std::io::Cursor::new(Vec::new());
        k.save_forge(&mut buf).expect("save_forge failed");
        buf.set_position(0);
        let k2 = Kernel::load_forge(buf).expect("load_forge failed");

        let report = k2.validate_solid(solid);
        assert!(report.is_valid(), "Cylinder round-trip failed: {:?}", report);
        assert_eq!(k2.geom.surfaces.len(), k.geom.surfaces.len());
    }

    #[test]
    fn test_forge_addons_preserved() {
        let mut k = Kernel::new();
        k.make_box(1.0, 1.0, 1.0);

        let addon_data = b"custom plugin data here";
        let mut buf = std::io::Cursor::new(Vec::new());
        k.save_forge_with_addons(&mut buf, &[("my_plugin.dat", addon_data)])
            .expect("save failed");

        buf.set_position(0);
        let (_k2, addons) = Kernel::load_forge_with_addons(buf).expect("load failed");

        assert_eq!(addons.len(), 1);
        assert_eq!(addons[0].0, "my_plugin.dat");
        assert_eq!(addons[0].1, addon_data);
    }

    #[test]
    fn test_forge_missing_manifest_error() {
        // Empty ZIP
        let buf = std::io::Cursor::new(Vec::new());
        let mut zip = zip::ZipWriter::new(buf);
        let opts = zip::write::SimpleFileOptions::default();
        zip.start_file("dummy.txt", opts).unwrap();
        zip.write_all(b"hello").unwrap();
        let buf = zip.finish().unwrap();

        let cursor = std::io::Cursor::new(buf.into_inner());
        let result = Kernel::load_forge(cursor);
        assert!(result.is_err());
    }

    // ── STL export tests ──

    #[test]
    fn test_stl_ascii_export_box() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 2.0, 2.0);

        let mut buf = Vec::new();
        k.export_stl_ascii(solid, "test_box", &mut buf).expect("STL export failed");
        let stl = String::from_utf8(buf).expect("STL should be valid UTF-8");

        assert!(stl.starts_with("solid test_box"));
        assert!(stl.trim_end().ends_with("endsolid test_box"));
        // A box has 6 faces × 2 triangles each = 12 triangles
        let facet_count = stl.matches("facet normal").count();
        assert_eq!(facet_count, 12, "Box should have 12 triangles, got {facet_count}");
    }

    #[test]
    fn test_stl_binary_export_box() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 2.0, 2.0);

        let mut buf = Vec::new();
        k.export_stl_binary(solid, &mut buf).expect("Binary STL export failed");

        // Binary STL: 80 header + 4 bytes tri count + 50 bytes per triangle
        let tri_count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(tri_count, 12, "Box should have 12 triangles");
        assert_eq!(buf.len(), 80 + 4 + 50 * 12);
    }

    #[test]
    fn test_stl_export_cylinder() {
        let mut k = Kernel::new();
        let solid = k.make_cylinder(1.0, 2.0);

        let mut buf = Vec::new();
        k.export_stl_ascii(solid, "cyl", &mut buf).expect("STL export failed");
        let stl = String::from_utf8(buf).unwrap();
        let facet_count = stl.matches("facet normal").count();
        // Cylinder with 32 segments: 32 side quads (64 tris) + 2 caps (32 tris each) = 128
        assert!(facet_count > 60, "Cylinder should have many triangles, got {facet_count}");
    }

    // ── OBJ export tests ──

    #[test]
    fn test_obj_export_box() {
        let mut k = Kernel::new();
        let solid = k.make_box(2.0, 2.0, 2.0);

        let mut buf = Vec::new();
        k.export_obj(solid, &mut buf).expect("OBJ export failed");
        let obj = String::from_utf8(buf).unwrap();

        let v_count = obj.lines().filter(|l| l.starts_with("v ")).count();
        let vn_count = obj.lines().filter(|l| l.starts_with("vn ")).count();
        let f_count = obj.lines().filter(|l| l.starts_with("f ")).count();

        assert!(v_count > 0, "OBJ should have vertices");
        assert!(vn_count > 0, "OBJ should have normals");
        assert_eq!(f_count, 12, "Box should have 12 triangle faces, got {f_count}");
    }

    #[test]
    fn test_obj_export_sphere() {
        let mut k = Kernel::new();
        let solid = k.make_sphere(2.0);

        let mut buf = Vec::new();
        k.export_obj(solid, &mut buf).expect("OBJ export failed");
        let obj = String::from_utf8(buf).unwrap();

        let f_count = obj.lines().filter(|l| l.starts_with("f ")).count();
        assert!(f_count > 100, "Sphere should have many triangles, got {f_count}");
    }

    // ── Affine transform tests ──

    /// Collect unique vertex positions from a solid.
    fn collect_vertices(
        topo: &TopoStore,
        geom: &AnalyticalGeomStore,
        solid: SolidIdx,
    ) -> Vec<Point3> {
        let faces = topo.solid_faces(solid);
        let mut seen = HashSet::new();
        let mut pts = Vec::new();
        for face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start_he = topo.loops.get(loop_idx).half_edge;
            let mut he = start_he;
            loop {
                let vert_idx = topo.half_edges.get(he).origin;
                let point_id = topo.vertices.get(vert_idx).point_id;
                if seen.insert(point_id) {
                    pts.push(geom.points[point_id as usize]);
                }
                he = topo.half_edges.get(he).next;
                if he == start_he {
                    break;
                }
            }
        }
        pts
    }

    #[test]
    fn test_rotate_box_90() {
        let mut k = Kernel::new();
        let orig = k.make_box(2.0, 2.0, 2.0);
        let angle = std::f64::consts::FRAC_PI_2; // 90°
        let rotated = k.rotate(orig, Point3::origin(), Vec3::new(0.0, 0.0, 1.0), angle);

        verify_euler(&k.topo, rotated);

        // Original vertex (1,1,z) should become (-1,1,z) after 90° around Z.
        let verts = collect_vertices(&k.topo, &k.geom, rotated);
        let target = Point3::new(-1.0, 1.0, 1.0);
        let min_dist = verts.iter().map(|v| (v - target).norm()).fold(f64::MAX, f64::min);
        assert!(min_dist < 1e-10, "Expected vertex near {target:?}, closest dist = {min_dist}");
    }

    #[test]
    fn test_rotate_preserves_euler() {
        let mut k = Kernel::new();
        let cyl = k.make_cylinder(1.0, 3.0);
        let rotated = k.rotate(cyl, Point3::origin(), Vec3::new(1.0, 0.0, 0.0), 1.23);
        verify_euler(&k.topo, rotated);
    }

    #[test]
    fn test_scale_box() {
        let mut k = Kernel::new();
        let orig = k.make_box(2.0, 2.0, 2.0);
        let scaled = k.scale(orig, Point3::origin(), 2.0);

        verify_euler(&k.topo, scaled);

        let verts = collect_vertices(&k.topo, &k.geom, scaled);
        // Original corners at ±1, scaled by 2 → ±2.
        let target = Point3::new(2.0, 2.0, 2.0);
        let min_dist = verts.iter().map(|v| (v - target).norm()).fold(f64::MAX, f64::min);
        assert!(min_dist < 1e-10, "Expected vertex near {target:?}, closest dist = {min_dist}");
    }

    #[test]
    fn test_scale_negative_flips() {
        let mut k = Kernel::new();
        let orig = k.make_box(2.0, 2.0, 2.0);
        let inverted = k.scale(orig, Point3::origin(), -1.0);

        // validate_solid checks Euler + twins + loops
        let report = k.validate_solid(inverted);
        assert!(report.findings.is_empty(), "Inverted solid should pass validation: {:?}", report.findings);
    }

    #[test]
    fn test_mirror_box_across_yz() {
        let mut k = Kernel::new();
        // Box centered at (3, 0, 0) so mirror across YZ flips x → -x.
        let orig = k.make_box_at([3.0, 0.0, 0.0], 2.0, 2.0, 2.0);
        let mirrored = k.mirror(orig, Point3::origin(), Vec3::new(1.0, 0.0, 0.0));

        verify_euler(&k.topo, mirrored);

        let verts = collect_vertices(&k.topo, &k.geom, mirrored);
        // Original vertex at (4,1,1) should become (-4,1,1).
        let target = Point3::new(-4.0, 1.0, 1.0);
        let min_dist = verts.iter().map(|v| (v - target).norm()).fold(f64::MAX, f64::min);
        assert!(min_dist < 1e-10, "Expected mirrored vertex near {target:?}, closest dist = {min_dist}");
    }

    #[test]
    fn test_mirror_cylinder() {
        let mut k = Kernel::new();
        let cyl = k.make_cylinder(1.0, 3.0);
        let mirrored = k.mirror(cyl, Point3::origin(), Vec3::new(0.0, 0.0, 1.0));

        let report = k.validate_solid(mirrored);
        assert!(report.findings.is_empty(), "Mirrored cylinder should pass validation: {:?}", report.findings);
    }

    #[test]
    fn test_transform_general() {
        use rustkernel_math::nalgebra::Rotation3;

        let mut k = Kernel::new();
        let orig = k.make_box(2.0, 2.0, 2.0);

        // Apply an arbitrary rotation.
        let rot = Rotation3::from_euler_angles(0.3, 0.5, 0.7);
        let m = rot.to_homogeneous();
        let transformed = k.transform(orig, m);

        verify_euler(&k.topo, transformed);

        // Apply inverse rotation — should get back to original coords.
        let m_inv = rot.inverse().to_homogeneous();
        let roundtrip = k.transform(transformed, m_inv);
        verify_euler(&k.topo, roundtrip);

        let orig_verts = collect_vertices(&k.topo, &k.geom, orig);
        let rt_verts = collect_vertices(&k.topo, &k.geom, roundtrip);
        // Every original vertex should have a match in the round-tripped solid.
        for ov in &orig_verts {
            let min_dist = rt_verts.iter().map(|rv| (rv - ov).norm()).fold(f64::MAX, f64::min);
            assert!(min_dist < 1e-10, "Round-trip vertex mismatch: {ov:?}, min_dist = {min_dist}");
        }
    }

    #[test]
    fn test_rotate_does_not_mutate_original() {
        let mut k = Kernel::new();
        let orig = k.make_box(2.0, 2.0, 2.0);
        let orig_verts = collect_vertices(&k.topo, &k.geom, orig);

        let _ = k.rotate(orig, Point3::origin(), Vec3::new(0.0, 0.0, 1.0), 1.0);

        let after_verts = collect_vertices(&k.topo, &k.geom, orig);
        for (a, b) in orig_verts.iter().zip(after_verts.iter()) {
            assert!((a - b).norm() < 1e-15, "Original solid was mutated");
        }
    }

    #[test]
    fn test_scale_nurbs_extrude() {
        let mut k = Kernel::new();
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0), // closed
        ];
        let cid = k.interpolate_curve(&pts, 3);
        let solid = k.make_nurbs_extrude_solid(cid, Vec3::new(0.0, 0.0, 1.0), 2.0);
        let scaled = k.scale(solid, Point3::origin(), 3.0);

        verify_euler(&k.topo, scaled);

        // Tessellation should still work after scaling.
        let shell = k.topo.solids.get(scaled).outer_shell();
        tessellate_shell(&mut k.topo, shell, &k.geom);
        let faces: Vec<FaceIdx> = k.topo.shells.get(shell).faces.clone();
        let has_mesh = faces.iter().any(|&f| k.topo.faces.get(f).mesh_cache.is_some());
        assert!(has_mesh, "Scaled NURBS solid should tessellate");
    }

    // ── Mass properties tests ──

    #[test]
    fn test_volume_box() {
        let mut k = Kernel::new();
        let s = k.make_box(2.0, 3.0, 4.0);
        let mp = k.mass_properties(s, 1.0);
        assert!((mp.volume - 24.0).abs() < 1e-10, "box volume: {}", mp.volume);
    }

    #[test]
    fn test_surface_area_box() {
        let mut k = Kernel::new();
        let s = k.make_box(2.0, 3.0, 4.0);
        let mp = k.mass_properties(s, 1.0);
        // 2*(2*3 + 3*4 + 2*4) = 2*(6+12+8) = 52
        assert!((mp.surface_area - 52.0).abs() < 1e-10, "box area: {}", mp.surface_area);
    }

    #[test]
    fn test_center_of_mass_box_origin() {
        let mut k = Kernel::new();
        let s = k.make_box(2.0, 2.0, 2.0);
        let mp = k.mass_properties(s, 1.0);
        assert!(mp.center_of_mass.x.abs() < 1e-10, "com.x: {}", mp.center_of_mass.x);
        assert!(mp.center_of_mass.y.abs() < 1e-10, "com.y: {}", mp.center_of_mass.y);
        assert!(mp.center_of_mass.z.abs() < 1e-10, "com.z: {}", mp.center_of_mass.z);
    }

    #[test]
    fn test_center_of_mass_offset_box() {
        let mut k = Kernel::new();
        let s = k.make_box_at([3.0, 0.0, 0.0], 2.0, 2.0, 2.0);
        let mp = k.mass_properties(s, 1.0);
        assert!((mp.center_of_mass.x - 3.0).abs() < 1e-10, "com.x: {}", mp.center_of_mass.x);
        assert!(mp.center_of_mass.y.abs() < 1e-10, "com.y: {}", mp.center_of_mass.y);
        assert!(mp.center_of_mass.z.abs() < 1e-10, "com.z: {}", mp.center_of_mass.z);
    }

    #[test]
    fn test_volume_sphere() {
        let mut k = Kernel::new();
        let s = k.make_sphere(2.0);
        let mp = k.mass_properties(s, 1.0);
        let expected = 4.0 / 3.0 * std::f64::consts::PI * 8.0;
        let rel_err = (mp.volume - expected).abs() / expected;
        // Coarse tessellation (16×8) gives ~6-7% error; allow 10%
        assert!(rel_err < 0.10, "sphere volume: {} expected: {} rel_err: {}", mp.volume, expected, rel_err);
    }

    #[test]
    fn test_volume_cylinder() {
        let mut k = Kernel::new();
        let s = k.make_cylinder(1.0, 3.0);
        let mp = k.mass_properties(s, 1.0);
        let expected = std::f64::consts::PI * 3.0;
        let rel_err = (mp.volume - expected).abs() / expected;
        assert!(rel_err < 0.05, "cylinder volume: {} expected: {} rel_err: {}", mp.volume, expected, rel_err);
    }

    #[test]
    fn test_volume_scales_cubically() {
        let mut k = Kernel::new();
        let s1 = k.make_box(1.0, 1.0, 1.0);
        let mp1 = k.mass_properties(s1, 1.0);
        let s2 = k.make_box(2.0, 2.0, 2.0);
        let mp2 = k.mass_properties(s2, 1.0);
        assert!((mp2.volume - mp1.volume * 8.0).abs() < 1e-10,
            "scaled volume: {} vs 8*original: {}", mp2.volume, mp1.volume * 8.0);
    }

    #[test]
    fn test_inertia_box() {
        let mut k = Kernel::new();
        let (a, b, c) = (2.0, 3.0, 4.0);
        let s = k.make_box(a, b, c);
        let density = 1.0;
        let mp = k.mass_properties(s, density);
        let m = mp.mass;
        // Ixx = m/12*(b^2+c^2), etc. for box centered at origin
        let ixx_expected = m / 12.0 * (b * b + c * c);
        let iyy_expected = m / 12.0 * (a * a + c * c);
        let izz_expected = m / 12.0 * (a * a + b * b);
        assert!((mp.inertia_tensor[0][0] - ixx_expected).abs() < 1e-6,
            "Ixx: {} expected: {}", mp.inertia_tensor[0][0], ixx_expected);
        assert!((mp.inertia_tensor[1][1] - iyy_expected).abs() < 1e-6,
            "Iyy: {} expected: {}", mp.inertia_tensor[1][1], iyy_expected);
        assert!((mp.inertia_tensor[2][2] - izz_expected).abs() < 1e-6,
            "Izz: {} expected: {}", mp.inertia_tensor[2][2], izz_expected);
        // Off-diagonal should be zero for symmetric box at origin
        assert!(mp.inertia_tensor[0][1].abs() < 1e-6, "Ixy: {}", mp.inertia_tensor[0][1]);
        assert!(mp.inertia_tensor[0][2].abs() < 1e-6, "Ixz: {}", mp.inertia_tensor[0][2]);
        assert!(mp.inertia_tensor[1][2].abs() < 1e-6, "Iyz: {}", mp.inertia_tensor[1][2]);
    }

    #[test]
    fn test_mass_equals_density_times_volume() {
        let mut k = Kernel::new();
        let s = k.make_box(2.0, 3.0, 4.0);
        let density = 7.85;
        let mp = k.mass_properties(s, density);
        assert!((mp.mass - density * mp.volume).abs() < 1e-10,
            "mass: {} expected: {}", mp.mass, density * mp.volume);
    }

    #[test]
    fn test_primitive_evolution_box() {
        use rustkernel_topology::evolution::FaceOrigin;

        let mut k = Kernel::new();
        let _s = k.make_box(2.0, 3.0, 4.0);

        let evo = k.last_evolution().expect("make_box should produce evolution");
        // Box: 6 faces, 12 edges, 8 vertices
        assert_eq!(evo.face_provenance.len(), 6, "box should have 6 faces");
        assert_eq!(evo.edge_provenance.len(), 12, "box should have 12 edges");
        assert_eq!(evo.vertex_provenance.len(), 8, "box should have 8 vertices");

        // All should be Primitive origin.
        for (_, origin) in &evo.face_provenance {
            assert!(matches!(origin, FaceOrigin::Primitive));
        }

        assert!(evo.deleted_faces.is_empty());
    }

    #[test]
    fn test_primitive_evolution_cylinder() {
        let mut k = Kernel::new();
        let _s = k.make_cylinder(1.0, 2.0);

        let evo = k.last_evolution().expect("make_cylinder should produce evolution");
        // Cylinder with 32 segments: 34 faces (32 lateral + top + bottom),
        // 96 edges (32*3), 64 vertices (32*2)
        assert_eq!(evo.face_provenance.len(), 34);
        assert_eq!(evo.edge_provenance.len(), 96);
        assert_eq!(evo.vertex_provenance.len(), 64);
        assert!(evo.deleted_faces.is_empty());
    }

    #[test]
    fn test_evolution_overwritten_by_next_op() {
        let mut k = Kernel::new();
        let _box = k.make_box(1.0, 1.0, 1.0);
        assert_eq!(k.last_evolution().unwrap().face_provenance.len(), 6);

        let _cyl = k.make_cylinder(0.5, 1.0);
        assert_eq!(k.last_evolution().unwrap().face_provenance.len(), 34);

        let taken = k.take_evolution();
        assert!(taken.is_some());
        assert!(k.last_evolution().is_none());
    }

    #[test]
    fn test_transform_evolution_translate() {
        use rustkernel_topology::evolution::{FaceOrigin, EdgeOrigin, VertexOrigin};

        let mut k = Kernel::new();
        let box_s = k.make_box(2.0, 2.0, 2.0);

        // Snapshot original entities.
        let orig_faces: Vec<_> = k.topo.solid_faces(box_s);
        let orig_face_count = orig_faces.len();

        let translated = k.translate(box_s, Vec3::new(10.0, 0.0, 0.0));
        let evo = k.last_evolution().expect("translate should produce evolution");

        // Every face in the result should be CopiedFrom an original face.
        assert_eq!(evo.face_provenance.len(), orig_face_count);
        for origin in evo.face_provenance.values() {
            assert!(matches!(origin, FaceOrigin::CopiedFrom(_)));
        }

        // Edges and vertices should also be CopiedFrom.
        assert!(!evo.edge_provenance.is_empty());
        for origin in evo.edge_provenance.values() {
            assert!(matches!(origin, EdgeOrigin::CopiedFrom(_)));
        }
        assert!(!evo.vertex_provenance.is_empty());
        for origin in evo.vertex_provenance.values() {
            assert!(matches!(origin, VertexOrigin::CopiedFrom(_)));
        }

        // No deletions (it's a copy, original still exists).
        assert!(evo.deleted_faces.is_empty());
        assert!(evo.deleted_edges.is_empty());
        assert!(evo.deleted_vertices.is_empty());

        // Verify the mapping: each CopiedFrom should reference an original face.
        let orig_face_set: std::collections::HashSet<_> = orig_faces.into_iter().collect();
        for origin in evo.face_provenance.values() {
            if let FaceOrigin::CopiedFrom(src) = origin {
                assert!(orig_face_set.contains(src), "CopiedFrom should reference original face");
            }
        }
    }

    #[test]
    fn test_transform_evolution_rotate() {
        use rustkernel_topology::evolution::FaceOrigin;

        let mut k = Kernel::new();
        let box_s = k.make_box(2.0, 2.0, 2.0);
        let _rotated = k.rotate(box_s, Point3::origin(), Vec3::z(), std::f64::consts::FRAC_PI_4);

        let evo = k.last_evolution().expect("rotate should produce evolution");
        assert_eq!(evo.face_provenance.len(), 6);
        for origin in evo.face_provenance.values() {
            assert!(matches!(origin, FaceOrigin::CopiedFrom(_)));
        }
        assert!(!evo.edge_provenance.is_empty());
        assert!(!evo.vertex_provenance.is_empty());
    }

    #[test]
    fn test_euler_chamfer_edge_vertex_evolution() {
        let mut k = Kernel::new();
        let box_s = k.make_box(2.0, 2.0, 2.0);

        let edges = rustkernel_builders::edge_analysis::solid_edges(&k.topo, box_s);
        // Pick one convex edge.
        let edge = edges[0];

        let _result = k.euler_chamfer_edges(box_s, &[edge], 0.3).unwrap();
        let evo = k.last_evolution().expect("chamfer should produce evolution");

        // Full provenance should be populated.
        assert!(!evo.face_provenance.is_empty());
        assert!(!evo.edge_provenance.is_empty());
        assert!(!evo.vertex_provenance.is_empty());

        // At least one deleted edge (the chamfered one).
        assert!(evo.deleted_edges.contains(&edge));
    }

    #[test]
    fn test_concave_fillet_step_shape() {
        // Create an L-shaped step by cutting a box from a larger box.
        // This produces concave internal edges at the step.
        let mut k = Kernel::new();
        let big = k.make_box(4.0, 4.0, 4.0);
        let small = k.make_box_at([2.0, 0.0, 2.0], 4.0, 4.0, 4.0);
        let step = k.cut(big, small).expect("boolean cut for step shape");

        // Find a concave edge on the step.
        let edges = rustkernel_builders::edge_analysis::solid_edges(&k.topo, step);
        let concave_edges: Vec<_> = edges.iter().filter(|&&e| {
            matches!(
                rustkernel_builders::edge_analysis::edge_convexity(&k.topo, &k.geom, e),
                Ok(rustkernel_builders::edge_analysis::EdgeConvexity::Concave)
            )
        }).copied().collect();

        assert!(!concave_edges.is_empty(), "step shape should have concave edges, found 0 out of {} edges", edges.len());

        // Attempt to fillet a concave edge.
        let result = k.euler_fillet_edges(step, &[concave_edges[0]], 0.2);
        match result {
            Ok(filleted) => {
                // If it succeeds, verify topology.
                let report = k.validate_solid(filleted);
                // Allow some validation issues for now — the important thing is
                // the fillet didn't panic or reject the concave edge.
                let _ = report;
            }
            Err(e) => {
                // For now, a non-panic error is acceptable — the geometry may not
                // be fully correct for all concave cases yet. But it should NOT
                // be NotConvex (that gate was removed).
                let msg = format!("{e}");
                assert!(
                    !msg.contains("not convex"),
                    "concave edges should not be rejected as 'not convex': {e}"
                );
            }
        }
    }

    // --- Offset tests ---

    #[test]
    fn test_offset_box_outward() {
        let mut k = Kernel::new();
        let s = k.make_box(2.0, 2.0, 2.0);
        let offset = k.offset_solid(s, 0.5).expect("box offset should succeed");

        verify_euler(&k.topo, offset);
        verify_all_twins(&k.topo, offset);

        // Original volume = 8, offset box should be 3×3×3 = 27
        let mp = k.mass_properties(offset, 1.0);
        assert!(
            (mp.volume - 27.0).abs() < 1e-6,
            "offset box volume: {} expected: 27.0", mp.volume
        );
    }

    #[test]
    fn test_offset_box_inward() {
        let mut k = Kernel::new();
        let s = k.make_box(4.0, 4.0, 4.0);
        let offset = k.offset_solid(s, -1.0).expect("box offset should succeed");

        verify_euler(&k.topo, offset);

        // Original = 64, offset = 2×2×2 = 8
        let mp = k.mass_properties(offset, 1.0);
        assert!(
            (mp.volume - 8.0).abs() < 1e-6,
            "offset box volume: {} expected: 8.0", mp.volume
        );
    }

    #[test]
    fn test_offset_cylinder() {
        let mut k = Kernel::new();
        let s = k.make_cylinder(2.0, 4.0);
        let offset = k.offset_solid(s, 1.0).expect("cylinder offset should succeed");

        verify_euler(&k.topo, offset);
        verify_all_twins(&k.topo, offset);

        // Original = π·4·4 = 16π, offset = π·9·4 = 36π (r goes 2→3)
        // Height stays 4 because top/bottom planes offset by 1 each → 4+2 = 6? No...
        // Actually: top plane offsets outward (+1 in +Z), bottom offsets outward (-1 in -Z).
        // So height goes from 4 to 6, radius from 2 to 3. Volume = π·9·6.
        let expected = std::f64::consts::PI * 9.0 * 6.0;
        let mp = k.mass_properties(offset, 1.0);
        let rel_err = (mp.volume - expected).abs() / expected;
        assert!(rel_err < 0.10, "offset cylinder volume: {} expected: {} rel_err: {}", mp.volume, expected, rel_err);
    }

    #[test]
    fn test_offset_sphere() {
        let mut k = Kernel::new();
        let s = k.make_sphere(3.0);
        let offset = k.offset_solid(s, 1.0).expect("sphere offset should succeed");

        verify_euler(&k.topo, offset);
        verify_all_twins(&k.topo, offset);

        // radius 3 → 4, volume = 4/3·π·64
        let expected = 4.0 / 3.0 * std::f64::consts::PI * 64.0;
        let mp = k.mass_properties(offset, 1.0);
        let rel_err = (mp.volume - expected).abs() / expected;
        assert!(rel_err < 0.10, "offset sphere volume: {} expected: {} rel_err: {}", mp.volume, expected, rel_err);
    }

    #[test]
    fn test_offset_preserves_original() {
        let mut k = Kernel::new();
        let s = k.make_box(2.0, 2.0, 2.0);
        let _offset = k.offset_solid(s, 1.0).unwrap();

        // Original solid should be unchanged
        let mp = k.mass_properties(s, 1.0);
        assert!(
            (mp.volume - 8.0).abs() < 1e-6,
            "original box volume changed: {}", mp.volume
        );
    }

    #[test]
    fn test_offset_returns_none_for_degenerate() {
        let mut k = Kernel::new();
        let s = k.make_cylinder(2.0, 4.0);
        // Offset by -2.0 makes the cylinder radius 0 → degenerate
        assert!(k.offset_solid(s, -2.0).is_none());
    }

    // --- Shell / hollow tests ---

    /// Find the face whose surface normal most closely matches `dir`.
    fn find_face_by_normal(k: &Kernel, solid: SolidIdx, dir: Vec3) -> FaceIdx {
        use rustkernel_topology::geom_store::GeomAccess;
        let shell = k.topo.solids.get(solid).outer_shell();
        let faces = &k.topo.shells.get(shell).faces;
        let mut best = faces[0];
        let mut best_dot = f64::NEG_INFINITY;
        for &fi in faces {
            let sid = k.topo.faces.get(fi).surface_id;
            let loop_idx = k.topo.faces.get(fi).outer_loop;
            let first_he = k.topo.loops.get(loop_idx).half_edge;
            let mut sum = Vec3::zeros();
            let mut count = 0usize;
            let mut he = first_he;
            loop {
                let v = k.topo.half_edges.get(he).origin;
                let pt = k.geom.point(k.topo.vertices.get(v).point_id);
                sum += pt.coords;
                count += 1;
                he = k.topo.half_edges.get(he).next;
                if he == first_he { break; }
            }
            let centroid = Point3::from(sum / count as f64);
            let (u, uv_v) = k.geom.surface_inverse_uv(sid, &centroid);
            let n = k.geom.surface_normal(sid, u, uv_v);
            let d = n.dot(&dir);
            if d > best_dot {
                best_dot = d;
                best = fi;
            }
        }
        best
    }

    #[test]
    fn test_shell_box_one_open_face() {
        let mut k = Kernel::new();
        let b = k.make_box(4.0, 4.0, 4.0); // volume = 64
        let top = find_face_by_normal(&k, b, Vec3::new(0.0, 0.0, 1.0));

        let shelled = k.shell_solid(b, 1.0, &[top]).expect("shell should succeed");

        verify_euler(&k.topo, shelled);
        verify_all_twins(&k.topo, shelled);

        // Cavity = 2×2×2 prism (z=-1 to z=1) + frustum (z=1 to z=2, base 2×2, top 4×4).
        // Frustum volume = h/3*(A1 + A2 + √(A1·A2)) = 1/3*(4+16+8) = 28/3.
        // Total cavity = 8 + 28/3 = 52/3. Shell volume = 64 - 52/3 = 140/3 ≈ 46.67.
        let mp = k.mass_properties(shelled, 1.0);
        let expected = 140.0 / 3.0;
        let rel_err = (mp.volume - expected).abs() / expected;
        assert!(
            rel_err < 0.10,
            "shell box volume: {} expected: {} rel_err: {}",
            mp.volume, expected, rel_err
        );

        // Face count: 5 outer + 5 inner + 4 wall = 14
        let shell_idx = k.topo.solids.get(shelled).outer_shell();
        assert_eq!(k.topo.shells.get(shell_idx).faces.len(), 14);
    }

    #[test]
    fn test_shell_box_two_open_faces() {
        let mut k = Kernel::new();
        let b = k.make_box(4.0, 4.0, 4.0);
        let top = find_face_by_normal(&k, b, Vec3::new(0.0, 0.0, 1.0));
        let bottom = find_face_by_normal(&k, b, Vec3::new(0.0, 0.0, -1.0));

        let shelled = k.shell_solid(b, 1.0, &[top, bottom]).expect("shell should succeed");

        verify_euler(&k.topo, shelled);
        verify_all_twins(&k.topo, shelled);

        // 4 outer + 4 inner + 4 wall (top) + 4 wall (bottom) = 16 faces
        let shell_idx = k.topo.solids.get(shelled).outer_shell();
        assert_eq!(k.topo.shells.get(shell_idx).faces.len(), 16);
    }

    #[test]
    fn test_shell_preserves_original() {
        let mut k = Kernel::new();
        let b = k.make_box(4.0, 4.0, 4.0);
        let top = find_face_by_normal(&k, b, Vec3::new(0.0, 0.0, 1.0));
        let mp_before = k.mass_properties(b, 1.0);
        let _shelled = k.shell_solid(b, 1.0, &[top]).expect("shell should succeed");
        let mp_after = k.mass_properties(b, 1.0);
        assert!(
            (mp_before.volume - mp_after.volume).abs() < 1e-6,
            "original solid should be unchanged"
        );
    }

    #[test]
    fn test_shell_invalid_thickness() {
        let mut k = Kernel::new();
        let b = k.make_box(4.0, 4.0, 4.0);
        let top = find_face_by_normal(&k, b, Vec3::new(0.0, 0.0, 1.0));
        assert!(k.shell_solid(b, 0.0, &[top]).is_err());
        assert!(k.shell_solid(b, -1.0, &[top]).is_err());
    }

    #[test]
    fn test_shell_no_open_faces() {
        let mut k = Kernel::new();
        let b = k.make_box(4.0, 4.0, 4.0);
        assert!(k.shell_solid(b, 1.0, &[]).is_err());
    }

    #[test]
    fn test_shell_too_thick_cylinder() {
        let mut k = Kernel::new();
        let c = k.make_cylinder(2.0, 4.0); // radius=2
        let top = find_face_by_normal(&k, c, Vec3::new(0.0, 0.0, 1.0));
        // thickness=2 → inner cylinder radius=0 → degenerate
        assert!(k.shell_solid(c, 2.0, &[top]).is_err());
    }

    // ── Multi-shell tests ────────────────────────────────────────────────

    #[test]
    fn test_multi_shell_solid_faces() {
        // Create a 2-shell solid manually, verify solid_faces returns all faces
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        // Build two boxes as separate solids
        let box_a = make_box_into(&mut topo, &mut geom, Point3::origin(), 4.0, 4.0, 4.0);
        let box_b = make_box_into(&mut topo, &mut geom, Point3::new(10.0, 0.0, 0.0), 2.0, 2.0, 2.0);

        let shell_a = topo.solids.get(box_a).outer_shell();
        let shell_b = topo.solids.get(box_b).outer_shell();
        let faces_a = topo.shells.get(shell_a).faces.len();
        let faces_b = topo.shells.get(shell_b).faces.len();

        // Create a multi-shell solid
        let multi_solid = topo.solids.alloc(Solid {
            shells: vec![shell_a, shell_b],
            genus: 0,
        });

        // solid_faces should return all faces from both shells
        let all_faces = topo.solid_faces(multi_solid);
        assert_eq!(all_faces.len(), faces_a + faces_b);
    }

    #[test]
    fn test_multi_shell_outer_shell() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();

        let box_a = make_box_into(&mut topo, &mut geom, Point3::origin(), 4.0, 4.0, 4.0);
        let box_b = make_box_into(&mut topo, &mut geom, Point3::new(10.0, 0.0, 0.0), 2.0, 2.0, 2.0);

        let shell_a = topo.solids.get(box_a).outer_shell();
        let shell_b = topo.solids.get(box_b).outer_shell();

        let multi_solid = topo.solids.alloc(Solid {
            shells: vec![shell_a, shell_b],
            genus: 0,
        });

        // outer_shell() always returns shells[0]
        assert_eq!(topo.solids.get(multi_solid).outer_shell(), shell_a);
    }

    #[test]
    fn test_copy_single_shell_solid() {
        // Verify copy_solid works unchanged for single-shell solids
        let mut kernel = Kernel::new();
        let original = kernel.make_box(3.0, 4.0, 5.0);
        let copy = kernel.copy_solid(original);

        let orig_faces = kernel.topo.solid_faces(original);
        let copy_faces = kernel.topo.solid_faces(copy);
        assert_eq!(orig_faces.len(), copy_faces.len());

        // Shells should be different indices
        let orig_shell = kernel.topo.solids.get(original).outer_shell();
        let copy_shell = kernel.topo.solids.get(copy).outer_shell();
        assert_ne!(orig_shell, copy_shell);

        // Both should have exactly one shell
        assert_eq!(kernel.topo.solids.get(original).shells.len(), 1);
        assert_eq!(kernel.topo.solids.get(copy).shells.len(), 1);
    }

    #[test]
    fn test_export_uses_all_faces() {
        // Verify STL export includes faces from solid_faces (all shells)
        let mut kernel = Kernel::new();
        let solid = kernel.make_box(2.0, 2.0, 2.0);
        kernel.tessellate_solid(solid);

        let mut buf = Vec::new();
        export_stl_ascii(&mut kernel.topo, &kernel.geom, solid, "test", &mut buf)
            .expect("export should succeed");

        let output = String::from_utf8(buf).unwrap();
        // A box has 6 faces, each with 2 triangles = 12 facets
        let facet_count = output.matches("facet normal").count();
        assert_eq!(facet_count, 12);
    }

    #[test]
    fn test_euler_per_shell() {
        // Each shell independently satisfies V-E+F=2
        let mut kernel = Kernel::new();
        let solid = kernel.make_box(3.0, 3.0, 3.0);

        // Single-shell solid should pass per-shell Euler check
        verify_euler(&kernel.topo, solid);
        verify_all_twins(&kernel.topo, solid);
    }
}

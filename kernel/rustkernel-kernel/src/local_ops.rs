use std::collections::HashMap;
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::arena::Idx;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use rustkernel_geom::{AnalyticalGeomStore, CurveDef};
use rustkernel_builders::chamfer_builder::{self, chamfer_edges_into};
use rustkernel_builders::fillet_builder::{self, fillet_edges_into};
use rustkernel_builders::euler_chamfer;
use rustkernel_builders::euler_fillet;
use rustkernel_builders::cylinder_builder::{build_face_from_vert_idxs, match_twins_from_map};
use tracing::info_span;

use crate::Kernel;

/// Error type for the shell/hollow operation.
#[derive(Debug)]
pub enum ShellError {
    /// Wall thickness must be positive.
    InvalidThickness,
    /// At least one face must be open.
    NoOpenFaces,
    /// An open face does not belong to the specified solid.
    FaceNotInSolid,
    /// Surface offset failed (e.g., NURBS not supported, or thickness collapses geometry).
    OffsetFailed,
}

impl std::fmt::Display for ShellError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ShellError::InvalidThickness => write!(f, "wall thickness must be positive"),
            ShellError::NoOpenFaces => write!(f, "at least one face must be open"),
            ShellError::FaceNotInSolid => write!(f, "open face does not belong to solid"),
            ShellError::OffsetFailed => write!(f, "surface offset failed"),
        }
    }
}

impl std::error::Error for ShellError {}

/// Offset all geometry of a solid in-place: surfaces, vertex positions, and edge curves.
fn offset_solid_geometry(
    topo: &mut TopoStore,
    geom: &mut AnalyticalGeomStore,
    solid: SolidIdx,
    distance: f64,
) {
    use std::collections::HashSet;
    use rustkernel_geom::LineSegment;
    use rustkernel_topology::geom_store::GeomAccess;

    let faces: Vec<FaceIdx> = topo.solid_faces(solid);

    // Phase 1: Offset surfaces.
    let mut visited_surfaces: HashSet<u32> = HashSet::new();
    for &face_idx in &faces {
        let surface_id = topo.faces.get(face_idx).surface_id;
        if visited_surfaces.insert(surface_id) {
            // Pre-check guarantees this succeeds.
            let offset = geom.surfaces[surface_id as usize]
                .offset(distance)
                .expect("pre-checked");
            geom.surfaces[surface_id as usize] = offset;
        }
    }

    // Phase 2: Move vertices along the combined normal from all adjacent faces.
    // A corner vertex shared by N faces must move along the intersection of
    // the N offset surfaces.  For analytical surfaces the exact formula is:
    //   displacement = d * N_faces * n_sum / |n_sum|²
    // where n_sum = Σ nᵢ (surface normals at the vertex from each adjacent face).
    // This correctly handles box corners (3 orthogonal normals), cylinder rims
    // (1 cap + 1 lateral), and single-surface vertices (spheres).

    // Pass 1: collect normal sum and face count per vertex.
    let mut vert_normal_sum: HashMap<u32, (Vec3, usize)> = HashMap::new();
    for &face_idx in &faces {
        let surface_id = topo.faces.get(face_idx).surface_id;
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let vert_idx = topo.half_edges.get(he).origin;
            let point_id = topo.vertices.get(vert_idx).point_id;
            let pos = geom.points[point_id as usize];
            let (u, v) = geom.surface_inverse_uv(surface_id, &pos);
            let normal = geom.surface_normal(surface_id, u, v);

            let entry = vert_normal_sum.entry(point_id).or_insert((Vec3::zeros(), 0));
            entry.0 += normal;
            entry.1 += 1;

            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }

    // Pass 2: displace each vertex.
    for (&point_id, &(n_sum, k)) in &vert_normal_sum {
        let mag_sq = n_sum.norm_squared();
        if mag_sq < 1e-24 {
            continue; // Degenerate — normals cancel out
        }
        let displacement = distance * (k as f64) * n_sum / mag_sq;
        geom.points[point_id as usize] += displacement;
    }

    // Phase 3: Recompute edge curves from updated vertex positions.
    let mut visited_edges: HashSet<u32> = HashSet::new();
    for &face_idx in &faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;
        let mut he = start_he;
        loop {
            let edge_idx = topo.half_edges.get(he).edge;
            if visited_edges.insert(edge_idx.raw()) {
                let origin = topo.half_edges.get(he).origin;
                let next_he = topo.half_edges.get(he).next;
                let dest = topo.half_edges.get(next_he).origin;

                let start_pt = geom.points[topo.vertices.get(origin).point_id as usize];
                let end_pt = geom.points[topo.vertices.get(dest).point_id as usize];

                let curve_id = topo.edges.get(edge_idx).curve_id;
                geom.curves[curve_id as usize] =
                    CurveDef::LineSegment(LineSegment { start: start_pt, end: end_pt });
            }

            he = topo.half_edges.get(he).next;
            if he == start_he {
                break;
            }
        }
    }
}

impl Kernel {
    // --- Offset ---

    /// Deep-copy a solid and offset all faces by `distance` along their outward normals.
    ///
    /// Positive distance = outward (enlarge), negative = inward (shrink).
    /// Returns `None` if any face's surface cannot be offset (e.g., NURBS).
    ///
    /// This produces a new solid whose surfaces are parallel to the original at
    /// distance `d`.  Vertex positions are moved along the local surface normal.
    /// Edge curves are recomputed as line segments between new vertex positions.
    ///
    /// Note: this is a *surface-level* offset — it does not recompute edge
    /// intersections, so results are exact for prismatic solids and approximate
    /// for solids with high curvature near edges.  The shell/hollow operation
    /// (future) will use this + booleans for a fully robust result.
    pub fn offset_solid(&mut self, solid: SolidIdx, distance: f64) -> Option<SolidIdx> {
        let _span = info_span!("kernel.offset_solid", solid = solid.raw(), distance).entered();

        // Pre-check: verify all surfaces can be offset before copying anything.
        {
            let faces: Vec<FaceIdx> = self.topo.solid_faces(solid);
            let mut checked: std::collections::HashSet<u32> = std::collections::HashSet::new();
            for &face_idx in &faces {
                let sid = self.topo.faces.get(face_idx).surface_id;
                if checked.insert(sid) {
                    if self.geom.surfaces[sid as usize].offset(distance).is_none() {
                        return None;
                    }
                }
            }
        }

        let new_solid = self.copy_solid(solid);
        offset_solid_geometry(&mut self.topo, &mut self.geom, new_solid, distance);
        Some(new_solid)
    }

    /// Create a hollow (shelled) solid by removing `open_faces` and adding
    /// uniform-thickness walls on all remaining faces.
    ///
    /// The result is a new solid. The original is not modified.
    ///
    /// Algorithm: for each kept face, emit an outer copy (original geometry) and
    /// an inner copy (offset inward by `thickness`, normals flipped). At each
    /// opening, wall faces connect the outer boundary to the inner boundary.
    pub fn shell_solid(
        &mut self,
        solid: SolidIdx,
        thickness: f64,
        open_faces: &[FaceIdx],
    ) -> Result<SolidIdx, ShellError> {
        let _span = info_span!("kernel.shell_solid", solid = solid.raw(), thickness).entered();

        if thickness <= 0.0 {
            return Err(ShellError::InvalidThickness);
        }
        if open_faces.is_empty() {
            return Err(ShellError::NoOpenFaces);
        }

        let all_faces: Vec<FaceIdx> = self.topo.solid_faces(solid);
        let open_set: std::collections::HashSet<FaceIdx> =
            open_faces.iter().cloned().collect();

        // Validate open faces belong to this solid.
        for &f in open_faces {
            if !all_faces.contains(&f) {
                return Err(ShellError::FaceNotInSolid);
            }
        }

        // Pre-check: all kept surfaces can be offset.
        let mut checked: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for &fi in &all_faces {
            if open_set.contains(&fi) {
                continue;
            }
            let sid = self.topo.faces.get(fi).surface_id;
            if checked.insert(sid) {
                if self.geom.surfaces[sid as usize].offset(-thickness).is_none() {
                    return Err(ShellError::OffsetFailed);
                }
            }
        }

        let kept_faces: Vec<FaceIdx> = all_faces
            .iter()
            .filter(|f| !open_set.contains(f))
            .cloned()
            .collect();

        // --- Collect all unique vertices ---
        let mut all_verts: Vec<VertexIdx> = Vec::new();
        for &fi in &all_faces {
            let loop_idx = self.topo.faces.get(fi).outer_loop;
            let first_he = self.topo.loops.get(loop_idx).half_edge;
            let mut he = first_he;
            loop {
                let v = self.topo.half_edges.get(he).origin;
                if !all_verts.contains(&v) {
                    all_verts.push(v);
                }
                he = self.topo.half_edges.get(he).next;
                if he == first_he {
                    break;
                }
            }
        }

        // --- Compute inner vertex positions (multi-face normal formula) ---
        use rustkernel_topology::geom_store::GeomAccess;

        let mut vert_normals: HashMap<VertexIdx, Vec<Vec3>> = HashMap::new();
        for &fi in &all_faces {
            let sid = self.topo.faces.get(fi).surface_id;
            let loop_idx = self.topo.faces.get(fi).outer_loop;
            let first_he = self.topo.loops.get(loop_idx).half_edge;
            let mut he = first_he;
            loop {
                let v = self.topo.half_edges.get(he).origin;
                let pt = self.geom.point(self.topo.vertices.get(v).point_id);
                let (u, uv_v) = self.geom.surface_inverse_uv(sid, &pt);
                let n = self.geom.surface_normal(sid, u, uv_v);
                vert_normals.entry(v).or_default().push(n);
                he = self.topo.half_edges.get(he).next;
                if he == first_he {
                    break;
                }
            }
        }

        let mut inner_positions: HashMap<VertexIdx, Point3> = HashMap::new();
        for &v in &all_verts {
            let pt = self.geom.point(self.topo.vertices.get(v).point_id);
            let normals = &vert_normals[&v];
            let n_sum: Vec3 = normals.iter().copied().fold(Vec3::zeros(), |a, b| a + b);
            let k = normals.len() as f64;
            let denom = n_sum.norm_squared();
            if denom < 1e-12 {
                inner_positions.insert(v, pt);
            } else {
                let displacement = (-thickness) * k * n_sum / denom;
                inner_positions.insert(v, pt + displacement);
            }
        }

        // --- Collect boundary edge loops (one loop per open face) ---
        struct BoundaryEdge {
            v0: VertexIdx,
            v1: VertexIdx,
        }
        let mut boundary_loops: Vec<Vec<BoundaryEdge>> = Vec::new();
        for &open_fi in open_faces {
            let loop_idx = self.topo.faces.get(open_fi).outer_loop;
            let first_he = self.topo.loops.get(loop_idx).half_edge;
            let mut he = first_he;
            let mut edges = Vec::new();
            loop {
                let v0 = self.topo.half_edges.get(he).origin;
                let next_he = self.topo.half_edges.get(he).next;
                let v1 = self.topo.half_edges.get(next_he).origin;
                edges.push(BoundaryEdge { v0, v1 });
                he = next_he;
                if he == first_he {
                    break;
                }
            }
            boundary_loops.push(edges);
        }

        // --- Build new solid ---
        let new_shell = self.topo.shells.alloc(Shell {
            faces: Vec::new(),
            solid: Idx::from_raw(0),
        });
        // Genus: shelling with N open faces creates a surface with genus max(0, N-1).
        // 1 opening → bowl (genus 0), 2 openings → tube (genus 1), etc.
        let genus = if open_faces.len() > 1 {
            (open_faces.len() - 1) as u32
        } else {
            0
        };
        let new_solid = self.topo.solids.alloc(Solid { shells: vec![new_shell], genus });
        self.topo.shells.get_mut(new_shell).solid = new_solid;

        // Create outer vertices (copies of originals).
        let mut outer_vmap: HashMap<VertexIdx, VertexIdx> = HashMap::new();
        for &v in &all_verts {
            let pt = self.geom.point(self.topo.vertices.get(v).point_id);
            let pid = self.geom.add_point(pt);
            let nv = self.topo.vertices.alloc(Vertex { point_id: pid });
            outer_vmap.insert(v, nv);
        }

        // Create inner vertices (offset positions).
        let mut inner_vmap: HashMap<VertexIdx, VertexIdx> = HashMap::new();
        for &v in &all_verts {
            let pt = inner_positions[&v];
            let pid = self.geom.add_point(pt);
            let nv = self.topo.vertices.alloc(Vertex { point_id: pid });
            inner_vmap.insert(v, nv);
        }

        let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

        // Outer faces (kept faces, original geometry).
        for &fi in &kept_faces {
            let sid = self.topo.faces.get(fi).surface_id;
            let new_surf = self.geom.surfaces[sid as usize].clone();
            let new_sid = self.geom.add_surface(new_surf);

            let verts = self.face_vertex_loop(fi);
            let mapped: Vec<VertexIdx> = verts.iter().map(|v| outer_vmap[v]).collect();

            build_face_from_vert_idxs(
                &mut self.topo,
                &mut self.geom,
                &mapped,
                new_sid,
                new_shell,
                &mut he_map,
            );
        }

        // Inner faces (kept faces, offset + flipped, reversed winding).
        for &fi in &kept_faces {
            let sid = self.topo.faces.get(fi).surface_id;
            let mut new_surf = self.geom.surfaces[sid as usize]
                .offset(-thickness)
                .unwrap(); // pre-checked
            new_surf.flip_normal();
            let new_sid = self.geom.add_surface(new_surf);

            let mut verts = self.face_vertex_loop(fi);
            verts.reverse(); // flip winding for inner face
            let mapped: Vec<VertexIdx> = verts.iter().map(|v| inner_vmap[v]).collect();

            build_face_from_vert_idxs(
                &mut self.topo,
                &mut self.geom,
                &mapped,
                new_sid,
                new_shell,
                &mut he_map,
            );
        }

        // Wall faces at each opening.
        for bloop in &boundary_loops {
            for be in bloop {
                let ov0 = outer_vmap[&be.v0];
                let ov1 = outer_vmap[&be.v1];
                let iv1 = inner_vmap[&be.v1];
                let iv0 = inner_vmap[&be.v0];

                let p0 = self.geom.point(self.topo.vertices.get(ov0).point_id);
                let p1 = self.geom.point(self.topo.vertices.get(ov1).point_id);
                let p3 = self.geom.point(self.topo.vertices.get(iv0).point_id);
                let e1 = p1 - p0;
                let e2 = p3 - p0;
                let normal = e1.cross(&e2).normalize();
                let p2 = self.geom.point(self.topo.vertices.get(iv1).point_id);
                let origin = Point3::from(
                    (p0.coords + p1.coords + p2.coords + p3.coords) * 0.25,
                );

                let plane = rustkernel_geom::Plane { origin, normal };
                let wall_sid =
                    self.geom.add_surface(rustkernel_geom::SurfaceDef::Plane(plane));

                build_face_from_vert_idxs(
                    &mut self.topo,
                    &mut self.geom,
                    &[ov0, ov1, iv1, iv0],
                    wall_sid,
                    new_shell,
                    &mut he_map,
                );
            }
        }

        // Twin matching.
        match_twins_from_map(&mut self.topo, &he_map);

        Ok(new_solid)
    }

    /// Walk a face's outer loop and return the vertex indices in order.
    fn face_vertex_loop(&self, face: FaceIdx) -> Vec<VertexIdx> {
        let loop_idx = self.topo.faces.get(face).outer_loop;
        let first_he = self.topo.loops.get(loop_idx).half_edge;
        let mut result = Vec::new();
        let mut he = first_he;
        loop {
            result.push(self.topo.half_edges.get(he).origin);
            he = self.topo.half_edges.get(he).next;
            if he == first_he {
                break;
            }
        }
        result
    }

    // --- Chamfer & Fillet ---

    /// Apply a constant-distance chamfer to straight edges between planar faces.
    /// All edges must be convex and no two edges may share a vertex.
    pub fn chamfer_edges(
        &mut self,
        solid: SolidIdx,
        edges: &[EdgeIdx],
        distance: f64,
    ) -> Result<SolidIdx, chamfer_builder::ChamferError> {
        let _span = info_span!("kernel.chamfer_edges", solid = solid.raw(), n_edges = edges.len(), distance).entered();
        let (result, evo) = chamfer_edges_into(&mut self.topo, &mut self.geom, solid, edges, distance)?;
        self.last_evolution = Some(evo);
        Ok(result)
    }

    /// Apply a constant-distance Euler-based chamfer (in-place mutation).
    /// All edges must be convex and no two edges may share a vertex.
    pub fn euler_chamfer_edges(
        &mut self,
        solid: SolidIdx,
        edges: &[EdgeIdx],
        distance: f64,
    ) -> Result<SolidIdx, euler_chamfer::EulerChamferError> {
        let _span = info_span!("kernel.euler_chamfer_edges", solid = solid.raw(), n_edges = edges.len(), distance).entered();
        let (result, evo) = euler_chamfer::euler_chamfer_edges(&mut self.topo, &mut self.geom, solid, edges, distance)?;
        self.last_evolution = Some(evo);
        Ok(result)
    }

    /// Apply a constant-radius Euler-based fillet (in-place mutation).
    /// Supports edges meeting at corner vertices. Uses 8 arc segments by default.
    pub fn euler_fillet_edges(
        &mut self,
        solid: SolidIdx,
        edges: &[EdgeIdx],
        radius: f64,
    ) -> Result<SolidIdx, euler_fillet::EulerFilletError> {
        let _span = info_span!("kernel.euler_fillet_edges", solid = solid.raw(), n_edges = edges.len(), radius).entered();
        let (result, evo) = euler_fillet::euler_fillet_edges(&mut self.topo, &mut self.geom, solid, edges, radius, 8)?;
        self.last_evolution = Some(evo);
        Ok(result)
    }

    /// Apply per-edge-radius Euler-based fillets (in-place mutation).
    /// `radii[i]` is the fillet radius for `edges[i]`. Dissimilar radii at a corner
    /// vertex produce a NURBS blend patch. Uses 8 arc segments by default.
    pub fn euler_fillet_edges_with_radii(
        &mut self,
        solid: SolidIdx,
        edges: &[EdgeIdx],
        radii: &[f64],
    ) -> Result<SolidIdx, euler_fillet::EulerFilletError> {
        let _span = info_span!("kernel.euler_fillet_edges_with_radii", solid = solid.raw(), n_edges = edges.len()).entered();
        let (result, evo) = euler_fillet::euler_fillet_edges_with_radii(&mut self.topo, &mut self.geom, solid, edges, radii, 8)?;
        self.last_evolution = Some(evo);
        Ok(result)
    }

    /// Apply a constant-radius fillet to straight edges between planar faces (rebuild-based fallback).
    /// All edges must be convex and no two edges may share a vertex.
    /// Uses 8 arc segments by default.
    pub fn fillet_edges(
        &mut self,
        solid: SolidIdx,
        edges: &[EdgeIdx],
        radius: f64,
    ) -> Result<SolidIdx, fillet_builder::FilletError> {
        let _span = info_span!("kernel.fillet_edges", solid = solid.raw(), n_edges = edges.len(), radius).entered();
        let (result, evo) = fillet_edges_into(&mut self.topo, &mut self.geom, solid, edges, radius, 8)?;
        self.last_evolution = Some(evo);
        Ok(result)
    }
}

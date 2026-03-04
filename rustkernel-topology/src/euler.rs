//! Euler operators for B-Rep topology mutation.
//!
//! Four atomic topology operators that preserve V - E + F = 2 - 2g by construction:
//! - **SEMV**: Split Edge, Make Vertex (V+1, E+1)
//! - **MEF**: Make Edge, Face (E+1, F+1)
//! - **KEV**: Kill Edge, Vertex — inverse of SEMV (V-1, E-1)
//! - **KEF**: Kill Edge, Face — inverse of MEF (E-1, F-1)
//!
//! Arena is append-only: kill operators leave dead elements as unreachable orphans.
//! `validate_solid()` walks from Solid→Shell→faces so dead slots are invisible.

use crate::store::TopoStore;
use crate::topo::*;
use tracing::debug_span;

// ─── Error type ──────────────────────────────────────────────────────────────

#[derive(Debug)]
pub enum EulerError {
    /// h1 and h2 are the same half-edge (MEF requires distinct).
    SameHalfEdge,
    /// h1 and h2 are not in the same loop.
    DifferentLoops,
    /// The edge does not have the specified vertex as an endpoint.
    VertexNotOnEdge(EdgeIdx, VertexIdx),
    /// The specified face is not adjacent to the edge (KEF).
    FaceNotAdjacentToEdge(EdgeIdx, FaceIdx),
    /// Half-edge has no twin (needed for operation).
    NoTwin(HalfEdgeIdx),
}

impl std::fmt::Display for EulerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EulerError::SameHalfEdge => write!(f, "h1 and h2 must be distinct"),
            EulerError::DifferentLoops => write!(f, "h1 and h2 must be in the same loop"),
            EulerError::VertexNotOnEdge(e, v) => {
                write!(f, "vertex {:?} not on edge {:?}", v, e)
            }
            EulerError::FaceNotAdjacentToEdge(e, face) => {
                write!(f, "face {:?} not adjacent to edge {:?}", face, e)
            }
            EulerError::NoTwin(he) => write!(f, "half-edge {:?} has no twin", he),
        }
    }
}

impl std::error::Error for EulerError {}

// ─── Result types ────────────────────────────────────────────────────────────

pub struct SemvResult {
    pub new_vertex: VertexIdx,
    pub new_edge: EdgeIdx,
    /// New half-edge in h's loop (V_new → V_b, on edge E_new).
    pub new_he_a: HalfEdgeIdx,
    /// New half-edge in twin's loop (V_new → V_a, on edge E), if twin exists.
    pub new_he_b: Option<HalfEdgeIdx>,
}

pub struct MefResult {
    pub new_edge: EdgeIdx,
    pub new_face: FaceIdx,
    pub new_loop: LoopIdx,
    /// Half-edge in old loop (V2 → V1).
    pub he_a: HalfEdgeIdx,
    /// Half-edge in new loop (V1 → V2).
    pub he_b: HalfEdgeIdx,
}

pub struct KevResult {
    pub killed_vertex: VertexIdx,
    pub killed_edge: EdgeIdx,
}

pub struct KefResult {
    pub killed_face: FaceIdx,
    pub killed_edge: EdgeIdx,
}

// ─── Utilities ───────────────────────────────────────────────────────────────

/// Find the half-edge that precedes `target_he` in the given loop.
/// O(n) walk; loop sizes are typically 3-6 edges.
///
/// Panics if `target_he` is not found in the loop.
pub fn find_prev_he(topo: &TopoStore, loop_idx: LoopIdx, target_he: HalfEdgeIdx) -> HalfEdgeIdx {
    let start = topo.loops.get(loop_idx).half_edge;
    let mut he = start;
    loop {
        let next = topo.half_edges.get(he).next;
        if next == target_he {
            return he;
        }
        he = next;
        if he == start {
            panic!(
                "target_he Idx({}) not found in loop Idx({})",
                target_he.raw(),
                loop_idx.raw()
            );
        }
    }
}

/// Collect all half-edge indices in a loop into a Vec.
pub fn collect_loop_hes(topo: &TopoStore, loop_idx: LoopIdx) -> Vec<HalfEdgeIdx> {
    let start = topo.loops.get(loop_idx).half_edge;
    let mut result = Vec::new();
    let mut he = start;
    loop {
        result.push(he);
        he = topo.half_edges.get(he).next;
        if he == start {
            break;
        }
    }
    result
}

/// Count the number of half-edges (edges) in a loop.
pub fn loop_edge_count(topo: &TopoStore, loop_idx: LoopIdx) -> usize {
    let start = topo.loops.get(loop_idx).half_edge;
    let mut count = 0;
    let mut he = start;
    loop {
        count += 1;
        he = topo.half_edges.get(he).next;
        if he == start {
            break;
        }
    }
    count
}

/// Find the half-edge originating from `vertex` in the outer loop of `face`.
/// Returns `None` if no such half-edge exists.
pub fn find_he_from_vertex(
    topo: &TopoStore,
    face: FaceIdx,
    vertex: VertexIdx,
) -> Option<HalfEdgeIdx> {
    let loop_idx = topo.faces.get(face).outer_loop;
    let start = topo.loops.get(loop_idx).half_edge;
    let mut he = start;
    loop {
        if topo.half_edges.get(he).origin == vertex {
            return Some(he);
        }
        he = topo.half_edges.get(he).next;
        if he == start {
            return None;
        }
    }
}

// ─── SEMV: Split Edge, Make Vertex ───────────────────────────────────────────

/// Split the edge at half-edge `h`, inserting a new vertex.
///
/// `h` is a half-edge from V_a to V_b. After the operation:
/// - V_new is inserted between V_a and V_b
/// - Edge E (h's edge) becomes V_a -- V_new
/// - New edge E_new represents V_new -- V_b
///
/// Euler: V+1, E+1, F unchanged.
pub fn semv(
    topo: &mut TopoStore,
    h: HalfEdgeIdx,
    point_id: u32,
    curve_a_id: u32,
    curve_b_id: u32,
) -> SemvResult {
    let _span = debug_span!("semv", he = h.raw(), point_id).entered();

    // Create new vertex
    let v_new = topo.vertices.alloc(Vertex { point_id });

    // Read h's data
    let old_h_next = topo.half_edges.get(h).next;
    let old_edge = topo.half_edges.get(h).edge;
    let h_loop = topo.half_edges.get(h).loop_ref;
    let h_twin = topo.half_edges.get(h).twin;

    // Update original edge curve to curve_a (V_a -- V_new)
    topo.edges.get_mut(old_edge).curve_id = curve_a_id;

    // Create new edge for V_new -- V_b segment
    let e_new = topo.edges.alloc(Edge {
        half_edges: [crate::arena::Idx::from_raw(0); 2],
        curve_id: curve_b_id,
    });

    // Create h1_new in h's loop: origin=V_new, next=old_h_next (edge E_new)
    let h1_new = topo.half_edges.alloc(HalfEdge {
        origin: v_new,
        twin: None,
        next: old_h_next,
        edge: e_new,
        loop_ref: h_loop,
    });

    // Wire h → h1_new
    topo.half_edges.get_mut(h).next = h1_new;

    // Handle twin side
    let new_he_b;
    if let Some(t) = h_twin {
        let old_t_next = topo.half_edges.get(t).next;
        let t_loop = topo.half_edges.get(t).loop_ref;

        // Create h2_new in twin's loop: origin=V_new, next=old_t_next (edge E: V_a--V_new)
        let h2_new = topo.half_edges.alloc(HalfEdge {
            origin: v_new,
            twin: Some(h),
            next: old_t_next,
            edge: old_edge,
            loop_ref: t_loop,
        });

        // Wire t → h2_new
        topo.half_edges.get_mut(t).next = h2_new;

        // Set twins: h <-> h2_new (Edge E), h1_new <-> t (Edge E_new)
        topo.half_edges.get_mut(h).twin = Some(h2_new);
        topo.half_edges.get_mut(h1_new).twin = Some(t);
        topo.half_edges.get_mut(t).twin = Some(h1_new);
        topo.half_edges.get_mut(t).edge = e_new;

        // Update edge records
        topo.edges.get_mut(old_edge).half_edges = [h, h2_new];
        topo.edges.get_mut(e_new).half_edges = [h1_new, t];

        new_he_b = Some(h2_new);
    } else {
        // No twin: boundary edge, only split one side.
        topo.edges.get_mut(old_edge).half_edges[0] = h;
        topo.edges.get_mut(e_new).half_edges[0] = h1_new;
        new_he_b = None;
    }

    SemvResult {
        new_vertex: v_new,
        new_edge: e_new,
        new_he_a: h1_new,
        new_he_b,
    }
}

// ─── MEF: Make Edge, Face ────────────────────────────────────────────────────

/// Connect h1.origin (V1) and h2.origin (V2) in the same face loop with a new
/// edge, splitting the face into two.
///
/// After the operation:
/// - Old loop: h1 -> ... -> prev(h2) -> he_a -> h1  (he_a goes V2->V1)
/// - New loop: h2 -> ... -> prev(h1) -> he_b -> h2  (he_b goes V1->V2)
///
/// Euler: E+1, F+1, V unchanged.
pub fn mef(
    topo: &mut TopoStore,
    h1: HalfEdgeIdx,
    h2: HalfEdgeIdx,
    surface_id: u32,
    curve_id: u32,
) -> Result<MefResult, EulerError> {
    let _span = debug_span!("mef", h1 = h1.raw(), h2 = h2.raw()).entered();

    if h1 == h2 {
        return Err(EulerError::SameHalfEdge);
    }

    // Verify same loop
    let loop1 = topo.half_edges.get(h1).loop_ref;
    let loop2 = topo.half_edges.get(h2).loop_ref;
    if loop1 != loop2 {
        return Err(EulerError::DifferentLoops);
    }

    let old_loop = loop1;
    let old_face = topo.loops.get(old_loop).face;
    let shell = topo.faces.get(old_face).shell;

    let v1 = topo.half_edges.get(h1).origin;
    let v2 = topo.half_edges.get(h2).origin;

    // Find predecessors
    let prev_h1 = find_prev_he(topo, old_loop, h1);
    let prev_h2 = find_prev_he(topo, old_loop, h2);

    // Create new edge
    let new_edge = topo.edges.alloc(Edge {
        half_edges: [crate::arena::Idx::from_raw(0); 2],
        curve_id,
    });

    // Create new face and loop for the split-off portion
    let new_face = topo.faces.alloc(Face {
        outer_loop: crate::arena::Idx::from_raw(0),
        surface_id,
        mesh_cache: None,
        shell,
    });
    let new_loop = topo.loops.alloc(Loop {
        half_edge: crate::arena::Idx::from_raw(0),
        face: new_face,
    });
    topo.faces.get_mut(new_face).outer_loop = new_loop;

    // Create two new half-edges:
    // he_a: origin=V2, in old loop, goes V2->V1, next=h1
    let he_a = topo.half_edges.alloc(HalfEdge {
        origin: v2,
        twin: None,
        next: h1,
        edge: new_edge,
        loop_ref: old_loop,
    });
    // he_b: origin=V1, in new loop, goes V1->V2, next=h2
    let he_b = topo.half_edges.alloc(HalfEdge {
        origin: v1,
        twin: None,
        next: h2,
        edge: new_edge,
        loop_ref: new_loop,
    });

    // Set twins
    topo.half_edges.get_mut(he_a).twin = Some(he_b);
    topo.half_edges.get_mut(he_b).twin = Some(he_a);

    // Update edge record
    topo.edges.get_mut(new_edge).half_edges = [he_a, he_b];

    // Rewire loop chains:
    // prev(h2).next = he_a  (old loop ends with: ...prev_h2 -> he_a -> h1 -> ...)
    topo.half_edges.get_mut(prev_h2).next = he_a;
    // prev(h1).next = he_b  (new loop ends with: ...prev_h1 -> he_b -> h2 -> ...)
    topo.half_edges.get_mut(prev_h1).next = he_b;

    // Update loop entries
    topo.loops.get_mut(old_loop).half_edge = he_a;
    topo.loops.get_mut(new_loop).half_edge = he_b;

    // Walk new loop and update loop_ref on all half-edges
    {
        let mut cur = he_b;
        loop {
            topo.half_edges.get_mut(cur).loop_ref = new_loop;
            cur = topo.half_edges.get(cur).next;
            if cur == he_b {
                break;
            }
        }
    }

    // Walk old loop to ensure loop_ref is correct
    {
        let mut cur = he_a;
        loop {
            topo.half_edges.get_mut(cur).loop_ref = old_loop;
            cur = topo.half_edges.get(cur).next;
            if cur == he_a {
                break;
            }
        }
    }

    // Add new face to shell
    topo.shells.get_mut(shell).faces.push(new_face);

    Ok(MefResult {
        new_edge,
        new_face,
        new_loop,
        he_a,
        he_b,
    })
}

// ─── KEV: Kill Edge, Vertex ──────────────────────────────────────────────────

/// Remove edge E and merge V_kill into V_keep. `v_keep` must be one of the
/// edge's endpoints; the other endpoint is killed.
///
/// Euler: V-1, E-1, F unchanged.
pub fn kev(
    topo: &mut TopoStore,
    edge: EdgeIdx,
    v_keep: VertexIdx,
) -> Result<KevResult, EulerError> {
    let _span = debug_span!("kev", edge = edge.raw(), v_keep = v_keep.raw()).entered();

    let [he_0, he_1] = topo.edges.get(edge).half_edges;

    // Determine which vertex to kill
    let origin_0 = topo.half_edges.get(he_0).origin;
    let origin_1 = topo.half_edges.get(he_1).origin;

    let (v_kill, he_from_keep, he_from_kill) = if origin_0 == v_keep {
        // he_0 goes v_keep -> v_kill, he_1 goes v_kill -> v_keep
        (origin_1, he_0, he_1)
    } else if origin_1 == v_keep {
        // he_1 goes v_keep -> v_kill, he_0 goes v_kill -> v_keep
        (origin_0, he_1, he_0)
    } else {
        return Err(EulerError::VertexNotOnEdge(edge, v_keep));
    };

    // Get loop info for both half-edges
    let loop_from_keep = topo.half_edges.get(he_from_keep).loop_ref;
    let loop_from_kill = topo.half_edges.get(he_from_kill).loop_ref;
    let next_from_keep = topo.half_edges.get(he_from_keep).next;
    let next_from_kill = topo.half_edges.get(he_from_kill).next;

    // Bypass he_from_keep: prev.next = he_from_keep.next
    let prev_keep = find_prev_he(topo, loop_from_keep, he_from_keep);
    topo.half_edges.get_mut(prev_keep).next = next_from_keep;

    // Bypass he_from_kill: prev.next = he_from_kill.next
    let prev_kill = find_prev_he(topo, loop_from_kill, he_from_kill);
    topo.half_edges.get_mut(prev_kill).next = next_from_kill;

    // Fix loop entries if they pointed to dead HEs
    if topo.loops.get(loop_from_keep).half_edge == he_from_keep {
        topo.loops.get_mut(loop_from_keep).half_edge = next_from_keep;
    }
    if topo.loops.get(loop_from_kill).half_edge == he_from_kill {
        topo.loops.get_mut(loop_from_kill).half_edge = next_from_kill;
    }

    // Rewrite all HE origins from v_kill -> v_keep across the entire shell
    let face_of_loop = topo.loops.get(loop_from_keep).face;
    let shell_idx = topo.faces.get(face_of_loop).shell;
    let faces: Vec<FaceIdx> = topo.shells.get(shell_idx).faces.clone();

    for &face_idx in &faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start = topo.loops.get(loop_idx).half_edge;
        let mut cur = start;
        loop {
            if topo.half_edges.get(cur).origin == v_kill {
                topo.half_edges.get_mut(cur).origin = v_keep;
            }
            cur = topo.half_edges.get(cur).next;
            if cur == start {
                break;
            }
        }
    }

    Ok(KevResult {
        killed_vertex: v_kill,
        killed_edge: edge,
    })
}

// ─── KEF: Kill Edge, Face ────────────────────────────────────────────────────

/// Remove edge E between f_survive and f_kill, merging the two face loops.
///
/// Euler: E-1, F-1, V unchanged.
pub fn kef(
    topo: &mut TopoStore,
    edge: EdgeIdx,
    f_survive: FaceIdx,
) -> Result<KefResult, EulerError> {
    let _span = debug_span!("kef", edge = edge.raw(), f_survive = f_survive.raw()).entered();

    let [he_0, he_1] = topo.edges.get(edge).half_edges;

    // Determine which half-edge is in f_survive vs f_kill
    let loop_0 = topo.half_edges.get(he_0).loop_ref;
    let loop_1 = topo.half_edges.get(he_1).loop_ref;
    let face_0 = topo.loops.get(loop_0).face;
    let face_1 = topo.loops.get(loop_1).face;

    let (he_survive, he_kill, l_survive, l_kill, f_kill) = if face_0 == f_survive {
        (he_0, he_1, loop_0, loop_1, face_1)
    } else if face_1 == f_survive {
        (he_1, he_0, loop_1, loop_0, face_0)
    } else {
        return Err(EulerError::FaceNotAdjacentToEdge(edge, f_survive));
    };

    let he_survive_next = topo.half_edges.get(he_survive).next;
    let he_kill_next = topo.half_edges.get(he_kill).next;

    // Bypass both half-edges by cross-connecting the loops:
    // prev_survive -> he_kill.next  (skip he_kill, enter kill loop's continuation)
    // prev_kill -> he_survive.next  (skip he_survive, enter survive loop's continuation)
    let prev_survive = find_prev_he(topo, l_survive, he_survive);
    let prev_kill = find_prev_he(topo, l_kill, he_kill);

    topo.half_edges.get_mut(prev_survive).next = he_kill_next;
    topo.half_edges.get_mut(prev_kill).next = he_survive_next;

    // Fix loop entry for survive loop (avoid dead HE)
    if topo.loops.get(l_survive).half_edge == he_survive {
        topo.loops.get_mut(l_survive).half_edge = he_kill_next;
    }

    // Walk merged loop and update all loop_ref -> l_survive
    {
        let start = topo.loops.get(l_survive).half_edge;
        let mut cur = start;
        loop {
            topo.half_edges.get_mut(cur).loop_ref = l_survive;
            cur = topo.half_edges.get(cur).next;
            if cur == start {
                break;
            }
        }
    }

    // Remove f_kill from shell
    let shell_idx = topo.faces.get(f_survive).shell;
    topo.shells.get_mut(shell_idx).faces.retain(|&f| f != f_kill);

    Ok(KefResult {
        killed_face: f_kill,
        killed_edge: edge,
    })
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arena::Idx;
    use crate::diagnostics::validate_solid;
    use std::collections::{HashMap, HashSet};

    /// Build a valid tetrahedron: 4V, 6E, 4F, genus=0, Euler=2.
    /// Edge.half_edges are set correctly for Euler operator use.
    fn make_test_tetrahedron(topo: &mut TopoStore) -> SolidIdx {
        let v: Vec<VertexIdx> = (0..4)
            .map(|i| topo.vertices.alloc(Vertex { point_id: i }))
            .collect();

        let solid_idx = topo.solids.alloc(Solid {
            shell: Idx::from_raw(0),
            genus: 0,
        });
        let shell_idx = topo.shells.alloc(Shell {
            faces: Vec::new(),
            solid: solid_idx,
        });
        topo.solids.get_mut(solid_idx).shell = shell_idx;

        let triangles: [(usize, usize, usize); 4] = [
            (0, 1, 2),
            (0, 3, 1),
            (1, 3, 2),
            (0, 2, 3),
        ];

        let mut face_idxs = Vec::new();
        let mut he_map: HashMap<(u32, u32), HalfEdgeIdx> = HashMap::new();

        for &(a, b, c) in &triangles {
            let face_idx = topo.faces.alloc(Face {
                outer_loop: Idx::from_raw(0),
                surface_id: 0,
                mesh_cache: None,
                shell: shell_idx,
            });
            let loop_idx = topo.loops.alloc(Loop {
                half_edge: Idx::from_raw(0),
                face: face_idx,
            });
            topo.faces.get_mut(face_idx).outer_loop = loop_idx;

            let verts_in_face = [v[a], v[b], v[c]];
            let mut hes = Vec::new();
            for i in 0..3 {
                let origin = verts_in_face[i];
                let edge_idx = topo.edges.alloc(Edge {
                    half_edges: [Idx::from_raw(0); 2],
                    curve_id: 0,
                });
                let he_idx = topo.half_edges.alloc(HalfEdge {
                    origin,
                    twin: None,
                    next: Idx::from_raw(0),
                    edge: edge_idx,
                    loop_ref: loop_idx,
                });
                topo.edges.get_mut(edge_idx).half_edges[0] = he_idx;
                hes.push(he_idx);
            }

            // Wire next pointers
            for i in 0..3 {
                topo.half_edges.get_mut(hes[i]).next = hes[(i + 1) % 3];
            }
            topo.loops.get_mut(loop_idx).half_edge = hes[0];

            // Record for twin matching
            let tri = [a as u32, b as u32, c as u32];
            for i in 0..3 {
                he_map.insert((tri[i], tri[(i + 1) % 3]), hes[i]);
            }

            face_idxs.push(face_idx);
        }

        // Match twins and set edge.half_edges properly
        let keys: Vec<(u32, u32)> = he_map.keys().cloned().collect();
        for (o, d) in keys {
            if let (Some(&he), Some(&twin)) = (he_map.get(&(o, d)), he_map.get(&(d, o))) {
                if topo.half_edges.get(he).twin.is_some() {
                    continue; // already matched
                }
                topo.half_edges.get_mut(he).twin = Some(twin);
                topo.half_edges.get_mut(twin).twin = Some(he);
                // Share edge: twin adopts he's edge
                let edge = topo.half_edges.get(he).edge;
                topo.half_edges.get_mut(twin).edge = edge;
                topo.edges.get_mut(edge).half_edges = [he, twin];
            }
        }

        topo.shells.get_mut(shell_idx).faces = face_idxs;
        solid_idx
    }

    /// Count V, E, F by walking the shell.
    fn count_vef(topo: &TopoStore, solid: SolidIdx) -> (usize, usize, usize) {
        let shell = topo.solids.get(solid).shell;
        let faces = &topo.shells.get(shell).faces;
        let mut verts = HashSet::new();
        let mut edges = HashSet::new();
        for &face_idx in faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                verts.insert(topo.half_edges.get(he).origin.raw());
                edges.insert(topo.half_edges.get(he).edge.raw());
                he = topo.half_edges.get(he).next;
                if he == start {
                    break;
                }
            }
        }
        (verts.len(), edges.len(), faces.len())
    }

    /// Find a half-edge from vertex `from` to vertex `to` (by vertex raw index)
    /// in the given solid.
    fn find_he(topo: &TopoStore, solid: SolidIdx, from: u32, to: u32) -> HalfEdgeIdx {
        let shell = topo.solids.get(solid).shell;
        for &face_idx in &topo.shells.get(shell).faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                let he_data = topo.half_edges.get(he);
                let origin = he_data.origin.raw();
                let dest = topo.half_edges.get(he_data.next).origin.raw();
                if origin == from && dest == to {
                    return he;
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
        panic!("No half-edge from vertex {} to vertex {} found", from, to);
    }

    // ── SEMV tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_semv_tetrahedron() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // Before: 4V, 6E, 4F
        assert_eq!(count_vef(&topo, solid), (4, 6, 4));

        // Split edge 0->1 (half-edge from vertex 0 to vertex 1)
        let h = find_he(&topo, solid, 0, 1);
        let result = semv(&mut topo, h, 99, 0, 0);

        // After: 5V, 7E, 4F
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (5, 7, 4), "SEMV: V+1, E+1, F unchanged");

        // Euler check
        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "After SEMV: {:?}", report.errors());
        assert_eq!(report.euler_value, 2);

        // New vertex exists
        assert_eq!(topo.vertices.get(result.new_vertex).point_id, 99);
    }

    #[test]
    fn test_semv_preserves_twins() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        let h = find_he(&topo, solid, 0, 1);
        let _result = semv(&mut topo, h, 99, 0, 0);

        // Every half-edge in the solid must have a twin
        let shell = topo.solids.get(solid).shell;
        for &face_idx in &topo.shells.get(shell).faces {
            let loop_idx = topo.faces.get(face_idx).outer_loop;
            let start = topo.loops.get(loop_idx).half_edge;
            let mut he = start;
            loop {
                let he_data = topo.half_edges.get(he);
                assert!(
                    he_data.twin.is_some(),
                    "Half-edge Idx({}) has no twin after SEMV",
                    he.raw()
                );
                // Twin consistency
                let twin = he_data.twin.unwrap();
                assert_eq!(
                    topo.half_edges.get(twin).twin,
                    Some(he),
                    "Twin asymmetry at Idx({})",
                    he.raw()
                );
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
    }

    #[test]
    fn test_semv_loop_sizes() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // Edge 0->1 is shared by faces (0,1,2) and (0,3,1).
        // After SEMV, those two faces become quads (4 edges), the other two stay triangles.
        let h = find_he(&topo, solid, 0, 1);
        let h_loop = topo.half_edges.get(h).loop_ref;
        let twin = topo.half_edges.get(h).twin.unwrap();
        let t_loop = topo.half_edges.get(twin).loop_ref;

        let _result = semv(&mut topo, h, 99, 0, 0);

        assert_eq!(loop_edge_count(&topo, h_loop), 4, "h's face becomes a quad");
        assert_eq!(loop_edge_count(&topo, t_loop), 4, "twin's face becomes a quad");
    }

    // ── MEF tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_mef_splits_face() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // First SEMV to create a quad face: split edge 0->1
        let h = find_he(&topo, solid, 0, 1);
        let semv_result = semv(&mut topo, h, 99, 0, 0);

        // After SEMV: 5V, 7E, 4F
        assert_eq!(count_vef(&topo, solid), (5, 7, 4));

        // The face containing h is now a quad: 0 -> V_new -> 1 -> 2 -> (back to 0)
        // Split it diagonally: connect V_new and vertex 2.
        // h1 = semv_result.new_he_a (origin=V_new), find h2 with origin=vertex 2 in same loop
        let h_loop = topo.half_edges.get(h).loop_ref;
        let hes_in_loop = collect_loop_hes(&topo, h_loop);
        let h2 = hes_in_loop.iter().copied().find(|&he| {
            topo.half_edges.get(he).origin.raw() == 2 // vertex 2
        }).expect("vertex 2 should be in this loop");

        let h1 = semv_result.new_he_a; // origin = V_new
        let mef_result = mef(&mut topo, h1, h2, 0, 0).expect("MEF should succeed");

        // After MEF: 5V, 8E, 5F
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (5, 8, 5), "MEF: E+1, F+1, V unchanged");

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "After MEF: {:?}", report.errors());
        assert_eq!(report.euler_value, 2);

        // Both resulting faces should be triangles
        assert_eq!(loop_edge_count(&topo, h_loop), 3);
        assert_eq!(loop_edge_count(&topo, mef_result.new_loop), 3);
    }

    #[test]
    fn test_mef_error_different_loops() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // Pick h1 from face 0 and h2 from face 1 — different loops
        let shell = topo.solids.get(solid).shell;
        let face0 = topo.shells.get(shell).faces[0];
        let face1 = topo.shells.get(shell).faces[1];
        let loop0 = topo.faces.get(face0).outer_loop;
        let loop1 = topo.faces.get(face1).outer_loop;
        let h1 = topo.loops.get(loop0).half_edge;
        let h2 = topo.loops.get(loop1).half_edge;

        let err = mef(&mut topo, h1, h2, 0, 0);
        assert!(matches!(err, Err(EulerError::DifferentLoops)));
    }

    #[test]
    fn test_mef_error_same_he() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        let shell = topo.solids.get(solid).shell;
        let face0 = topo.shells.get(shell).faces[0];
        let loop0 = topo.faces.get(face0).outer_loop;
        let h1 = topo.loops.get(loop0).half_edge;

        let err = mef(&mut topo, h1, h1, 0, 0);
        assert!(matches!(err, Err(EulerError::SameHalfEdge)));
    }

    // ── KEV roundtrip ────────────────────────────────────────────────────────

    #[test]
    fn test_kev_roundtrip() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // Before: 4V, 6E, 4F
        let before = count_vef(&topo, solid);
        assert_eq!(before, (4, 6, 4));

        // SEMV on edge 0->1
        let h = find_he(&topo, solid, 0, 1);
        let semv_result = semv(&mut topo, h, 99, 0, 0);

        // After SEMV: 5V, 7E, 4F
        assert_eq!(count_vef(&topo, solid), (5, 7, 4));

        // KEV on the new edge, keeping V_b (the original destination).
        // The new_he_a (V_new -> V_b) is on edge E_new. V_b's vertex index is 1 (raw).
        // We need to find the vertex index for V_b = original vertex 1.
        // After SEMV, h1_new goes V_new -> V_b. V_b = h1_new.next.origin in the loop? No,
        // V_b is the origin of the next HE after h1_new. But we stored V_b as vertex raw=1.
        // Find a vertex on E_new that isn't V_new:
        let e_new = semv_result.new_edge;
        let [he_a, he_b] = topo.edges.get(e_new).half_edges;
        let origin_a = topo.half_edges.get(he_a).origin;
        let origin_b = topo.half_edges.get(he_b).origin;
        let v_keep = if origin_a == semv_result.new_vertex {
            origin_b
        } else {
            origin_a
        };

        kev(&mut topo, e_new, v_keep).expect("KEV should succeed");

        // After KEV: back to 4V, 6E, 4F
        let after = count_vef(&topo, solid);
        assert_eq!(after, (4, 6, 4), "KEV roundtrip restores counts");

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "After KEV roundtrip: {:?}", report.errors());
    }

    // ── KEF roundtrip ────────────────────────────────────────────────────────

    #[test]
    fn test_kef_roundtrip() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // SEMV to create a quad
        let h = find_he(&topo, solid, 0, 1);
        let semv_result = semv(&mut topo, h, 99, 0, 0);
        assert_eq!(count_vef(&topo, solid), (5, 7, 4));

        // MEF to split the quad into two triangles
        let h_loop = topo.half_edges.get(h).loop_ref;
        let hes_in_loop = collect_loop_hes(&topo, h_loop);
        let h2 = hes_in_loop.iter().copied().find(|&he| {
            topo.half_edges.get(he).origin.raw() == 2
        }).unwrap();
        let h1 = semv_result.new_he_a;

        let mef_result = mef(&mut topo, h1, h2, 0, 0).expect("MEF should succeed");

        // After MEF: 5V, 8E, 5F
        assert_eq!(count_vef(&topo, solid), (5, 8, 5));

        // KEF on the new edge, surviving the old face
        let old_face = topo.loops.get(h_loop).face;
        kef(&mut topo, mef_result.new_edge, old_face).expect("KEF should succeed");

        // After KEF: back to 5V, 7E, 4F
        let after = count_vef(&topo, solid);
        assert_eq!(after, (5, 7, 4), "KEF roundtrip restores counts");

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "After KEF roundtrip: {:?}", report.errors());
    }

    // ── Utility tests ────────────────────────────────────────────────────────

    #[test]
    fn test_collect_loop_hes() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        let shell = topo.solids.get(solid).shell;
        let face0 = topo.shells.get(shell).faces[0];
        let loop0 = topo.faces.get(face0).outer_loop;

        let hes = collect_loop_hes(&topo, loop0);
        assert_eq!(hes.len(), 3, "Triangle face has 3 half-edges");
    }

    #[test]
    fn test_loop_edge_count() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        let shell = topo.solids.get(solid).shell;
        let face0 = topo.shells.get(shell).faces[0];
        let loop0 = topo.faces.get(face0).outer_loop;

        assert_eq!(loop_edge_count(&topo, loop0), 3);
    }

    #[test]
    fn test_find_prev_he() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        let shell = topo.solids.get(solid).shell;
        let face0 = topo.shells.get(shell).faces[0];
        let loop0 = topo.faces.get(face0).outer_loop;
        let start = topo.loops.get(loop0).half_edge;
        let second = topo.half_edges.get(start).next;

        assert_eq!(find_prev_he(&topo, loop0, second), start);
    }

    // ── Compound operations ──────────────────────────────────────────────────

    #[test]
    fn test_compound_semv_mef() {
        // SEMV two edges on the same face, then MEF between the two new vertices
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        // Face (0,1,2): edges 0->1, 1->2, 2->0
        // SEMV on edge 0->1
        let h01 = find_he(&topo, solid, 0, 1);
        let r1 = semv(&mut topo, h01, 90, 0, 0);
        // 5V, 7E, 4F

        // SEMV on edge 1->2
        let h12 = find_he(&topo, solid, 1, 2);
        let r2 = semv(&mut topo, h12, 91, 0, 0);
        // 6V, 8E, 4F
        assert_eq!(count_vef(&topo, solid), (6, 8, 4));

        // The face containing both new vertices is now a pentagon:
        // 0 -> V_new1 -> 1 -> V_new2 -> 2 -> 0
        // MEF between V_new1 and V_new2
        let h1 = r1.new_he_a; // origin = V_new1
        let h_loop = topo.half_edges.get(h1).loop_ref;

        // Find the HE with origin = V_new2 in the same loop
        let hes = collect_loop_hes(&topo, h_loop);
        let h2 = hes.iter().copied().find(|&he| {
            topo.half_edges.get(he).origin == r2.new_vertex
        }).expect("V_new2 should be in the loop");

        mef(&mut topo, h1, h2, 0, 0).expect("MEF should succeed");

        // 6V, 9E, 5F
        let (v, e, f) = count_vef(&topo, solid);
        assert_eq!((v, e, f), (6, 9, 5));

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "Compound SEMV+MEF: {:?}", report.errors());
    }

    #[test]
    fn test_tetrahedron_baseline() {
        let mut topo = TopoStore::new();
        let solid = make_test_tetrahedron(&mut topo);

        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "Baseline tetrahedron: {:?}", report.errors());
        assert_eq!(report.vertex_count, 4);
        assert_eq!(report.edge_count, 6);
        assert_eq!(report.face_count, 4);
        assert_eq!(report.euler_value, 2);
    }
}

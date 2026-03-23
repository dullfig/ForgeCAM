//! Shape evolution tracking for persistent entity naming.
//!
//! Every modeling operation (boolean, chamfer, fillet, extrude, etc.) produces
//! a `ShapeEvolution` that records what happened to each topological entity:
//! which faces/edges/vertices were created, modified, split, or deleted.
//!
//! This is the foundation for persistent entity IDs. The parameter-registry
//! crate will consume these evolution records to maintain a naming journal
//! that tracks entities across feature tree replays.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::topo::{EdgeIdx, FaceIdx, VertexIdx};

/// Where a face in the result solid came from.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FaceOrigin {
    /// Newly created, no predecessor (primitive builder, extrude cap, etc.).
    Primitive,
    /// Same face with modified boundary (e.g., boolean trimmed it, chamfer
    /// changed an adjacent edge). The inner value is the input FaceIdx.
    Modified(FaceIdx),
    /// This face is one piece of a face that was split (e.g., boolean
    /// intersection curve divided a face). The inner value is the original.
    SplitFrom(FaceIdx),
    /// Generated from an edge (e.g., fillet/chamfer creates a blending face
    /// where the edge used to be).
    FromEdge(EdgeIdx),
    /// Generated at the intersection of two faces from different solids
    /// during a boolean operation (the "seam" face along the cut).
    FromIntersection(FaceIdx, FaceIdx),
    /// Copied from another solid without modification (e.g., boolean
    /// keeps a face from solid_a or solid_b unchanged).
    CopiedFrom(FaceIdx),
}

/// Where an edge in the result solid came from.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EdgeOrigin {
    /// Newly created, no predecessor.
    Primitive,
    /// Same edge with modified endpoints or geometry.
    Modified(EdgeIdx),
    /// Created by splitting an existing edge.
    SplitFrom(EdgeIdx),
    /// Created at a face-face intersection during boolean.
    FromIntersection(FaceIdx, FaceIdx),
    /// Copied from another solid without modification.
    CopiedFrom(EdgeIdx),
}

/// Where a vertex in the result solid came from.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum VertexOrigin {
    /// Newly created.
    Primitive,
    /// Same vertex, possibly moved.
    Modified(VertexIdx),
    /// Copied from another solid.
    CopiedFrom(VertexIdx),
}

/// Complete record of what a modeling operation did to the topology.
///
/// For each entity in the **result** solid, the provenance map says where it
/// came from. The deleted lists record which **input** entities no longer
/// exist in the result.
///
/// # Usage
///
/// ```ignore
/// // After a boolean cut:
/// let evo = kernel.last_evolution().unwrap();
///
/// // "What happened to the face my dimension was on?"
/// if evo.deleted_faces.contains(&my_face) {
///     // Face was consumed — annotation is broken
/// } else if let Some(origin) = evo.face_provenance.values()
///     .find(|o| matches!(o, FaceOrigin::Modified(f) if *f == my_face))
/// {
///     // Face survived with modified boundary — annotation is fine
/// }
/// ```
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ShapeEvolution {
    /// For each face in the result solid, where it came from.
    pub face_provenance: HashMap<FaceIdx, FaceOrigin>,
    /// For each edge in the result solid, where it came from.
    pub edge_provenance: HashMap<EdgeIdx, EdgeOrigin>,
    /// For each vertex in the result solid, where it came from.
    pub vertex_provenance: HashMap<VertexIdx, VertexOrigin>,

    /// Faces from the input solid(s) that were removed.
    pub deleted_faces: Vec<FaceIdx>,
    /// Edges from the input solid(s) that were removed.
    pub deleted_edges: Vec<EdgeIdx>,
    /// Vertices from the input solid(s) that were removed.
    pub deleted_vertices: Vec<VertexIdx>,
}

impl ShapeEvolution {
    pub fn new() -> Self {
        Self::default()
    }

    /// Create evolution for a primitive builder — every entity is new.
    pub fn all_primitive(
        faces: impl IntoIterator<Item = FaceIdx>,
        edges: impl IntoIterator<Item = EdgeIdx>,
        vertices: impl IntoIterator<Item = VertexIdx>,
    ) -> Self {
        let mut evo = Self::new();
        for f in faces {
            evo.face_provenance.insert(f, FaceOrigin::Primitive);
        }
        for e in edges {
            evo.edge_provenance.insert(e, EdgeOrigin::Primitive);
        }
        for v in vertices {
            evo.vertex_provenance.insert(v, VertexOrigin::Primitive);
        }
        evo
    }

    /// Record that a result face came from the given origin.
    pub fn record_face(&mut self, result: FaceIdx, origin: FaceOrigin) {
        self.face_provenance.insert(result, origin);
    }

    /// Record that a result edge came from the given origin.
    pub fn record_edge(&mut self, result: EdgeIdx, origin: EdgeOrigin) {
        self.edge_provenance.insert(result, origin);
    }

    /// Record that a result vertex came from the given origin.
    pub fn record_vertex(&mut self, result: VertexIdx, origin: VertexOrigin) {
        self.vertex_provenance.insert(result, origin);
    }

    /// Record that an input face was deleted (no longer in result).
    pub fn record_deleted_face(&mut self, face: FaceIdx) {
        self.deleted_faces.push(face);
    }

    /// Record that an input edge was deleted.
    pub fn record_deleted_edge(&mut self, edge: EdgeIdx) {
        self.deleted_edges.push(edge);
    }

    /// Number of entities tracked.
    pub fn entity_count(&self) -> usize {
        self.face_provenance.len()
            + self.edge_provenance.len()
            + self.vertex_provenance.len()
    }

    /// True if no evolution was recorded (empty/unused).
    pub fn is_empty(&self) -> bool {
        self.face_provenance.is_empty()
            && self.edge_provenance.is_empty()
            && self.vertex_provenance.is_empty()
            && self.deleted_faces.is_empty()
            && self.deleted_edges.is_empty()
            && self.deleted_vertices.is_empty()
    }

    /// Find all result faces that originated from a given input face.
    ///
    /// Returns faces marked as Modified, SplitFrom, or CopiedFrom the input.
    pub fn descendants_of_face(&self, input: FaceIdx) -> Vec<FaceIdx> {
        self.face_provenance
            .iter()
            .filter(|(_, origin)| match origin {
                FaceOrigin::Modified(f) | FaceOrigin::SplitFrom(f) | FaceOrigin::CopiedFrom(f) => {
                    *f == input
                }
                _ => false,
            })
            .map(|(&result, _)| result)
            .collect()
    }

    /// Find all result edges that originated from a given input edge.
    pub fn descendants_of_edge(&self, input: EdgeIdx) -> Vec<EdgeIdx> {
        self.edge_provenance
            .iter()
            .filter(|(_, origin)| match origin {
                EdgeOrigin::Modified(e) | EdgeOrigin::SplitFrom(e) | EdgeOrigin::CopiedFrom(e) => {
                    *e == input
                }
                _ => false,
            })
            .map(|(&result, _)| result)
            .collect()
    }

    /// Check whether an input face survived (has any descendant in the result).
    pub fn face_survived(&self, input: FaceIdx) -> bool {
        !self.descendants_of_face(input).is_empty()
    }

    /// Check whether an input face was deleted.
    pub fn face_was_deleted(&self, input: FaceIdx) -> bool {
        self.deleted_faces.contains(&input)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primitive_evolution() {
        let f1 = FaceIdx::from_raw(0);
        let f2 = FaceIdx::from_raw(1);
        let e1 = EdgeIdx::from_raw(0);
        let v1 = VertexIdx::from_raw(0);

        let evo = ShapeEvolution::all_primitive([f1, f2], [e1], [v1]);
        assert_eq!(evo.face_provenance.len(), 2);
        assert_eq!(evo.edge_provenance.len(), 1);
        assert_eq!(evo.vertex_provenance.len(), 1);
        assert!(matches!(evo.face_provenance[&f1], FaceOrigin::Primitive));
        assert!(evo.deleted_faces.is_empty());
    }

    #[test]
    fn test_split_tracking() {
        let mut evo = ShapeEvolution::new();
        let original = FaceIdx::from_raw(5);
        let split_a = FaceIdx::from_raw(10);
        let split_b = FaceIdx::from_raw(11);

        evo.record_face(split_a, FaceOrigin::SplitFrom(original));
        evo.record_face(split_b, FaceOrigin::SplitFrom(original));
        evo.record_deleted_face(original);

        assert!(evo.face_was_deleted(original));
        let desc = evo.descendants_of_face(original);
        assert_eq!(desc.len(), 2);
        assert!(desc.contains(&split_a));
        assert!(desc.contains(&split_b));
    }

    #[test]
    fn test_modification_tracking() {
        let mut evo = ShapeEvolution::new();
        let old_face = FaceIdx::from_raw(3);
        let new_face = FaceIdx::from_raw(20);

        evo.record_face(new_face, FaceOrigin::Modified(old_face));

        assert!(evo.face_survived(old_face));
        assert!(!evo.face_was_deleted(old_face));
        assert_eq!(evo.descendants_of_face(old_face), vec![new_face]);
    }

    #[test]
    fn test_fillet_evolution() {
        let mut evo = ShapeEvolution::new();
        let filleted_edge = EdgeIdx::from_raw(7);
        let fillet_face = FaceIdx::from_raw(15);
        let adj_face_1 = FaceIdx::from_raw(2);
        let adj_face_1_new = FaceIdx::from_raw(16);
        let adj_face_2 = FaceIdx::from_raw(4);
        let adj_face_2_new = FaceIdx::from_raw(17);

        // Fillet creates a new face from the edge
        evo.record_face(fillet_face, FaceOrigin::FromEdge(filleted_edge));
        // Adjacent faces get modified boundaries
        evo.record_face(adj_face_1_new, FaceOrigin::Modified(adj_face_1));
        evo.record_face(adj_face_2_new, FaceOrigin::Modified(adj_face_2));
        // The edge is consumed
        evo.record_deleted_edge(filleted_edge);

        assert!(evo.face_survived(adj_face_1));
        assert!(evo.face_survived(adj_face_2));
        assert!(evo.deleted_edges.contains(&filleted_edge));
        assert!(matches!(
            evo.face_provenance[&fillet_face],
            FaceOrigin::FromEdge(e) if e == filleted_edge
        ));
    }

    #[test]
    fn test_boolean_evolution() {
        let mut evo = ShapeEvolution::new();

        // Faces from solid A that survived (copied into result)
        let a_top = FaceIdx::from_raw(1);
        let a_top_result = FaceIdx::from_raw(30);
        evo.record_face(a_top_result, FaceOrigin::CopiedFrom(a_top));

        // Face from solid A that was split by intersection
        let a_front = FaceIdx::from_raw(2);
        let a_front_outer = FaceIdx::from_raw(31);
        let a_front_inner = FaceIdx::from_raw(32);
        evo.record_face(a_front_outer, FaceOrigin::SplitFrom(a_front));
        evo.record_face(a_front_inner, FaceOrigin::SplitFrom(a_front));

        // Face from solid B that was deleted (inside solid A, cut away)
        let b_top = FaceIdx::from_raw(10);
        evo.record_deleted_face(b_top);

        assert!(evo.face_survived(a_top));
        assert!(evo.face_survived(a_front));
        assert!(!evo.face_survived(b_top));
        assert!(evo.face_was_deleted(b_top));

        let front_desc = evo.descendants_of_face(a_front);
        assert_eq!(front_desc.len(), 2);
    }
}

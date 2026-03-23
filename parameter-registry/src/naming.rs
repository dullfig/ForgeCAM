//! Naming journal — maps persistent `EntityRef` handles to live arena indices.
//!
//! The journal processes [`ShapeEvolution`] records from kernel operations and
//! maintains stable identity for every topological entity across feature-tree replays.

use std::collections::HashMap;

use forgecam_geometry::entity_ref::{EntityKind, EntityRef};
use rustkernel_topology::evolution::{
    EdgeOrigin, FaceOrigin, ShapeEvolution, VertexOrigin,
};
use rustkernel_topology::topo::{EdgeIdx, FaceIdx, VertexIdx};
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// FeatureId
// ---------------------------------------------------------------------------

/// Persistent identifier for a feature in the parametric tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct FeatureId(pub u32);

// ---------------------------------------------------------------------------
// NamingKey
// ---------------------------------------------------------------------------

/// Deterministic name derived from the creation context of an entity.
///
/// Two replays of the same feature tree produce identical `NamingKey` values,
/// which lets the journal map them to the same `EntityRef`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NamingKey {
    /// Entity created by a primitive builder (box, cylinder, etc.).
    Primitive {
        feature: FeatureId,
        kind: EntityKind,
        ordinal: u32,
    },
    /// Face generated from an edge (fillet/chamfer blending face).
    FromEdge {
        feature: FeatureId,
        source_edge_key: Box<NamingKey>,
    },
    /// Face generated at the intersection of two faces during a boolean.
    FromIntersection {
        feature: FeatureId,
        face_a_key: Box<NamingKey>,
        face_b_key: Box<NamingKey>,
    },
    /// One piece of a face that was split (e.g., boolean intersection curve).
    SplitFrom {
        feature: FeatureId,
        parent_key: Box<NamingKey>,
        split_ordinal: u32,
    },
}

// ---------------------------------------------------------------------------
// ResolveResult
// ---------------------------------------------------------------------------

/// Result of resolving an `EntityRef` to a live arena index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResolveResult<Idx> {
    /// The entity exists at this index.
    Live(Idx),
    /// The entity was split by a feature into multiple children.
    Split {
        by_feature: FeatureId,
        children: Vec<(EntityRef, Idx)>,
    },
    /// The entity was consumed/deleted by a feature.
    Deleted { by_feature: FeatureId },
    /// The reference is not known to this journal (stale or from another model).
    Unknown,
}

// ---------------------------------------------------------------------------
// NamingJournal
// ---------------------------------------------------------------------------

/// The runtime mapping engine that maintains stable entity identity.
///
/// After every kernel operation, the document layer calls [`process_evolution`] with
/// the operation's [`ShapeEvolution`]. The journal allocates or reuses [`EntityRef`]
/// handles and updates the live index mappings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamingJournal {
    /// Source of truth: deterministic naming key → stable handle.
    key_to_ref: HashMap<NamingKey, EntityRef>,

    // Forward maps: EntityRef → current live index (None = deleted).
    face_map: HashMap<EntityRef, Option<FaceIdx>>,
    edge_map: HashMap<EntityRef, Option<EdgeIdx>>,
    vertex_map: HashMap<EntityRef, Option<VertexIdx>>,

    // Reverse maps: live index → EntityRef.
    face_reverse: HashMap<FaceIdx, EntityRef>,
    edge_reverse: HashMap<EdgeIdx, EntityRef>,
    vertex_reverse: HashMap<VertexIdx, EntityRef>,

    // Deletion tracking: which feature deleted an entity.
    deleted_by: HashMap<EntityRef, FeatureId>,

    // Split tracking: parent → (feature, children).
    split_info: HashMap<EntityRef, (FeatureId, Vec<EntityRef>)>,

    /// Monotonic counter for allocating new EntityRef handles.
    next_counter: u64,
}

impl NamingJournal {
    /// Create an empty journal.
    pub fn new() -> Self {
        Self {
            key_to_ref: HashMap::new(),
            face_map: HashMap::new(),
            edge_map: HashMap::new(),
            vertex_map: HashMap::new(),
            face_reverse: HashMap::new(),
            edge_reverse: HashMap::new(),
            vertex_reverse: HashMap::new(),
            deleted_by: HashMap::new(),
            split_info: HashMap::new(),
            next_counter: 0,
        }
    }

    /// Allocate a new `EntityRef`, or reuse one if the key already exists.
    fn alloc_or_reuse(&mut self, kind: EntityKind, key: NamingKey) -> EntityRef {
        if let Some(&existing) = self.key_to_ref.get(&key) {
            return existing;
        }
        let er = EntityRef::new(kind, self.next_counter);
        self.next_counter += 1;
        self.key_to_ref.insert(key, er);
        er
    }

    /// Allocate a brand-new `EntityRef` (no key reuse — used for splits where
    /// the key includes a split_ordinal that is always new).
    fn alloc_fresh(&mut self, kind: EntityKind, key: NamingKey) -> EntityRef {
        // For splits, the key includes the ordinal, so alloc_or_reuse works correctly.
        self.alloc_or_reuse(kind, key)
    }

    /// Process a [`ShapeEvolution`] from a kernel operation.
    ///
    /// This is the core method. Call it after every kernel operation with the
    /// feature ID and the evolution record. It updates all internal mappings.
    pub fn process_evolution(&mut self, feature: FeatureId, evo: &ShapeEvolution) {
        // Phase 1: Process face provenance.
        // We need to build a temporary reverse map of old FaceIdx → EntityRef
        // from the current state, so we can look up parents.
        let old_face_reverse = self.face_reverse.clone();
        let old_edge_reverse = self.edge_reverse.clone();
        let old_vertex_reverse = self.vertex_reverse.clone();

        // Remove reverse map entries touched by this evolution.
        // We do NOT clear the entire map — other solids may have live entries
        // that are unrelated to this operation.
        for (&face_idx, origin) in &evo.face_provenance {
            self.face_reverse.remove(&face_idx);
            // Also remove old source entries (parent has moved to new index).
            match origin {
                FaceOrigin::Modified(old) | FaceOrigin::CopiedFrom(old) | FaceOrigin::SplitFrom(old) => {
                    self.face_reverse.remove(old);
                }
                _ => {}
            }
        }
        for &deleted_face in &evo.deleted_faces {
            self.face_reverse.remove(&deleted_face);
        }
        for (&edge_idx, origin) in &evo.edge_provenance {
            self.edge_reverse.remove(&edge_idx);
            match origin {
                EdgeOrigin::Modified(old) | EdgeOrigin::CopiedFrom(old) | EdgeOrigin::SplitFrom(old) => {
                    self.edge_reverse.remove(old);
                }
                _ => {}
            }
        }
        for &deleted_edge in &evo.deleted_edges {
            self.edge_reverse.remove(&deleted_edge);
        }
        for (&vert_idx, origin) in &evo.vertex_provenance {
            self.vertex_reverse.remove(&vert_idx);
            match origin {
                VertexOrigin::Modified(old) | VertexOrigin::CopiedFrom(old) => {
                    self.vertex_reverse.remove(old);
                }
                _ => {}
            }
        }
        for &deleted_vert in &evo.deleted_vertices {
            self.vertex_reverse.remove(&deleted_vert);
        }

        // Track which EntityRefs got split so we can record split_info.
        let mut split_children: HashMap<EntityRef, Vec<(EntityRef, FaceIdx)>> = HashMap::new();

        // Assign ordinals to primitive faces for deterministic naming.
        let mut face_prim_ordinal: u32 = 0;
        let mut edge_prim_ordinal: u32 = 0;
        let mut vertex_prim_ordinal: u32 = 0;

        // Count splits per parent for ordinal assignment.
        let mut split_counts: HashMap<FaceIdx, u32> = HashMap::new();
        for origin in evo.face_provenance.values() {
            if let FaceOrigin::SplitFrom(parent) = origin {
                *split_counts.entry(*parent).or_insert(0) += 1;
            }
        }
        // Track ordinal assignment per parent during iteration.
        let mut split_ordinal_counters: HashMap<FaceIdx, u32> = HashMap::new();

        // Sort face provenance by FaceIdx raw value for deterministic ordering.
        let mut face_entries: Vec<_> = evo.face_provenance.iter().collect();
        face_entries.sort_by_key(|(idx, _)| idx.raw());

        for (&result_idx, origin) in &face_entries {
            match origin {
                FaceOrigin::Primitive => {
                    let key = NamingKey::Primitive {
                        feature,
                        kind: EntityKind::Face,
                        ordinal: face_prim_ordinal,
                    };
                    face_prim_ordinal += 1;
                    let er = self.alloc_or_reuse(EntityKind::Face, key);
                    self.face_map.insert(er, Some(result_idx));
                    self.face_reverse.insert(result_idx, er);
                }
                FaceOrigin::Modified(old_idx) | FaceOrigin::CopiedFrom(old_idx) => {
                    // Inherit parent's EntityRef.
                    if let Some(&parent_ref) = old_face_reverse.get(old_idx) {
                        self.face_map.insert(parent_ref, Some(result_idx));
                        self.face_reverse.insert(result_idx, parent_ref);
                    }
                    // If parent not found in reverse map, this is a cross-solid copy
                    // (boolean). We skip — the entity keeps whatever ref it had.
                }
                FaceOrigin::SplitFrom(old_idx) => {
                    if let Some(&parent_ref) = old_face_reverse.get(old_idx) {
                        let ordinal = split_ordinal_counters.entry(*old_idx).or_insert(0);
                        let parent_key = self.find_key_for_ref(parent_ref).cloned();
                        let key = if let Some(pk) = parent_key {
                            NamingKey::SplitFrom {
                                feature,
                                parent_key: Box::new(pk),
                                split_ordinal: *ordinal,
                            }
                        } else {
                            // Fallback: use a primitive-style key.
                            NamingKey::SplitFrom {
                                feature,
                                parent_key: Box::new(NamingKey::Primitive {
                                    feature: FeatureId(u32::MAX),
                                    kind: EntityKind::Face,
                                    ordinal: old_idx.raw(),
                                }),
                                split_ordinal: *ordinal,
                            }
                        };
                        *ordinal += 1;

                        let child_ref = self.alloc_fresh(EntityKind::Face, key);
                        self.face_map.insert(child_ref, Some(result_idx));
                        self.face_reverse.insert(result_idx, child_ref);

                        // Mark parent as dead (split) and record children.
                        self.face_map.insert(parent_ref, None);
                        split_children
                            .entry(parent_ref)
                            .or_default()
                            .push((child_ref, result_idx));
                    }
                }
                FaceOrigin::FromEdge(edge_idx) => {
                    let source_edge_key = old_edge_reverse
                        .get(edge_idx)
                        .and_then(|er| self.find_key_for_ref(*er).cloned());
                    let key = if let Some(ek) = source_edge_key {
                        NamingKey::FromEdge {
                            feature,
                            source_edge_key: Box::new(ek),
                        }
                    } else {
                        // Edge not tracked yet — use a primitive-style fallback.
                        NamingKey::FromEdge {
                            feature,
                            source_edge_key: Box::new(NamingKey::Primitive {
                                feature: FeatureId(u32::MAX),
                                kind: EntityKind::Edge,
                                ordinal: edge_idx.raw(),
                            }),
                        }
                    };
                    let er = self.alloc_or_reuse(EntityKind::Face, key);
                    self.face_map.insert(er, Some(result_idx));
                    self.face_reverse.insert(result_idx, er);
                }
                FaceOrigin::FromIntersection(fa, fb) => {
                    let key_a = old_face_reverse
                        .get(fa)
                        .and_then(|er| self.find_key_for_ref(*er).cloned());
                    let key_b = old_face_reverse
                        .get(fb)
                        .and_then(|er| self.find_key_for_ref(*er).cloned());
                    let key = NamingKey::FromIntersection {
                        feature,
                        face_a_key: Box::new(key_a.unwrap_or(NamingKey::Primitive {
                            feature: FeatureId(u32::MAX),
                            kind: EntityKind::Face,
                            ordinal: fa.raw(),
                        })),
                        face_b_key: Box::new(key_b.unwrap_or(NamingKey::Primitive {
                            feature: FeatureId(u32::MAX),
                            kind: EntityKind::Face,
                            ordinal: fb.raw(),
                        })),
                    };
                    let er = self.alloc_or_reuse(EntityKind::Face, key);
                    self.face_map.insert(er, Some(result_idx));
                    self.face_reverse.insert(result_idx, er);
                }
            }
        }

        // Record split info.
        for (parent_ref, children) in split_children {
            self.split_info.insert(parent_ref, (feature, children.iter().map(|(er, _)| *er).collect()));
        }

        // Process deleted faces.
        for &deleted_idx in &evo.deleted_faces {
            if let Some(&dead_ref) = old_face_reverse.get(&deleted_idx) {
                // Only mark as deleted if not already handled by split.
                if !self.split_info.contains_key(&dead_ref) {
                    self.face_map.insert(dead_ref, None);
                    self.deleted_by.insert(dead_ref, feature);
                }
            }
        }

        // Phase 2: Process edge provenance.
        let mut edge_entries: Vec<_> = evo.edge_provenance.iter().collect();
        edge_entries.sort_by_key(|(idx, _)| idx.raw());

        for (&result_idx, origin) in &edge_entries {
            match origin {
                EdgeOrigin::Primitive => {
                    let key = NamingKey::Primitive {
                        feature,
                        kind: EntityKind::Edge,
                        ordinal: edge_prim_ordinal,
                    };
                    edge_prim_ordinal += 1;
                    let er = self.alloc_or_reuse(EntityKind::Edge, key);
                    self.edge_map.insert(er, Some(result_idx));
                    self.edge_reverse.insert(result_idx, er);
                }
                EdgeOrigin::Modified(old_idx) | EdgeOrigin::CopiedFrom(old_idx) => {
                    if let Some(&parent_ref) = old_edge_reverse.get(old_idx) {
                        self.edge_map.insert(parent_ref, Some(result_idx));
                        self.edge_reverse.insert(result_idx, parent_ref);
                    }
                }
                EdgeOrigin::SplitFrom(old_idx) => {
                    if let Some(&parent_ref) = old_edge_reverse.get(old_idx) {
                        let parent_key = self.find_key_for_ref(parent_ref).cloned();
                        // Use a simple ordinal based on result index for determinism.
                        let key = NamingKey::SplitFrom {
                            feature,
                            parent_key: Box::new(parent_key.unwrap_or(NamingKey::Primitive {
                                feature: FeatureId(u32::MAX),
                                kind: EntityKind::Edge,
                                ordinal: old_idx.raw(),
                            })),
                            split_ordinal: result_idx.raw(),
                        };
                        let child_ref = self.alloc_fresh(EntityKind::Edge, key);
                        self.edge_map.insert(child_ref, Some(result_idx));
                        self.edge_reverse.insert(result_idx, child_ref);
                        self.edge_map.insert(parent_ref, None);
                    }
                }
                EdgeOrigin::FromIntersection(fa, fb) => {
                    let key_a = old_face_reverse
                        .get(fa)
                        .and_then(|er| self.find_key_for_ref(*er).cloned());
                    let key_b = old_face_reverse
                        .get(fb)
                        .and_then(|er| self.find_key_for_ref(*er).cloned());
                    let key = NamingKey::FromIntersection {
                        feature,
                        face_a_key: Box::new(key_a.unwrap_or(NamingKey::Primitive {
                            feature: FeatureId(u32::MAX),
                            kind: EntityKind::Face,
                            ordinal: fa.raw(),
                        })),
                        face_b_key: Box::new(key_b.unwrap_or(NamingKey::Primitive {
                            feature: FeatureId(u32::MAX),
                            kind: EntityKind::Face,
                            ordinal: fb.raw(),
                        })),
                    };
                    let er = self.alloc_or_reuse(EntityKind::Edge, key);
                    self.edge_map.insert(er, Some(result_idx));
                    self.edge_reverse.insert(result_idx, er);
                }
            }
        }

        for &deleted_idx in &evo.deleted_edges {
            if let Some(&dead_ref) = old_edge_reverse.get(&deleted_idx) {
                self.edge_map.insert(dead_ref, None);
                self.deleted_by.insert(dead_ref, feature);
            }
        }

        // Phase 3: Process vertex provenance.
        let mut vertex_entries: Vec<_> = evo.vertex_provenance.iter().collect();
        vertex_entries.sort_by_key(|(idx, _)| idx.raw());

        for (&result_idx, origin) in &vertex_entries {
            match origin {
                VertexOrigin::Primitive => {
                    let key = NamingKey::Primitive {
                        feature,
                        kind: EntityKind::Vertex,
                        ordinal: vertex_prim_ordinal,
                    };
                    vertex_prim_ordinal += 1;
                    let er = self.alloc_or_reuse(EntityKind::Vertex, key);
                    self.vertex_map.insert(er, Some(result_idx));
                    self.vertex_reverse.insert(result_idx, er);
                }
                VertexOrigin::Modified(old_idx) | VertexOrigin::CopiedFrom(old_idx) => {
                    if let Some(&parent_ref) = old_vertex_reverse.get(old_idx) {
                        self.vertex_map.insert(parent_ref, Some(result_idx));
                        self.vertex_reverse.insert(result_idx, parent_ref);
                    }
                }
            }
        }

        for &deleted_idx in &evo.deleted_vertices {
            if let Some(&dead_ref) = old_vertex_reverse.get(&deleted_idx) {
                self.vertex_map.insert(dead_ref, None);
                self.deleted_by.insert(dead_ref, feature);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Resolution API
    // -----------------------------------------------------------------------

    /// Resolve an `EntityRef` for a face to its current live `FaceIdx`.
    pub fn resolve_face(&self, entity: EntityRef) -> ResolveResult<FaceIdx> {
        match self.face_map.get(&entity) {
            Some(Some(idx)) => ResolveResult::Live(*idx),
            Some(None) => {
                // Check if it was split.
                if let Some((by_feature, children)) = self.split_info.get(&entity) {
                    let resolved: Vec<_> = children
                        .iter()
                        .filter_map(|&child_ref| {
                            self.face_map
                                .get(&child_ref)
                                .and_then(|opt| opt.map(|idx| (child_ref, idx)))
                        })
                        .collect();
                    if !resolved.is_empty() {
                        return ResolveResult::Split {
                            by_feature: *by_feature,
                            children: resolved,
                        };
                    }
                }
                // Deleted.
                if let Some(&by_feature) = self.deleted_by.get(&entity) {
                    ResolveResult::Deleted { by_feature }
                } else {
                    ResolveResult::Deleted {
                        by_feature: FeatureId(0),
                    }
                }
            }
            None => ResolveResult::Unknown,
        }
    }

    /// Resolve an `EntityRef` for an edge to its current live `EdgeIdx`.
    pub fn resolve_edge(&self, entity: EntityRef) -> ResolveResult<EdgeIdx> {
        match self.edge_map.get(&entity) {
            Some(Some(idx)) => ResolveResult::Live(*idx),
            Some(None) => {
                if let Some(&by_feature) = self.deleted_by.get(&entity) {
                    ResolveResult::Deleted { by_feature }
                } else {
                    ResolveResult::Deleted {
                        by_feature: FeatureId(0),
                    }
                }
            }
            None => ResolveResult::Unknown,
        }
    }

    /// Resolve an `EntityRef` for a vertex to its current live `VertexIdx`.
    pub fn resolve_vertex(&self, entity: EntityRef) -> ResolveResult<VertexIdx> {
        match self.vertex_map.get(&entity) {
            Some(Some(idx)) => ResolveResult::Live(*idx),
            Some(None) => {
                if let Some(&by_feature) = self.deleted_by.get(&entity) {
                    ResolveResult::Deleted { by_feature }
                } else {
                    ResolveResult::Deleted {
                        by_feature: FeatureId(0),
                    }
                }
            }
            None => ResolveResult::Unknown,
        }
    }

    // -----------------------------------------------------------------------
    // Lookup helpers
    // -----------------------------------------------------------------------

    /// Look up the EntityRef for a live face index.
    pub fn entity_for_face(&self, idx: FaceIdx) -> Option<EntityRef> {
        self.face_reverse.get(&idx).copied()
    }

    /// Look up the EntityRef for a live edge index.
    pub fn entity_for_edge(&self, idx: EdgeIdx) -> Option<EntityRef> {
        self.edge_reverse.get(&idx).copied()
    }

    /// Look up the EntityRef for a live vertex index.
    pub fn entity_for_vertex(&self, idx: VertexIdx) -> Option<EntityRef> {
        self.vertex_reverse.get(&idx).copied()
    }

    /// Number of tracked face entities (live + dead).
    pub fn face_count(&self) -> usize {
        self.face_map.len()
    }

    /// Number of tracked edge entities (live + dead).
    pub fn edge_count(&self) -> usize {
        self.edge_map.len()
    }

    /// Number of tracked vertex entities (live + dead).
    pub fn vertex_count(&self) -> usize {
        self.vertex_map.len()
    }

    /// Number of live face mappings.
    pub fn live_face_count(&self) -> usize {
        self.face_map.values().filter(|v| v.is_some()).count()
    }

    /// Clear all state. Used before feature-tree replay.
    pub fn clear(&mut self) {
        // Keep key_to_ref — that's the stable mapping we want to survive replay.
        // Clear everything else so process_evolution rebuilds from scratch.
        self.face_map.clear();
        self.edge_map.clear();
        self.vertex_map.clear();
        self.face_reverse.clear();
        self.edge_reverse.clear();
        self.vertex_reverse.clear();
        self.deleted_by.clear();
        self.split_info.clear();
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Find the NamingKey for a given EntityRef (reverse lookup in key_to_ref).
    fn find_key_for_ref(&self, er: EntityRef) -> Option<&NamingKey> {
        self.key_to_ref
            .iter()
            .find(|(_, &v)| v == er)
            .map(|(k, _)| k)
    }
}

impl Default for NamingJournal {
    fn default() -> Self {
        Self::new()
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a ShapeEvolution for a box (6 faces, 12 edges, 8 vertices, all primitive).
    fn box_evolution() -> ShapeEvolution {
        let faces = (0..6).map(FaceIdx::from_raw);
        let edges = (0..12).map(EdgeIdx::from_raw);
        let vertices = (0..8).map(VertexIdx::from_raw);
        ShapeEvolution::all_primitive(faces, edges, vertices)
    }

    #[test]
    fn primitive_box_allocates_six_face_refs() {
        let mut journal = NamingJournal::new();
        let evo = box_evolution();
        journal.process_evolution(FeatureId(1), &evo);

        assert_eq!(journal.live_face_count(), 6);
        assert_eq!(journal.edge_count(), 12);
        assert_eq!(journal.vertex_count(), 8);

        // Every face should resolve.
        for i in 0..6 {
            let idx = FaceIdx::from_raw(i);
            let er = journal.entity_for_face(idx).expect("face should be tracked");
            assert_eq!(er.kind(), EntityKind::Face);
            match journal.resolve_face(er) {
                ResolveResult::Live(resolved) => assert_eq!(resolved, idx),
                other => panic!("expected Live, got {:?}", other),
            }
        }
    }

    #[test]
    fn modified_face_keeps_entity_ref() {
        let mut journal = NamingJournal::new();

        // Step 1: Create a box.
        let box_evo = box_evolution();
        journal.process_evolution(FeatureId(1), &box_evo);

        // Remember the EntityRef for face 0.
        let original_ref = journal
            .entity_for_face(FaceIdx::from_raw(0))
            .expect("face 0 tracked");

        // Step 2: Chamfer modifies face 0 (now at index 20) and creates a new face from edge 3.
        let mut chamfer_evo = ShapeEvolution::new();
        // Face 0 modified → now at index 20.
        chamfer_evo.record_face(FaceIdx::from_raw(20), FaceOrigin::Modified(FaceIdx::from_raw(0)));
        // Faces 1-5 copied.
        for i in 1..6 {
            chamfer_evo.record_face(
                FaceIdx::from_raw(20 + i),
                FaceOrigin::CopiedFrom(FaceIdx::from_raw(i)),
            );
        }
        // New chamfer face from edge 3.
        chamfer_evo.record_face(FaceIdx::from_raw(26), FaceOrigin::FromEdge(EdgeIdx::from_raw(3)));
        // Edges: most modified/copied, the chamfered edge is deleted.
        for i in 0..12 {
            if i == 3 {
                chamfer_evo.record_deleted_edge(EdgeIdx::from_raw(i));
            } else {
                chamfer_evo.record_edge(
                    EdgeIdx::from_raw(20 + i),
                    EdgeOrigin::Modified(EdgeIdx::from_raw(i)),
                );
            }
        }
        // Vertices copied.
        for i in 0..8 {
            chamfer_evo.record_vertex(
                VertexIdx::from_raw(20 + i),
                VertexOrigin::CopiedFrom(VertexIdx::from_raw(i)),
            );
        }

        journal.process_evolution(FeatureId(2), &chamfer_evo);

        // The modified face should keep the same EntityRef.
        let after_ref = journal
            .entity_for_face(FaceIdx::from_raw(20))
            .expect("modified face tracked");
        assert_eq!(original_ref, after_ref, "modified face must keep its EntityRef");

        // The new chamfer face should have a different EntityRef.
        let chamfer_face_ref = journal
            .entity_for_face(FaceIdx::from_raw(26))
            .expect("chamfer face tracked");
        assert_ne!(chamfer_face_ref, original_ref);

        // Resolving the original ref should give the new index.
        match journal.resolve_face(original_ref) {
            ResolveResult::Live(idx) => assert_eq!(idx, FaceIdx::from_raw(20)),
            other => panic!("expected Live, got {:?}", other),
        }
    }

    #[test]
    fn split_produces_distinct_children() {
        let mut journal = NamingJournal::new();

        // Step 1: Create a box.
        let box_evo = box_evolution();
        journal.process_evolution(FeatureId(1), &box_evo);

        let original_ref = journal
            .entity_for_face(FaceIdx::from_raw(2))
            .expect("face 2 tracked");

        // Step 2: Boolean splits face 2 into two pieces.
        let mut bool_evo = ShapeEvolution::new();
        // Faces 0,1 survive as copies.
        bool_evo.record_face(FaceIdx::from_raw(30), FaceOrigin::CopiedFrom(FaceIdx::from_raw(0)));
        bool_evo.record_face(FaceIdx::from_raw(31), FaceOrigin::CopiedFrom(FaceIdx::from_raw(1)));
        // Face 2 splits into 32 and 33.
        bool_evo.record_face(FaceIdx::from_raw(32), FaceOrigin::SplitFrom(FaceIdx::from_raw(2)));
        bool_evo.record_face(FaceIdx::from_raw(33), FaceOrigin::SplitFrom(FaceIdx::from_raw(2)));
        // Faces 3-5 survive (indices 34, 35, 36 — no overlap with splits).
        bool_evo.record_face(FaceIdx::from_raw(34), FaceOrigin::CopiedFrom(FaceIdx::from_raw(3)));
        bool_evo.record_face(FaceIdx::from_raw(35), FaceOrigin::CopiedFrom(FaceIdx::from_raw(4)));
        bool_evo.record_face(FaceIdx::from_raw(36), FaceOrigin::CopiedFrom(FaceIdx::from_raw(5)));
        // Intersection seam face.
        bool_evo.record_face(
            FaceIdx::from_raw(37),
            FaceOrigin::FromIntersection(FaceIdx::from_raw(2), FaceIdx::from_raw(10)),
        );

        journal.process_evolution(FeatureId(3), &bool_evo);

        // The two split children should have distinct EntityRefs.
        let child_a = journal.entity_for_face(FaceIdx::from_raw(32)).expect("split child a");
        let child_b = journal.entity_for_face(FaceIdx::from_raw(33)).expect("split child b");
        assert_ne!(child_a, child_b, "split children must have distinct refs");

        // Both should differ from the original.
        assert_ne!(child_a, original_ref);
        assert_ne!(child_b, original_ref);

        // Resolving the original should show it was split.
        match journal.resolve_face(original_ref) {
            ResolveResult::Split { by_feature, children } => {
                assert_eq!(by_feature, FeatureId(3));
                assert_eq!(children.len(), 2);
                let child_refs: Vec<EntityRef> = children.iter().map(|(er, _)| *er).collect();
                assert!(child_refs.contains(&child_a));
                assert!(child_refs.contains(&child_b));
            }
            other => panic!("expected Split, got {:?}", other),
        }
    }

    #[test]
    fn deleted_face_resolves_as_deleted() {
        let mut journal = NamingJournal::new();

        let box_evo = box_evolution();
        journal.process_evolution(FeatureId(1), &box_evo);

        let doomed_ref = journal
            .entity_for_face(FaceIdx::from_raw(4))
            .expect("face 4 tracked");

        // Delete face 4, keep others.
        let mut del_evo = ShapeEvolution::new();
        for i in 0..6 {
            if i == 4 {
                del_evo.record_deleted_face(FaceIdx::from_raw(i));
            } else {
                del_evo.record_face(
                    FaceIdx::from_raw(i),
                    FaceOrigin::CopiedFrom(FaceIdx::from_raw(i)),
                );
            }
        }

        journal.process_evolution(FeatureId(5), &del_evo);

        match journal.resolve_face(doomed_ref) {
            ResolveResult::Deleted { by_feature } => {
                assert_eq!(by_feature, FeatureId(5));
            }
            other => panic!("expected Deleted, got {:?}", other),
        }
    }

    #[test]
    fn replay_produces_same_entity_refs() {
        // First run.
        let mut journal = NamingJournal::new();
        let box_evo = box_evolution();
        journal.process_evolution(FeatureId(1), &box_evo);

        let refs_first: Vec<EntityRef> = (0..6)
            .map(|i| journal.entity_for_face(FaceIdx::from_raw(i)).unwrap())
            .collect();

        // Simulate replay: clear live state, re-process same evolution.
        journal.clear();
        journal.process_evolution(FeatureId(1), &box_evo);

        let refs_second: Vec<EntityRef> = (0..6)
            .map(|i| journal.entity_for_face(FaceIdx::from_raw(i)).unwrap())
            .collect();

        assert_eq!(refs_first, refs_second, "replay must produce identical EntityRefs");
    }

    #[test]
    fn replay_with_chamfer_is_stable() {
        // Full sequence: box + chamfer.
        let mut journal = NamingJournal::new();
        let box_evo = box_evolution();
        journal.process_evolution(FeatureId(1), &box_evo);

        let mut chamfer_evo = ShapeEvolution::new();
        chamfer_evo.record_face(FaceIdx::from_raw(10), FaceOrigin::Modified(FaceIdx::from_raw(0)));
        for i in 1..6 {
            chamfer_evo.record_face(
                FaceIdx::from_raw(10 + i),
                FaceOrigin::CopiedFrom(FaceIdx::from_raw(i)),
            );
        }
        chamfer_evo.record_face(FaceIdx::from_raw(16), FaceOrigin::FromEdge(EdgeIdx::from_raw(0)));
        journal.process_evolution(FeatureId(2), &chamfer_evo);

        let ref_modified = journal.entity_for_face(FaceIdx::from_raw(10)).unwrap();
        let ref_chamfer = journal.entity_for_face(FaceIdx::from_raw(16)).unwrap();

        // Replay.
        journal.clear();
        journal.process_evolution(FeatureId(1), &box_evo);
        journal.process_evolution(FeatureId(2), &chamfer_evo);

        let ref_modified_2 = journal.entity_for_face(FaceIdx::from_raw(10)).unwrap();
        let ref_chamfer_2 = journal.entity_for_face(FaceIdx::from_raw(16)).unwrap();

        assert_eq!(ref_modified, ref_modified_2);
        assert_eq!(ref_chamfer, ref_chamfer_2);
    }

    #[test]
    fn unknown_entity_ref() {
        let journal = NamingJournal::new();
        let bogus = EntityRef::new(EntityKind::Face, 999);
        assert_eq!(journal.resolve_face(bogus), ResolveResult::Unknown);
    }

    #[test]
    fn edge_and_vertex_tracking() {
        let mut journal = NamingJournal::new();
        let box_evo = box_evolution();
        journal.process_evolution(FeatureId(1), &box_evo);

        // All 12 edges should resolve.
        for i in 0..12 {
            let idx = EdgeIdx::from_raw(i);
            let er = journal.entity_for_edge(idx).expect("edge should be tracked");
            assert_eq!(er.kind(), EntityKind::Edge);
            match journal.resolve_edge(er) {
                ResolveResult::Live(resolved) => assert_eq!(resolved, idx),
                other => panic!("expected Live, got {:?}", other),
            }
        }

        // All 8 vertices should resolve.
        for i in 0..8 {
            let idx = VertexIdx::from_raw(i);
            let er = journal.entity_for_vertex(idx).expect("vertex should be tracked");
            assert_eq!(er.kind(), EntityKind::Vertex);
            match journal.resolve_vertex(er) {
                ResolveResult::Live(resolved) => assert_eq!(resolved, idx),
                other => panic!("expected Live, got {:?}", other),
            }
        }
    }
}

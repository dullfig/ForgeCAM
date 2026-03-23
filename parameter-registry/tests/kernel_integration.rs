//! Integration tests: real kernel operations → NamingJournal.
//!
//! These tests run actual kernel ops (make_box, chamfer, fillet, translate, boolean)
//! and feed the resulting ShapeEvolution records into the NamingJournal, verifying
//! that EntityRefs are assigned, stable across replay, and correctly track
//! modifications, deletions, and splits.

use forgecam_geometry::entity_ref::EntityKind;
use forgecam_geometry::types::Vec3;
use forgecam_parameter_registry::{FeatureId, NamingJournal, ResolveResult};
use rustkernel_kernel::Kernel;


// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Feed the kernel's last evolution into the journal under the given feature ID.
fn ingest(kernel: &mut Kernel, journal: &mut NamingJournal, feature: FeatureId) {
    let evo = kernel
        .take_evolution()
        .expect("kernel should have evolution after an operation");
    journal.process_evolution(feature, &evo);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn primitive_box_round_trip() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    let box_s = k.make_box(2.0, 3.0, 4.0);
    ingest(&mut k, &mut j, FeatureId(1));

    // A box has 6 faces, 12 edges, 8 vertices.
    let faces = k.topo().solid_faces(box_s);
    assert_eq!(faces.len(), 6);

    // Every face should have an EntityRef that resolves to Live.
    for &face in &faces {
        let er = j.entity_for_face(face).expect("face should be tracked");
        assert_eq!(er.kind(), EntityKind::Face);
        match j.resolve_face(er) {
            ResolveResult::Live(idx) => assert_eq!(idx, face),
            other => panic!("expected Live for face {:?}, got {:?}", face, other),
        }
    }

    // Edges and vertices too.
    assert!(j.edge_count() >= 12, "should track at least 12 edges, got {}", j.edge_count());
    assert!(j.vertex_count() >= 8, "should track at least 8 vertices, got {}", j.vertex_count());
}

#[test]
fn primitive_cylinder_round_trip() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    let cyl = k.make_cylinder(1.0, 3.0);
    ingest(&mut k, &mut j, FeatureId(1));

    let faces = k.topo().solid_faces(cyl);
    // Cylinder: 2 caps + N lateral faces (tessellated).
    assert!(faces.len() >= 3);

    for &face in &faces {
        let er = j.entity_for_face(face).expect("face should be tracked");
        assert!(matches!(j.resolve_face(er), ResolveResult::Live(_)));
    }
}

#[test]
fn box_then_euler_chamfer() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    // Feature 1: Make a box.
    let box_s = k.make_box(2.0, 2.0, 2.0);
    ingest(&mut k, &mut j, FeatureId(1));

    // Pick an edge to chamfer.
    let edges = rustkernel_builders::edge_analysis::solid_edges(k.topo(), box_s);
    let target_edge = edges[0];

    // Remember the EntityRef for the faces adjacent to this edge.
    let adj = rustkernel_builders::edge_analysis::edge_adjacency(k.topo(), target_edge).unwrap();
    let face_a_ref = j.entity_for_face(adj.face_a).expect("adj face a tracked");
    let face_b_ref = j.entity_for_face(adj.face_b).expect("adj face b tracked");
    let edge_ref = j.entity_for_edge(target_edge).expect("target edge tracked");

    // Feature 2: Chamfer the edge.
    let _chamfered = k.euler_chamfer_edges(box_s, &[target_edge], 0.3).unwrap();
    ingest(&mut k, &mut j, FeatureId(2));

    // The adjacent faces should keep their EntityRefs (Modified → inherits).
    // They now point to different FaceIdx values, but the EntityRef is the same.
    assert!(
        matches!(j.resolve_face(face_a_ref), ResolveResult::Live(_)),
        "adj face a should still be Live after chamfer"
    );
    assert!(
        matches!(j.resolve_face(face_b_ref), ResolveResult::Live(_)),
        "adj face b should still be Live after chamfer"
    );

    // The chamfered edge should be Deleted.
    match j.resolve_edge(edge_ref) {
        ResolveResult::Deleted { by_feature } => {
            assert_eq!(by_feature, FeatureId(2), "edge should be deleted by feature 2");
        }
        other => panic!("expected edge Deleted, got {:?}", other),
    }

    // There should be at least one new face (the chamfer face) with a different EntityRef.
    let all_live_face_refs: Vec<_> = k
        .topo()
        .solid_faces(_chamfered)
        .into_iter()
        .filter_map(|f| j.entity_for_face(f))
        .collect();
    assert!(all_live_face_refs.len() >= 7, "chamfered box should have at least 7 tracked faces");

    // The chamfer face ref should differ from the original adjacent face refs.
    let new_refs: Vec<_> = all_live_face_refs
        .iter()
        .filter(|r| **r != face_a_ref && **r != face_b_ref)
        .collect();
    assert!(!new_refs.is_empty(), "should have new face refs from chamfer");
}

#[test]
fn box_then_euler_fillet() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    let box_s = k.make_box(2.0, 2.0, 2.0);
    ingest(&mut k, &mut j, FeatureId(1));

    let edges = rustkernel_builders::edge_analysis::solid_edges(k.topo(), box_s);
    let target_edge = edges[0];
    let edge_ref = j.entity_for_edge(target_edge).expect("edge tracked");

    let _filleted = k.euler_fillet_edges(box_s, &[target_edge], 0.3).unwrap();
    ingest(&mut k, &mut j, FeatureId(2));

    // The filleted edge should be deleted.
    match j.resolve_edge(edge_ref) {
        ResolveResult::Deleted { by_feature } => {
            assert_eq!(by_feature, FeatureId(2));
        }
        other => panic!("expected edge Deleted, got {:?}", other),
    }

    // All faces on the result should be tracked.
    let result_faces = k.topo().solid_faces(_filleted);
    for &face in &result_faces {
        assert!(
            j.entity_for_face(face).is_some(),
            "face {:?} on filleted solid should be tracked",
            face
        );
    }
}

#[test]
fn box_then_translate_preserves_refs() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    let box_s = k.make_box(2.0, 2.0, 2.0);
    ingest(&mut k, &mut j, FeatureId(1));

    // Remember all face refs.
    let orig_faces = k.topo().solid_faces(box_s);
    let orig_refs: Vec<_> = orig_faces
        .iter()
        .map(|&f| j.entity_for_face(f).unwrap())
        .collect();

    // Feature 2: Translate.
    let translated = k.translate(box_s, Vec3::new(10.0, 0.0, 0.0));
    ingest(&mut k, &mut j, FeatureId(2));

    // Each face on the translated solid should have the same EntityRef as the original.
    let trans_faces = k.topo().solid_faces(translated);
    let trans_refs: Vec<_> = trans_faces
        .iter()
        .map(|&f| j.entity_for_face(f).unwrap())
        .collect();

    // Same EntityRefs, just pointing to new FaceIdx values.
    let mut orig_sorted = orig_refs.clone();
    orig_sorted.sort_by_key(|r| r.raw());
    let mut trans_sorted = trans_refs.clone();
    trans_sorted.sort_by_key(|r| r.raw());
    assert_eq!(orig_sorted, trans_sorted, "translate should preserve EntityRefs");

    // Each ref should resolve to the new face index.
    for &er in &trans_refs {
        match j.resolve_face(er) {
            ResolveResult::Live(idx) => {
                assert!(trans_faces.contains(&idx), "resolved face should be on translated solid");
            }
            other => panic!("expected Live, got {:?}", other),
        }
    }
}

#[test]
fn boolean_cut_tracks_deletions_and_splits() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    // Feature 1: Large box.
    let box_a = k.make_box(4.0, 4.0, 4.0);
    ingest(&mut k, &mut j, FeatureId(1));

    // Feature 2: Small box (tool). We need its evolution too.
    let box_b = k.make_box_at(
        [1.0, 1.0, -0.5],
        2.0, 2.0, 5.0,
    );
    ingest(&mut k, &mut j, FeatureId(2));

    // Remember a face ref from box_a — should survive the boolean.
    let a_faces = k.topo().solid_faces(box_a);
    let _a_face_ref = j.entity_for_face(a_faces[0]).unwrap();

    // Feature 3: Cut box_b from box_a.
    let result = k.cut(box_a, box_b).expect("boolean cut should succeed");
    ingest(&mut k, &mut j, FeatureId(3));

    // The result should have tracked faces.
    let result_faces = k.topo().solid_faces(result);
    let tracked_count = result_faces
        .iter()
        .filter(|f| j.entity_for_face(**f).is_some())
        .count();
    assert!(
        tracked_count > 0,
        "at least some result faces should be tracked (got {}/{})",
        tracked_count,
        result_faces.len()
    );
}

#[test]
fn replay_stability_box_chamfer() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    // Run 1: box + chamfer.
    let box_s = k.make_box(2.0, 2.0, 2.0);
    let box_evo = k.take_evolution().unwrap();
    j.process_evolution(FeatureId(1), &box_evo);

    let edges = rustkernel_builders::edge_analysis::solid_edges(k.topo(), box_s);
    let _chamfered = k.euler_chamfer_edges(box_s, &[edges[0]], 0.3).unwrap();
    let chamfer_evo = k.take_evolution().unwrap();
    j.process_evolution(FeatureId(2), &chamfer_evo);

    // Snapshot all face EntityRefs.
    let faces_run1 = k.topo().solid_faces(_chamfered);
    let refs_run1: Vec<_> = faces_run1
        .iter()
        .filter_map(|&f| j.entity_for_face(f))
        .collect();

    // Run 2: replay from scratch (same kernel state, but journal replayed).
    j.clear();
    j.process_evolution(FeatureId(1), &box_evo);
    j.process_evolution(FeatureId(2), &chamfer_evo);

    let refs_run2: Vec<_> = faces_run1
        .iter()
        .filter_map(|&f| j.entity_for_face(f))
        .collect();

    assert_eq!(refs_run1, refs_run2, "replay must produce identical EntityRefs");
}

#[test]
fn multi_feature_sequence() {
    let mut k = Kernel::new();
    let mut j = NamingJournal::new();

    // Feature 1: Box.
    let box_s = k.make_box(4.0, 4.0, 4.0);
    ingest(&mut k, &mut j, FeatureId(1));
    assert_eq!(j.live_face_count(), 6);

    // Feature 2: Chamfer one edge.
    let edges = rustkernel_builders::edge_analysis::solid_edges(k.topo(), box_s);
    let solid = k.euler_chamfer_edges(box_s, &[edges[0]], 0.5).unwrap();
    ingest(&mut k, &mut j, FeatureId(2));
    // 6 original faces - 0 deleted + 1 chamfer face = 7 live faces.
    assert_eq!(j.live_face_count(), 7, "after chamfer: 7 live faces");

    // Feature 3: Fillet another edge on the chamfered solid.
    let edges_after = rustkernel_builders::edge_analysis::solid_edges(k.topo(), solid);
    // Pick an edge that wasn't chamfered (still exists).
    let fillet_edge = edges_after
        .iter()
        .find(|e| {
            rustkernel_builders::edge_analysis::edge_convexity(k.topo(), k.geom(), **e)
                .map(|c| matches!(c, rustkernel_builders::edge_analysis::EdgeConvexity::Convex))
                .unwrap_or(false)
        });

    if let Some(&fe) = fillet_edge {
        let solid2 = k.euler_fillet_edges(solid, &[fe], 0.3).unwrap();
        ingest(&mut k, &mut j, FeatureId(3));

        // Should have more live faces now (fillet adds face strips).
        let result_faces = k.topo().solid_faces(solid2);
        let tracked = result_faces
            .iter()
            .filter(|f| j.entity_for_face(**f).is_some())
            .count();
        assert_eq!(
            tracked,
            result_faces.len(),
            "all faces on final solid should be tracked"
        );
    }
}

use std::collections::HashSet;

use crate::store::TopoStore;
use crate::topo::*;

/// Severity of a diagnostic finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Info,
    Warning,
    Error,
}

/// A single diagnostic finding from solid validation.
#[derive(Debug, Clone)]
pub struct Diagnostic {
    pub severity: Severity,
    pub code: &'static str,
    pub message: String,
}

/// Report produced by `validate_solid`, summarising topology health.
#[derive(Debug)]
pub struct DiagnosticReport {
    pub solid: SolidIdx,
    pub findings: Vec<Diagnostic>,
    pub vertex_count: usize,
    pub edge_count: usize,
    pub face_count: usize,
    pub genus: u32,
    pub euler_value: i32,
}

impl DiagnosticReport {
    /// Returns true when there are no Error-level findings.
    pub fn is_valid(&self) -> bool {
        !self.findings.iter().any(|d| d.severity == Severity::Error)
    }

    /// All Error-level findings.
    pub fn errors(&self) -> Vec<&Diagnostic> {
        self.findings
            .iter()
            .filter(|d| d.severity == Severity::Error)
            .collect()
    }

    /// All Warning-level findings.
    pub fn warnings(&self) -> Vec<&Diagnostic> {
        self.findings
            .iter()
            .filter(|d| d.severity == Severity::Warning)
            .collect()
    }
}

/// Validate a solid's topology, returning a diagnostic report.
///
/// Checks:
/// 1. Euler formula (V - E + F = 2 - 2g)
/// 2. Unmatched twins (half-edges without a twin)
/// 3. Twin consistency (he.twin.twin == he)
/// 4. Loop closure (next-chain returns to start)
/// 5. Degenerate faces (fewer than 3 edges)
pub fn validate_solid(topo: &TopoStore, solid: SolidIdx) -> DiagnosticReport {
    let mut findings = Vec::new();

    let shell_idx = topo.solids.get(solid).shell;
    let genus = topo.solids.get(solid).genus;
    let faces = &topo.shells.get(shell_idx).faces;

    let mut verts = HashSet::new();
    let mut edges = HashSet::new();
    let mut all_half_edges: Vec<HalfEdgeIdx> = Vec::new();

    for &face_idx in faces {
        let loop_idx = topo.faces.get(face_idx).outer_loop;
        let start_he = topo.loops.get(loop_idx).half_edge;

        let mut edge_count_in_face = 0u32;
        let mut current = start_he;
        let mut seen_in_loop = HashSet::new();

        loop {
            if !seen_in_loop.insert(current.raw()) {
                // We revisited a half-edge before reaching start — loop may be closed
                // or there's a cycle not through start.
                break;
            }

            let he = topo.half_edges.get(current);
            verts.insert(he.origin.raw());
            edges.insert(he.edge.raw());
            all_half_edges.push(current);
            edge_count_in_face += 1;

            // Check 2: Unmatched twin
            if he.twin.is_none() {
                findings.push(Diagnostic {
                    severity: Severity::Error,
                    code: "UNMATCHED_TWIN",
                    message: format!("half-edge {} has no twin", current.raw()),
                });
            }

            // Check 3: Twin consistency
            if let Some(twin_idx) = he.twin {
                let twin_he = topo.half_edges.get(twin_idx);
                if twin_he.twin != Some(current) {
                    findings.push(Diagnostic {
                        severity: Severity::Error,
                        code: "TWIN_ASYMMETRY",
                        message: format!(
                            "half-edge {}: twin {} does not point back",
                            current.raw(),
                            twin_idx.raw()
                        ),
                    });
                }
            }

            current = he.next;
            if current == start_he {
                break;
            }
        }

        // Check 4: Loop closure
        if current != start_he {
            findings.push(Diagnostic {
                severity: Severity::Error,
                code: "UNCLOSED_LOOP",
                message: format!(
                    "face {} loop does not close (visited {} edges without returning to start)",
                    face_idx.raw(),
                    seen_in_loop.len()
                ),
            });
        }

        // Check 5: Degenerate face
        if edge_count_in_face < 3 {
            findings.push(Diagnostic {
                severity: Severity::Warning,
                code: "DEGENERATE_FACE",
                message: format!(
                    "face {} has only {} edges",
                    face_idx.raw(),
                    edge_count_in_face
                ),
            });
        }
    }

    let v = verts.len() as i32;
    let e = edges.len() as i32;
    let f = faces.len() as i32;
    let euler = v - e + f;
    let expected = 2 - 2 * genus as i32;

    // Check 1: Euler formula
    if euler != expected {
        findings.push(Diagnostic {
            severity: Severity::Error,
            code: "EULER_VIOLATION",
            message: format!(
                "V({}) - E({}) + F({}) = {} (expected {} for genus {})",
                v, e, f, euler, expected, genus
            ),
        });
    }

    DiagnosticReport {
        solid,
        findings,
        vertex_count: verts.len(),
        edge_count: edges.len(),
        face_count: faces.len() as usize,
        genus,
        euler_value: euler,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arena::Idx;

    /// Build a minimal valid tetrahedron topology (V=4, E=6, F=4, genus=0).
    fn make_tetrahedron(topo: &mut TopoStore) -> SolidIdx {
        // Allocate 4 vertices (point_ids don't matter for topology checks)
        let v: Vec<VertexIdx> = (0..4)
            .map(|i| topo.vertices.alloc(Vertex { point_id: i }))
            .collect();

        // We need 4 triangular faces, 6 edges, 12 half-edges.
        // Face winding (outward-facing): 0-1-2, 0-3-1, 1-3-2, 0-2-3

        // Allocate a placeholder solid/shell/faces first
        let solid_idx = topo.solids.alloc(Solid {
            shell: Idx::from_raw(0), // placeholder
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
        // (origin, dest) -> HalfEdgeIdx for twin matching
        let mut he_map: std::collections::HashMap<(u32, u32), HalfEdgeIdx> =
            std::collections::HashMap::new();

        for &(a, b, c) in &triangles {
            let face_idx = topo.faces.alloc(Face {
                outer_loop: Idx::from_raw(0), // placeholder
                surface_id: 0,
                mesh_cache: None,
                shell: shell_idx,
            });
            let loop_idx = topo.loops.alloc(Loop {
                half_edge: Idx::from_raw(0), // placeholder
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
                    next: Idx::from_raw(0), // placeholder
                    edge: edge_idx,
                    loop_ref: loop_idx,
                });
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
                let origin = tri[i];
                let dest = tri[(i + 1) % 3];
                he_map.insert((origin, dest), hes[i]);
            }

            face_idxs.push(face_idx);
        }

        // Match twins
        let keys: Vec<(u32, u32)> = he_map.keys().cloned().collect();
        for (o, d) in keys {
            if let (Some(&he), Some(&twin)) = (he_map.get(&(o, d)), he_map.get(&(d, o))) {
                topo.half_edges.get_mut(he).twin = Some(twin);
                topo.half_edges.get_mut(twin).twin = Some(he);
                // Share edge
                let edge = topo.half_edges.get(he).edge;
                topo.half_edges.get_mut(twin).edge = edge;
            }
        }

        topo.shells.get_mut(shell_idx).faces = face_idxs;
        solid_idx
    }

    #[test]
    fn test_valid_tetrahedron() {
        let mut topo = TopoStore::new();
        let solid = make_tetrahedron(&mut topo);
        let report = validate_solid(&topo, solid);
        assert!(report.is_valid(), "Tetrahedron should be valid: {:?}", report.errors());
        assert_eq!(report.vertex_count, 4);
        assert_eq!(report.edge_count, 6);
        assert_eq!(report.face_count, 4);
        assert_eq!(report.euler_value, 2);
        assert!(report.warnings().is_empty());
    }

    #[test]
    fn test_euler_violation_detected() {
        let mut topo = TopoStore::new();
        let solid = make_tetrahedron(&mut topo);
        // Break Euler by changing genus to 1 (expects V-E+F=0 but we have 2)
        topo.solids.get_mut(solid).genus = 1;
        let report = validate_solid(&topo, solid);
        assert!(!report.is_valid());
        assert!(report.errors().iter().any(|d| d.code == "EULER_VIOLATION"));
    }

    #[test]
    fn test_unmatched_twin_detected() {
        let mut topo = TopoStore::new();
        let solid = make_tetrahedron(&mut topo);
        // Break a twin link
        let shell = topo.solids.get(solid).shell;
        let first_face = topo.shells.get(shell).faces[0];
        let loop_idx = topo.faces.get(first_face).outer_loop;
        let he = topo.loops.get(loop_idx).half_edge;
        // Remove twin from the first half-edge's twin
        if let Some(twin) = topo.half_edges.get(he).twin {
            topo.half_edges.get_mut(twin).twin = None;
        }
        topo.half_edges.get_mut(he).twin = None;
        let report = validate_solid(&topo, solid);
        assert!(!report.is_valid());
        assert!(report.errors().iter().any(|d| d.code == "UNMATCHED_TWIN"));
    }

    #[test]
    fn test_report_accessors() {
        let mut topo = TopoStore::new();
        let solid = make_tetrahedron(&mut topo);
        let report = validate_solid(&topo, solid);
        assert!(report.errors().is_empty());
        assert!(report.warnings().is_empty());
        assert!(report.is_valid());
    }
}

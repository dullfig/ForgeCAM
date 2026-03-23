use crate::arena::Arena;
use crate::topo::*;
use serde::{Serialize, Deserialize};

/// Central store holding all topological entities in contiguous arenas.
#[derive(Serialize, Deserialize)]
pub struct TopoStore {
    pub vertices: Arena<Vertex>,
    pub half_edges: Arena<HalfEdge>,
    pub edges: Arena<Edge>,
    pub loops: Arena<Loop>,
    pub faces: Arena<Face>,
    pub shells: Arena<Shell>,
    pub solids: Arena<Solid>,
}

impl TopoStore {
    pub fn new() -> Self {
        Self {
            vertices: Arena::new(),
            half_edges: Arena::new(),
            edges: Arena::new(),
            loops: Arena::new(),
            faces: Arena::new(),
            shells: Arena::new(),
            solids: Arena::new(),
        }
    }
}

impl TopoStore {
    /// All faces across all shells of a solid.
    pub fn solid_faces(&self, solid: SolidIdx) -> Vec<FaceIdx> {
        self.solids.get(solid).shells.iter()
            .flat_map(|&sh| self.shells.get(sh).faces.iter().copied())
            .collect()
    }
}

impl Default for TopoStore {
    fn default() -> Self {
        Self::new()
    }
}

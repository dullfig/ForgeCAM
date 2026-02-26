use crate::arena::Arena;
use crate::topo::*;

/// Central store holding all topological entities in contiguous arenas.
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

impl Default for TopoStore {
    fn default() -> Self {
        Self::new()
    }
}

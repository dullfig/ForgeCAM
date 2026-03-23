use crate::arena::Idx;
use crate::mesh_cache::FaceMesh;
use serde::{Serialize, Deserialize};

pub type VertexIdx = Idx<Vertex>;
pub type HalfEdgeIdx = Idx<HalfEdge>;
pub type EdgeIdx = Idx<Edge>;
pub type LoopIdx = Idx<Loop>;
pub type FaceIdx = Idx<Face>;
pub type ShellIdx = Idx<Shell>;
pub type SolidIdx = Idx<Solid>;

#[derive(Serialize, Deserialize)]
pub struct Vertex {
    /// Index into the geometry store for this vertex's point.
    pub point_id: u32,
}

#[derive(Serialize, Deserialize)]
pub struct HalfEdge {
    pub origin: VertexIdx,
    pub twin: Option<HalfEdgeIdx>,
    pub next: HalfEdgeIdx,
    pub edge: EdgeIdx,
    pub loop_ref: LoopIdx,
}

#[derive(Serialize, Deserialize)]
pub struct Edge {
    pub half_edges: [HalfEdgeIdx; 2],
    /// Index into the geometry store for this edge's curve.
    pub curve_id: u32,
}

#[derive(Serialize, Deserialize)]
pub struct Loop {
    /// Any half-edge in this loop (entry point for traversal).
    pub half_edge: HalfEdgeIdx,
    pub face: FaceIdx,
}

#[derive(Serialize, Deserialize)]
pub struct Face {
    /// The outer loop of this face.
    pub outer_loop: LoopIdx,
    /// Index into the geometry store for this face's surface.
    pub surface_id: u32,
    /// Cached tessellation mesh (dual representation). Skipped during serialization.
    #[serde(skip)]
    pub mesh_cache: Option<FaceMesh>,
    pub shell: ShellIdx,
}

#[derive(Serialize, Deserialize)]
pub struct Shell {
    pub faces: Vec<FaceIdx>,
    pub solid: SolidIdx,
}

#[derive(Serialize, Deserialize)]
pub struct Solid {
    /// Shells composing this solid. `shells[0]` is always the outer shell.
    /// Additional shells represent internal voids.
    pub shells: Vec<ShellIdx>,
    /// Genus of the solid (0 for sphere-like, 1 for torus-like).
    /// Used for Euler validation: V - E + F = 2 - 2g.
    pub genus: u32,
}

impl Solid {
    /// Returns the outer shell (always `shells[0]`).
    pub fn outer_shell(&self) -> ShellIdx {
        self.shells[0]
    }
}

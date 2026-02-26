use crate::arena::Idx;
use crate::mesh_cache::FaceMesh;

pub type VertexIdx = Idx<Vertex>;
pub type HalfEdgeIdx = Idx<HalfEdge>;
pub type EdgeIdx = Idx<Edge>;
pub type LoopIdx = Idx<Loop>;
pub type FaceIdx = Idx<Face>;
pub type ShellIdx = Idx<Shell>;
pub type SolidIdx = Idx<Solid>;

pub struct Vertex {
    /// Index into the geometry store for this vertex's point.
    pub point_id: u32,
}

pub struct HalfEdge {
    pub origin: VertexIdx,
    pub twin: Option<HalfEdgeIdx>,
    pub next: HalfEdgeIdx,
    pub edge: EdgeIdx,
    pub loop_ref: LoopIdx,
}

pub struct Edge {
    pub half_edges: [HalfEdgeIdx; 2],
    /// Index into the geometry store for this edge's curve.
    pub curve_id: u32,
}

pub struct Loop {
    /// Any half-edge in this loop (entry point for traversal).
    pub half_edge: HalfEdgeIdx,
    pub face: FaceIdx,
}

pub struct Face {
    /// The outer loop of this face.
    pub outer_loop: LoopIdx,
    /// Index into the geometry store for this face's surface.
    pub surface_id: u32,
    /// Cached tessellation mesh (dual representation).
    pub mesh_cache: Option<FaceMesh>,
    pub shell: ShellIdx,
}

pub struct Shell {
    pub faces: Vec<FaceIdx>,
    pub solid: SolidIdx,
}

pub struct Solid {
    pub shell: ShellIdx,
}

use rustkernel_math::{Point3, Vec3};

/// Cached tessellation mesh for a single face.
/// Stores both rendering data and UV parameters for exact-geometry projection.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct FaceMesh {
    /// Vertex positions in 3D world space.
    pub positions: Vec<Point3>,
    /// Per-vertex normals.
    pub normals: Vec<Vec3>,
    /// Triangle indices (every 3 consecutive u32s form one triangle).
    pub indices: Vec<u32>,
    /// UV parameters on the face's surface (for projecting mesh hits back to exact geometry).
    pub uvs: Vec<[f64; 2]>,
}

impl FaceMesh {
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
            uvs: Vec::new(),
        }
    }

    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
}

impl Default for FaceMesh {
    fn default() -> Self {
        Self::new()
    }
}

// TriMesh type and mesh operations
// See API.md Section 5

use serde::{Deserialize, Serialize};

use crate::surfaces::Surface;
use crate::types::*;

// ---------------------------------------------------------------------------
// TriMesh
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriMesh {
    pub positions: Vec<Point3>,
    pub normals: Vec<Vec3>,
    pub indices: Vec<[u32; 3]>,   // triangle index triples
    pub uvs: Vec<(f64, f64)>,     // per-vertex UV coordinates
}

impl TriMesh {
    pub fn new() -> Self {
        TriMesh {
            positions: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
            uvs: Vec::new(),
        }
    }

    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    pub fn triangle_count(&self) -> usize {
        self.indices.len()
    }
}

impl Default for TriMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// tessellate_surface_grid
// ---------------------------------------------------------------------------

/// Tessellate a surface over a uniform parameter grid.
/// Evaluates (divs_u+1)*(divs_v+1) grid points via surface.eval() and surface.normal().
/// Generates 2 triangles per quad cell.
pub fn tessellate_surface_grid(
    surface: &dyn Surface,
    divs_u: usize,
    divs_v: usize,
) -> TriMesh {
    let ((u_min, u_max), (v_min, v_max)) = surface.domain();
    let nu = divs_u + 1;
    let nv = divs_v + 1;

    let mut positions = Vec::with_capacity(nu * nv);
    let mut normals = Vec::with_capacity(nu * nv);
    let mut uvs = Vec::with_capacity(nu * nv);

    for i in 0..nu {
        let u = u_min + (u_max - u_min) * i as f64 / divs_u as f64;
        for j in 0..nv {
            let v = v_min + (v_max - v_min) * j as f64 / divs_v as f64;
            positions.push(surface.eval(u, v));
            normals.push(surface.normal(u, v));
            uvs.push((u, v));
        }
    }

    let mut indices = Vec::with_capacity(divs_u * divs_v * 2);
    for i in 0..divs_u {
        for j in 0..divs_v {
            let a = (i * nv + j) as u32;
            let b = (i * nv + j + 1) as u32;
            let c = ((i + 1) * nv + j) as u32;
            let d = ((i + 1) * nv + j + 1) as u32;
            // Two triangles per quad: (a, c, b) and (b, c, d)
            indices.push([a, c, b]);
            indices.push([b, c, d]);
        }
    }

    TriMesh {
        positions,
        normals,
        indices,
        uvs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surfaces::SphereSurface;

    #[test]
    fn test_trimesh_new() {
        let mesh = TriMesh::new();
        assert_eq!(mesh.vertex_count(), 0);
        assert_eq!(mesh.triangle_count(), 0);
    }

    #[test]
    fn test_tessellate_sphere() {
        let sphere = SphereSurface {
            center: Point3::origin(),
            radius: 1.0,
        };
        let mesh = tessellate_surface_grid(&sphere, 8, 4);
        // (8+1)*(4+1) = 45 vertices
        assert_eq!(mesh.vertex_count(), 45);
        // 8*4*2 = 64 triangles
        assert_eq!(mesh.triangle_count(), 64);
        assert_eq!(mesh.uvs.len(), 45);
    }
}

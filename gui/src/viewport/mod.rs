//! Bridge between kernel tessellation and wgpu vertex buffers.

pub mod camera;

use rustkernel_kernel::Kernel;
use rustkernel_topology::topo::SolidIdx;

/// Interleaved vertex for the wgpu pipeline.
#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
}

/// GPU-ready mesh data from a kernel solid.
pub struct SolidMesh {
    pub vertices: Vec<Vertex>,
    pub indices: Vec<u32>,
}

/// Tessellate a kernel solid and produce GPU-ready interleaved vertex data.
pub fn tessellate_solid(kernel: &mut Kernel, solid: SolidIdx) -> SolidMesh {
    // Tessellate with reasonable quality for interactive display.
    let opts = rustkernel_topology::tessellate::TessellationOptions {
        min_segments: 24,
        angular_tolerance: std::f64::consts::PI / 24.0,
        chordal_tolerance: 0.005,
    };
    kernel.tessellate_solid_with_options(solid, &opts);

    // Merge per-face meshes into one interleaved vertex/index buffer.
    let shell = kernel.topo().solids.get(solid).outer_shell();
    let faces = &kernel.topo().shells.get(shell).faces.clone();
    let mut vertices: Vec<Vertex> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    for &face_idx in faces {
        let face = kernel.topo().faces.get(face_idx);
        if let Some(ref mesh) = face.mesh_cache {
            let base = vertices.len() as u32;
            for (p, n) in mesh.positions.iter().zip(mesh.normals.iter()) {
                vertices.push(Vertex {
                    position: [p.x as f32, p.y as f32, p.z as f32],
                    normal: [n.x as f32, n.y as f32, n.z as f32],
                });
            }
            for &idx in &mesh.indices {
                indices.push(base + idx);
            }
        }
    }

    SolidMesh { vertices, indices }
}

/// Build a demo scene: a box with a chamfered edge and a filleted edge.
pub fn build_demo_scene(kernel: &mut Kernel) -> Vec<(SolidIdx, SolidMesh)> {
    use rustkernel_topology::topo::EdgeIdx;

    let mut solids = Vec::new();

    // A box with one chamfered edge and one filleted edge.
    let box_s = kernel.make_box(3.0, 2.0, 1.5);

    let edges = rustkernel_builders::edge_analysis::solid_edges(kernel.topo(), box_s);
    let convex_edges: Vec<EdgeIdx> = edges
        .iter()
        .filter(|&&e| {
            matches!(
                rustkernel_builders::edge_analysis::edge_convexity(
                    kernel.topo(),
                    kernel.geom(),
                    e
                ),
                Ok(rustkernel_builders::edge_analysis::EdgeConvexity::Convex)
            )
        })
        .copied()
        .collect();

    // Chamfer one edge.
    let solid = if convex_edges.len() >= 2 {
        let chamfered = kernel
            .euler_chamfer_edges(box_s, &[convex_edges[0]], 0.2)
            .unwrap_or(box_s);

        // Fillet another edge on the chamfered solid.
        let edges_after = rustkernel_builders::edge_analysis::solid_edges(kernel.topo(), chamfered);
        let second_convex: Vec<EdgeIdx> = edges_after
            .iter()
            .filter(|&&e| {
                matches!(
                    rustkernel_builders::edge_analysis::edge_convexity(
                        kernel.topo(),
                        kernel.geom(),
                        e
                    ),
                    Ok(rustkernel_builders::edge_analysis::EdgeConvexity::Convex)
                )
            })
            .copied()
            .collect();

        if !second_convex.is_empty() {
            kernel
                .euler_fillet_edges(chamfered, &[second_convex[0]], 0.3)
                .unwrap_or(chamfered)
        } else {
            chamfered
        }
    } else {
        box_s
    };

    let mesh = tessellate_solid(kernel, solid);
    solids.push((solid, mesh));

    // Add a cylinder off to the side.
    let cyl = kernel.make_cylinder(0.6, 2.0);
    let cyl_moved = kernel.translate(
        cyl,
        nalgebra::Vector3::new(4.0, 0.0, 0.0),
    );
    let cyl_mesh = tessellate_solid(kernel, cyl_moved);
    solids.push((cyl_moved, cyl_mesh));

    solids
}

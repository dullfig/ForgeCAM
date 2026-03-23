use std::io::Write as IoWrite;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use rustkernel_geom::AnalyticalGeomStore;

use crate::Kernel;

// ── STL / OBJ export ──

/// Error type for export operations.
#[derive(Debug)]
pub enum ExportError {
    Io(std::io::Error),
    NoMesh(String),
}

impl std::fmt::Display for ExportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExportError::Io(e) => write!(f, "IO error: {e}"),
            ExportError::NoMesh(msg) => write!(f, "no mesh data: {msg}"),
        }
    }
}

impl std::error::Error for ExportError {}

impl From<std::io::Error> for ExportError {
    fn from(e: std::io::Error) -> Self { ExportError::Io(e) }
}

/// Write ASCII STL for a solid to any writer. Tessellates if mesh_cache is missing.
pub fn export_stl_ascii<W: IoWrite>(
    topo: &mut TopoStore,
    geom: &AnalyticalGeomStore,
    solid: SolidIdx,
    name: &str,
    writer: &mut W,
) -> Result<(), ExportError> {
    use rustkernel_topology::tessellate::tessellate_shell;

    for &sh in &topo.solids.get(solid).shells.clone() {
        tessellate_shell(topo, sh, geom);
    }

    writeln!(writer, "solid {name}")?;
    let faces: Vec<FaceIdx> = topo.solid_faces(solid);
    for face_idx in faces {
        if let Some(mesh) = &topo.faces.get(face_idx).mesh_cache {
            for tri in 0..mesh.triangle_count() {
                let i0 = mesh.indices[tri * 3] as usize;
                let i1 = mesh.indices[tri * 3 + 1] as usize;
                let i2 = mesh.indices[tri * 3 + 2] as usize;
                let n = mesh.normals[i0];
                let p0 = mesh.positions[i0];
                let p1 = mesh.positions[i1];
                let p2 = mesh.positions[i2];
                writeln!(writer, "  facet normal {} {} {}", n.x, n.y, n.z)?;
                writeln!(writer, "    outer loop")?;
                writeln!(writer, "      vertex {} {} {}", p0.x, p0.y, p0.z)?;
                writeln!(writer, "      vertex {} {} {}", p1.x, p1.y, p1.z)?;
                writeln!(writer, "      vertex {} {} {}", p2.x, p2.y, p2.z)?;
                writeln!(writer, "    endloop")?;
                writeln!(writer, "  endfacet")?;
            }
        }
    }
    writeln!(writer, "endsolid {name}")?;
    Ok(())
}

/// Write binary STL for a solid to any writer. Tessellates if mesh_cache is missing.
pub fn export_stl_binary<W: IoWrite>(
    topo: &mut TopoStore,
    geom: &AnalyticalGeomStore,
    solid: SolidIdx,
    writer: &mut W,
) -> Result<(), ExportError> {
    use rustkernel_topology::tessellate::tessellate_shell;

    for &sh in &topo.solids.get(solid).shells.clone() {
        tessellate_shell(topo, sh, geom);
    }

    // Count total triangles
    let faces: Vec<FaceIdx> = topo.solid_faces(solid);
    let total_tris: u32 = faces.iter()
        .filter_map(|&f| topo.faces.get(f).mesh_cache.as_ref())
        .map(|m| m.triangle_count() as u32)
        .sum();

    // 80-byte header (zeros)
    writer.write_all(&[0u8; 80])?;
    writer.write_all(&total_tris.to_le_bytes())?;

    for face_idx in faces {
        if let Some(mesh) = &topo.faces.get(face_idx).mesh_cache {
            for tri in 0..mesh.triangle_count() {
                let i0 = mesh.indices[tri * 3] as usize;
                let i1 = mesh.indices[tri * 3 + 1] as usize;
                let i2 = mesh.indices[tri * 3 + 2] as usize;
                let n = mesh.normals[i0];
                let p0 = mesh.positions[i0];
                let p1 = mesh.positions[i1];
                let p2 = mesh.positions[i2];
                // Normal (3 × f32)
                writer.write_all(&(n.x as f32).to_le_bytes())?;
                writer.write_all(&(n.y as f32).to_le_bytes())?;
                writer.write_all(&(n.z as f32).to_le_bytes())?;
                // Vertex 1
                writer.write_all(&(p0.x as f32).to_le_bytes())?;
                writer.write_all(&(p0.y as f32).to_le_bytes())?;
                writer.write_all(&(p0.z as f32).to_le_bytes())?;
                // Vertex 2
                writer.write_all(&(p1.x as f32).to_le_bytes())?;
                writer.write_all(&(p1.y as f32).to_le_bytes())?;
                writer.write_all(&(p1.z as f32).to_le_bytes())?;
                // Vertex 3
                writer.write_all(&(p2.x as f32).to_le_bytes())?;
                writer.write_all(&(p2.y as f32).to_le_bytes())?;
                writer.write_all(&(p2.z as f32).to_le_bytes())?;
                // Attribute byte count
                writer.write_all(&0u16.to_le_bytes())?;
            }
        }
    }
    Ok(())
}

/// Write OBJ format for a solid. Vertices are shared across faces for smaller files.
pub fn export_obj<W: IoWrite>(
    topo: &mut TopoStore,
    geom: &AnalyticalGeomStore,
    solid: SolidIdx,
    writer: &mut W,
) -> Result<(), ExportError> {
    use rustkernel_topology::tessellate::tessellate_shell;

    for &sh in &topo.solids.get(solid).shells.clone() {
        tessellate_shell(topo, sh, geom);
    }

    writeln!(writer, "# ForgeCAM OBJ export")?;

    let faces: Vec<FaceIdx> = topo.solid_faces(solid);

    // Global vertex + normal offset for OBJ (1-indexed)
    let mut vert_offset: usize = 1;
    let mut norm_offset: usize = 1;

    for face_idx in faces {
        if let Some(mesh) = &topo.faces.get(face_idx).mesh_cache {
            // Emit vertices
            for p in &mesh.positions {
                writeln!(writer, "v {} {} {}", p.x, p.y, p.z)?;
            }
            // Emit normals
            for n in &mesh.normals {
                writeln!(writer, "vn {} {} {}", n.x, n.y, n.z)?;
            }
            // Emit faces (1-indexed)
            for tri in 0..mesh.triangle_count() {
                let i0 = mesh.indices[tri * 3] as usize + vert_offset;
                let i1 = mesh.indices[tri * 3 + 1] as usize + vert_offset;
                let i2 = mesh.indices[tri * 3 + 2] as usize + vert_offset;
                let n0 = mesh.indices[tri * 3] as usize + norm_offset;
                let n1 = mesh.indices[tri * 3 + 1] as usize + norm_offset;
                let n2 = mesh.indices[tri * 3 + 2] as usize + norm_offset;
                writeln!(writer, "f {i0}//{n0} {i1}//{n1} {i2}//{n2}")?;
            }
            vert_offset += mesh.positions.len();
            norm_offset += mesh.normals.len();
        }
    }
    Ok(())
}

impl Kernel {
    /// Export a solid to ASCII STL format.
    pub fn export_stl_ascii<W: IoWrite>(
        &mut self, solid: SolidIdx, name: &str, writer: &mut W,
    ) -> Result<(), ExportError> {
        export_stl_ascii(&mut self.topo, &self.geom, solid, name, writer)
    }

    /// Export a solid to binary STL format.
    pub fn export_stl_binary<W: IoWrite>(
        &mut self, solid: SolidIdx, writer: &mut W,
    ) -> Result<(), ExportError> {
        export_stl_binary(&mut self.topo, &self.geom, solid, writer)
    }

    /// Export a solid to OBJ format.
    pub fn export_obj<W: IoWrite>(
        &mut self, solid: SolidIdx, writer: &mut W,
    ) -> Result<(), ExportError> {
        export_obj(&mut self.topo, &self.geom, solid, writer)
    }

    /// Export a solid to an ASCII STL file.
    pub fn export_stl_file(
        &mut self, solid: SolidIdx, path: &std::path::Path,
    ) -> Result<(), ExportError> {
        let mut file = std::fs::File::create(path)?;
        let name = path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("solid");
        self.export_stl_ascii(solid, name, &mut file)
    }

    /// Export a solid to an OBJ file.
    pub fn export_obj_file(
        &mut self, solid: SolidIdx, path: &std::path::Path,
    ) -> Result<(), ExportError> {
        let mut file = std::fs::File::create(path)?;
        self.export_obj(solid, &mut file)
    }
}

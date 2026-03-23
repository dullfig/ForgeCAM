use std::io::{Read as IoRead, Seek, Write as IoWrite};
use rustkernel_topology::store::TopoStore;
use rustkernel_geom::AnalyticalGeomStore;
use rustkernel_solvers::default_pipeline;
use tracing::info_span;

use crate::Kernel;

// ── .forge file format ──

/// Current format version.
/// v1 → v2: Solid.shell (single ShellIdx) changed to Solid.shells (Vec<ShellIdx>).
pub(crate) const FORGE_VERSION: u32 = 2;

/// Manifest stored as `manifest.json` inside the .forge ZIP container.
#[derive(serde::Serialize, serde::Deserialize)]
struct ForgeManifest {
    version: u32,
    creator: String,
}

/// Errors that can occur during .forge save/load.
#[derive(Debug)]
pub enum ForgeError {
    Io(std::io::Error),
    Zip(zip::result::ZipError),
    Bincode(bincode::Error),
    Json(serde_json::Error),
    UnsupportedVersion(u32),
    MissingEntry(String),
}

impl std::fmt::Display for ForgeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ForgeError::Io(e) => write!(f, "IO error: {e}"),
            ForgeError::Zip(e) => write!(f, "ZIP error: {e}"),
            ForgeError::Bincode(e) => write!(f, "bincode error: {e}"),
            ForgeError::Json(e) => write!(f, "JSON error: {e}"),
            ForgeError::UnsupportedVersion(v) => write!(f, "unsupported .forge version: {v}"),
            ForgeError::MissingEntry(name) => write!(f, "missing entry in .forge: {name}"),
        }
    }
}

impl std::error::Error for ForgeError {}

impl From<std::io::Error> for ForgeError {
    fn from(e: std::io::Error) -> Self { ForgeError::Io(e) }
}
impl From<zip::result::ZipError> for ForgeError {
    fn from(e: zip::result::ZipError) -> Self { ForgeError::Zip(e) }
}
impl From<bincode::Error> for ForgeError {
    fn from(e: bincode::Error) -> Self { ForgeError::Bincode(e) }
}
impl From<serde_json::Error> for ForgeError {
    fn from(e: serde_json::Error) -> Self { ForgeError::Json(e) }
}

impl Kernel {
    /// Save the kernel state to a `.forge` file (ZIP container).
    ///
    /// Layout:
    /// - `manifest.json` — version + creator metadata
    /// - `topology.bin`  — bincode-encoded `TopoStore`
    /// - `geometry.bin`  — bincode-encoded `AnalyticalGeomStore`
    ///
    /// Additional entries in the `addons/` directory are preserved on round-trip
    /// when using `save_forge_with_addons`.
    pub fn save_forge<W: IoWrite + Seek>(&self, writer: W) -> Result<(), ForgeError> {
        self.save_forge_with_addons(writer, &[])
    }

    /// Save with optional addon entries (each is a `(name, data)` pair stored
    /// under `addons/<name>`).
    pub fn save_forge_with_addons<W: IoWrite + Seek>(
        &self,
        writer: W,
        addons: &[(&str, &[u8])],
    ) -> Result<(), ForgeError> {
        let _span = info_span!("save_forge").entered();

        let mut zip = zip::ZipWriter::new(writer);
        let options = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated);

        // manifest.json
        let manifest = ForgeManifest {
            version: FORGE_VERSION,
            creator: format!("ForgeCAM {}", env!("CARGO_PKG_VERSION")),
        };
        zip.start_file("manifest.json", options)?;
        let manifest_json = serde_json::to_vec_pretty(&manifest)?;
        zip.write_all(&manifest_json)?;

        // topology.bin
        zip.start_file("topology.bin", options)?;
        let topo_bytes = bincode::serialize(&self.topo)?;
        zip.write_all(&topo_bytes)?;

        // geometry.bin
        zip.start_file("geometry.bin", options)?;
        let geom_bytes = bincode::serialize(&self.geom)?;
        zip.write_all(&geom_bytes)?;

        // addon sections
        for &(name, data) in addons {
            zip.start_file(format!("addons/{name}"), options)?;
            zip.write_all(data)?;
        }

        zip.finish()?;
        Ok(())
    }

    /// Load kernel state from a `.forge` file (ZIP container).
    /// Unknown entries (including `addons/`) are silently ignored.
    pub fn load_forge<R: IoRead + Seek>(reader: R) -> Result<Self, ForgeError> {
        let (kernel, _) = Self::load_forge_with_addons(reader)?;
        Ok(kernel)
    }

    /// Load kernel state and return any addon entries found.
    pub fn load_forge_with_addons<R: IoRead + Seek>(
        reader: R,
    ) -> Result<(Self, Vec<(String, Vec<u8>)>), ForgeError> {
        let _span = info_span!("load_forge").entered();

        let mut archive = zip::ZipArchive::new(reader)?;

        // Read manifest
        let manifest: ForgeManifest = {
            let mut entry = archive.by_name("manifest.json")
                .map_err(|_| ForgeError::MissingEntry("manifest.json".into()))?;
            let mut buf = Vec::new();
            entry.read_to_end(&mut buf)?;
            serde_json::from_slice(&buf)?
        };

        if manifest.version > FORGE_VERSION {
            return Err(ForgeError::UnsupportedVersion(manifest.version));
        }

        // Read topology
        let topo: TopoStore = {
            let mut entry = archive.by_name("topology.bin")
                .map_err(|_| ForgeError::MissingEntry("topology.bin".into()))?;
            let mut buf = Vec::new();
            entry.read_to_end(&mut buf)?;
            bincode::deserialize(&buf)?
        };

        // Read geometry
        let geom: AnalyticalGeomStore = {
            let mut entry = archive.by_name("geometry.bin")
                .map_err(|_| ForgeError::MissingEntry("geometry.bin".into()))?;
            let mut buf = Vec::new();
            entry.read_to_end(&mut buf)?;
            bincode::deserialize(&buf)?
        };

        // Collect addon entries
        let mut addons = Vec::new();
        for i in 0..archive.len() {
            let mut entry = archive.by_index(i)?;
            let name = entry.name().to_string();
            if name.starts_with("addons/") && entry.is_file() {
                let addon_name = name.strip_prefix("addons/").unwrap().to_string();
                let mut buf = Vec::new();
                entry.read_to_end(&mut buf)?;
                addons.push((addon_name, buf));
            }
        }

        let kernel = Kernel {
            topo,
            geom,
            pipeline: default_pipeline(),
            last_evolution: None,
        };

        Ok((kernel, addons))
    }

    /// Save to a file path (convenience wrapper).
    pub fn save_forge_file(&self, path: &std::path::Path) -> Result<(), ForgeError> {
        let file = std::fs::File::create(path)?;
        let writer = std::io::BufWriter::new(file);
        self.save_forge(writer)
    }

    /// Load from a file path (convenience wrapper).
    pub fn load_forge_file(path: &std::path::Path) -> Result<Self, ForgeError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::load_forge(reader)
    }
}

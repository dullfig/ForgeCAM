//! forgecam-document — Document model for ForgeCAM.
//!
//! Each `Document` owns all per-part state: B-Rep solid, annotations, feature tree,
//! toolpaths, MRSEV results, and metadata.
//!
//! The `Session` holds all open documents. One is active (being edited / rendered by
//! the GUI). Others may be running background computation (feature recognition,
//! toolpath generation). No cross-document references.

pub mod metadata;
pub mod session;

use rustkernel_kernel::Kernel;
use metadata::{PartMetadata, Units};

/// The top-level document. Everything about "the current part" lives here.
///
/// Subsystems that don't exist yet (annotations, toolpaths, feature tree, mrsev)
/// will be added as their crates come online. For now the document owns the
/// kernel state and part metadata.
///
/// Serialize/Deserialize and Debug will be derived once Kernel gains those traits.
pub struct Document {
    /// Part metadata (name, units, material, stock).
    pub metadata: PartMetadata,

    /// The B-Rep solid modeling kernel state.
    pub kernel: Kernel,

    // Future fields (uncomment as crates become available):
    // pub annotations: AnnotationSet,
    // pub features: FeatureTree,
    // pub toolpaths: ToolpathSet,
    // pub mrsev: FeatureMap,
}

impl std::fmt::Debug for Document {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Document")
            .field("metadata", &self.metadata)
            .field("kernel", &"<Kernel>")
            .finish()
    }
}

impl Document {
    /// Create a new empty document with the given name and unit system.
    pub fn new(name: impl Into<String>, units: Units) -> Self {
        Self {
            metadata: PartMetadata::new(name, units),
            kernel: Kernel::new(),
        }
    }
}

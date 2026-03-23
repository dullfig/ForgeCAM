//! forgecam-mrsev: Material Removal Shape Element Volume recognition.
//!
//! Identifies machining features in the delta volume (stock minus part),
//! generates alternative decompositions as irredundant set covers, and
//! produces [`FeatureVolume`]s for toolpath consumption.
//!
//! Based on the MRSEV framework from Kramer, Gupta, and Nau (1991-1994).
//!
//! # Usage
//!
//! ```ignore
//! use forgecam_mrsev::*;
//!
//! let features = recognize_features(
//!     kernel.topo(), kernel.geom(), stock, part, &RecognizeConfig::default(),
//! )?;
//! let models = generate_models(&features, kernel.topo(), kernel.geom(), &ModelGenConfig::default())?;
//! let closures = generate_closures(&models[0], &features, kernel.topo(), kernel.geom())?;
//! let volumes = to_feature_volumes(&models[0], &features, &closures);
//! ```

pub mod access;
pub mod closure;
pub mod evaluate;
pub mod feature;
pub mod model;
pub mod recognize;
pub mod volume;

// ── Re-exports ──

pub use access::{check_accessibility, find_access_directions, AccessDirection};
pub use closure::{generate_closures, ClosureSurface};
pub use evaluate::{compute_metrics, estimate_feature_time};
pub use feature::{FeatureId, FeatureSet, RecognizedFeature};
pub use model::{
    generate_models, propose_sequence, recalculate, ModelGenConfig, ModelMetrics, MrsevModel,
    SequenceEntry,
};
pub use recognize::{recognize_features, RecognizeConfig};
pub use volume::*;

// ── Error type ──

/// Errors from the MRSEV recognition/modeling pipeline.
#[derive(Debug, thiserror::Error)]
pub enum MrsevError {
    #[error("boolean operation failed computing delta volume: {0}")]
    BooleanFailed(String),

    #[error("no features recognized in delta volume")]
    NoFeaturesFound,

    #[error("no valid MRSEV model covers the entire delta volume (covered {covered_fraction:.1}%)")]
    IncompleteCoverage { covered_fraction: f64 },

    #[error("feature {0:?} is not accessible from any direction")]
    Inaccessible(FeatureId),

    #[error("topology error: {0}")]
    Topology(String),

    #[error("geometry error: {0}")]
    Geometry(String),
}

// ── Toolpath bridge types ──

use rustkernel_math::Vec3;
use serde::{Deserialize, Serialize};

/// What toolpath consumes from mrsev. This is the bridge type.
///
/// Contains everything toolpath needs to generate cutting paths
/// without importing the kernel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureVolume {
    pub feature_type: FeatureType,
    pub volume: MrsevVolume,
    pub access_direction: Vec3,
    pub closure: Option<ClosureSurface>,
    pub depth: f64,
    pub corner_radii: Vec<f64>,
    pub draft_angle: Option<f64>,
}

/// High-level classification of a machining feature for toolpath strategy selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FeatureType {
    Pocket,
    ThroughPocket,
    Hole,
    BlindHole,
    Slot,
    Groove,
    Chamfer,
    Fillet,
    /// Slab / face mill.
    Face,
    /// Exterior profile.
    Contour,
}

/// Convert a recognized model into [`FeatureVolume`]s for toolpath.
pub fn to_feature_volumes(
    _model: &MrsevModel,
    _feature_set: &FeatureSet,
    _closures: &[ClosureSurface],
) -> Vec<FeatureVolume> {
    // TODO Phase 6: map each sequence entry to a FeatureVolume
    todo!("to_feature_volumes: not yet implemented")
}

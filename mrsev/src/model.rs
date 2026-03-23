//! MRSEV model generation and sequencing.
//!
//! Generates irredundant set covers from a feature set, proposes default
//! machining sequences, and recalculates metrics when the user reorders.

use rustkernel_geom::AnalyticalGeomStore;
use rustkernel_math::Vec3;
use rustkernel_topology::store::TopoStore;
use serde::{Deserialize, Serialize};

use crate::feature::{FeatureId, FeatureSet};
use crate::MrsevError;

/// An MRSEV Model: an irredundant set of features that covers the delta
/// volume, presented as an ordered machining sequence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MrsevModel {
    /// Ordered sequence of features. Position = machining order.
    /// User can reorder this in the GUI; call [`recalculate`] after.
    pub sequence: Vec<SequenceEntry>,
    pub metrics: ModelMetrics,
}

/// One step in the machining sequence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceEntry {
    pub feature_id: FeatureId,
    pub access_direction: Vec3,
    /// Intermediate workpiece volume remaining after this feature is removed.
    /// Recomputed by [`recalculate`].
    pub remaining_volume: f64,
}

/// Computed metrics for a model. Recalculated when user reorders.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ModelMetrics {
    pub total_volume_removed: f64,
    /// Number of unique orientations.
    pub setup_count: u32,
    pub tool_change_count: u32,
    /// Rough time estimate in seconds.
    pub estimated_time_secs: f64,
    /// 1.0 = complete coverage of delta volume.
    pub coverage_fraction: f64,
    /// Warnings, e.g. "feature X inaccessible at position 3".
    pub warnings: Vec<String>,
}

/// Configuration for model generation.
#[derive(Debug, Clone)]
pub struct ModelGenConfig {
    /// Limit number of alternative models explored.
    pub max_models: usize,
    /// Prune overly complex decompositions.
    pub max_features_per_model: usize,
}

impl Default for ModelGenConfig {
    fn default() -> Self {
        Self {
            max_models: 10,
            max_features_per_model: 50,
        }
    }
}

/// Generate candidate MRSEV models from a feature set.
///
/// Uses irredundant set cover to find valid decompositions, then
/// proposes a default sequence for each (grouped by orientation,
/// big features first). Returns multiple alternatives for user to choose from.
pub fn generate_models(
    _feature_set: &FeatureSet,
    _topo: &TopoStore,
    _geom: &AnalyticalGeomStore,
    _config: &ModelGenConfig,
) -> Result<Vec<MrsevModel>, MrsevError> {
    // TODO Phase 5: set cover algorithm (FIND_COVERS), sequence proposal
    todo!("generate_models: not yet implemented")
}

/// Recalculate metrics after user reorders the sequence.
///
/// Recomputes: intermediate workpiece states, accessibility at each step,
/// setup/tool-change counts, estimated time, and warnings for any
/// features that become inaccessible due to the new ordering.
pub fn recalculate(
    _model: &mut MrsevModel,
    _feature_set: &FeatureSet,
    _topo: &TopoStore,
    _geom: &AnalyticalGeomStore,
) -> Result<(), MrsevError> {
    // TODO Phase 5: intermediate workpiece states + accessibility warnings
    todo!("recalculate: not yet implemented")
}

/// Propose a default sequence for a set of features.
///
/// Heuristic: group by access direction (minimize setups), then order
/// large-to-small within each group (rough before finish).
pub fn propose_sequence(
    _feature_ids: &[FeatureId],
    _feature_set: &FeatureSet,
) -> Vec<SequenceEntry> {
    // TODO Phase 5: grouping + ordering heuristic
    todo!("propose_sequence: not yet implemented")
}

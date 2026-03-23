//! Metric computation for MRSEV models.
//!
//! Evaluation is metric computation rather than an optimization objective.
//! The user sees the numbers and decides. Built-in estimators feed ModelMetrics.

use crate::feature::{FeatureSet, RecognizedFeature};
use crate::model::{ModelMetrics, SequenceEntry};

/// Estimate machining time for a single feature (rough approximation).
///
/// Uses volume / MRR (material removal rate) for the feature type.
/// `mrr` is in volume units per second.
pub fn estimate_feature_time(feature: &RecognizedFeature, mrr: f64) -> f64 {
    if mrr <= 0.0 {
        return 0.0;
    }
    feature.effective_volume / mrr
}

/// Compute full [`ModelMetrics`] for a sequence.
///
/// Called by [`crate::model::recalculate`] internally; also available standalone.
///
/// * `mrr` — material removal rate (volume units / second)
/// * `setup_time` — time per setup change (seconds)
/// * `tool_change_time` — time per tool change (seconds)
pub fn compute_metrics(
    _sequence: &[SequenceEntry],
    _feature_set: &FeatureSet,
    _mrr: f64,
    _setup_time: f64,
    _tool_change_time: f64,
) -> ModelMetrics {
    // TODO Phase 5: compute setup count, tool changes, time estimates, coverage
    todo!("compute_metrics: not yet implemented")
}

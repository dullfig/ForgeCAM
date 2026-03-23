//! Closure surface generation.
//!
//! A closure surface "plugs" a feature opening to reconstruct the
//! stock volume above it. Needed by toolpath to know where material starts.

use rustkernel_geom::{AnalyticalGeomStore, SurfaceDef};
use rustkernel_topology::store::TopoStore;
use serde::{Deserialize, Serialize};

use crate::feature::FeatureId;
use crate::model::MrsevModel;
use crate::volume::Profile2D;
use crate::MrsevError;
use crate::FeatureSet;

/// A closure surface caps a feature opening to reconstruct the stock
/// volume above it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClosureSurface {
    pub feature_id: FeatureId,
    /// The cap surface definition (from kernel geometry).
    pub surface: SurfaceDef,
    /// Trim boundary of the closure surface.
    pub boundary: Profile2D,
}

/// Generate closure surfaces for all features in a model.
///
/// Each feature gets a cap surface at its entry face that represents
/// the stock surface before the feature was machined.
pub fn generate_closures(
    _model: &MrsevModel,
    _feature_set: &FeatureSet,
    _topo: &TopoStore,
    _geom: &AnalyticalGeomStore,
) -> Result<Vec<ClosureSurface>, MrsevError> {
    // TODO Phase 6: cap planes for holes, extension surfaces for pockets
    todo!("generate_closures: not yet implemented")
}

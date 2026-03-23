//! Feature recognition algorithm.
//!
//! Traverses faces of the delta volume, classifies by surface type,
//! constructs maximal primary MRSEVs, truncates to stock bounds,
//! and verifies accessibility.
//!
//! Algorithm (from Gupta/Kramer/Nau 1994):
//! 1. Compute delta volume: stock - part
//! 2. For each face of delta volume:
//!    - Concave cylindrical -> candidate hole side or round pocket wall
//!    - Convex cylindrical -> candidate through-pocket wall
//!    - Planar -> candidate pocket floor or through-pocket side
//!    - Concave conical -> candidate hole end
//! 3. Construct maximal primary MRSEV instance for each candidate
//! 4. Truncate to stock bounds
//! 5. Verify accessibility (tool can reach from outside stock)
//! 6. Return FeatureSet of all valid primary MRSEVs

use rustkernel_geom::AnalyticalGeomStore;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::SolidIdx;

use crate::feature::FeatureSet;
use crate::MrsevError;

/// Configuration for the recognizer.
#[derive(Debug, Clone)]
pub struct RecognizeConfig {
    /// Geometric tolerance for face classification.
    pub tolerance: f64,
    /// Ignore features smaller than this volume.
    pub min_feature_volume: f64,
    pub recognize_holes: bool,
    pub recognize_pockets: bool,
    pub recognize_grooves: bool,
    pub recognize_edge_cuts: bool,
}

impl Default for RecognizeConfig {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            min_feature_volume: 1e-9,
            recognize_holes: true,
            recognize_pockets: true,
            recognize_grooves: true,
            recognize_edge_cuts: true,
        }
    }
}

/// Main entry point: recognize all primary MRSEVs in the delta volume.
///
/// Takes the topology and geometry stores plus stock and part solid indices.
/// Returns a [`FeatureSet`] containing all recognized features.
///
/// # Note
/// The delta volume (stock minus part) is computed internally. The stock
/// and part solids must both exist in the provided stores.
pub fn recognize_features(
    _topo: &TopoStore,
    _geom: &AnalyticalGeomStore,
    _stock: SolidIdx,
    _part: SolidIdx,
    _config: &RecognizeConfig,
) -> Result<FeatureSet, MrsevError> {
    // TODO Phase 2+: implement face traversal, classification, MRSEV construction
    todo!("recognize_features: not yet implemented")
}

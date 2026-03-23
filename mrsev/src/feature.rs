//! Recognized feature types — output of the recognition phase.

use rustkernel_math::Vec3;
use rustkernel_topology::topo::{FaceIdx, SolidIdx};
use serde::{Deserialize, Serialize};

use crate::volume::MrsevVolume;

/// A recognized feature: an MRSEV volume + metadata from recognition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecognizedFeature {
    pub id: FeatureId,
    pub volume: MrsevVolume,
    /// B-Rep faces that generated this feature.
    pub source_faces: Vec<FaceIdx>,
    /// Direction tool enters.
    pub access_direction: Vec3,
    /// Volume of material actually removed (eff(m, S)).
    pub effective_volume: f64,
    /// Feature extends fully through stock.
    pub is_through: bool,
}

/// Opaque feature identifier, stable within a recognition session.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct FeatureId(pub u32);

/// Result of the recognition phase: all primary MRSEVs found.
#[derive(Debug, Clone)]
pub struct FeatureSet {
    pub features: Vec<RecognizedFeature>,
    pub stock_solid: SolidIdx,
    pub part_solid: SolidIdx,
}

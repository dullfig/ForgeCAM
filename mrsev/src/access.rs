//! Accessibility analysis — determine tool approach directions.

use rustkernel_geom::AnalyticalGeomStore;
use rustkernel_math::Vec3;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::SolidIdx;
use serde::{Deserialize, Serialize};

use crate::feature::RecognizedFeature;
use crate::MrsevError;

/// Direction from which a tool can reach a feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AccessDirection {
    /// Unit vector, tool approach direction.
    pub direction: Vec3,
    /// Distance from stock surface to feature entry.
    pub clearance: f64,
}

/// Check whether a feature is accessible from a given direction,
/// considering the current workpiece state.
pub fn check_accessibility(
    _feature: &RecognizedFeature,
    _workpiece: SolidIdx,
    _direction: &Vec3,
    _tool_radius: f64,
    _topo: &TopoStore,
    _geom: &AnalyticalGeomStore,
) -> Result<bool, MrsevError> {
    // TODO Phase 5: ray/cylinder casting against workpiece
    todo!("check_accessibility: not yet implemented")
}

/// Find all valid access directions for a feature.
///
/// For 3-axis: checks the 6 principal directions + feature orientation.
pub fn find_access_directions(
    _feature: &RecognizedFeature,
    _workpiece: SolidIdx,
    _tool_radius: f64,
    _topo: &TopoStore,
    _geom: &AnalyticalGeomStore,
) -> Result<Vec<AccessDirection>, MrsevError> {
    // TODO Phase 5: principal direction + feature orientation testing
    todo!("find_access_directions: not yet implemented")
}

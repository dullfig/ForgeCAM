use serde::{Deserialize, Serialize};

/// Unit system for the document.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Units {
    Millimeter,
    Inch,
}

/// Material specification for the part.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Material {
    /// Material name (e.g. "6061-T6", "304 Stainless").
    pub name: String,
    /// Density in kg/m^3 (or lb/in^3 depending on unit system — stored in SI).
    pub density: Option<f64>,
}

/// Stock definition — the raw material shape before machining.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum StockShape {
    /// Rectangular bar stock.
    Block { x: f64, y: f64, z: f64 },
    /// Round bar stock.
    Cylinder { diameter: f64, length: f64 },
    /// Stock is a specific solid in the kernel (e.g. a casting).
    FromSolid,
}

/// Part-level metadata that travels with the document.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PartMetadata {
    pub name: String,
    pub units: Units,
    pub material: Option<Material>,
    pub stock: Option<StockShape>,
}

impl PartMetadata {
    pub fn new(name: impl Into<String>, units: Units) -> Self {
        Self {
            name: name.into(),
            units,
            material: None,
            stock: None,
        }
    }
}

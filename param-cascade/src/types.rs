use serde::{Deserialize, Serialize};
use std::fmt;

/// Canonical parameter key — interned string for fast comparison.
///
/// The dictionary defines what keys exist. Toolpath types consume keys
/// from the dictionary rather than inventing their own.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ParamKey(pub String);

impl ParamKey {
    pub fn new(s: &str) -> Self {
        Self(s.to_string())
    }
}

impl fmt::Display for ParamKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl<S: Into<String>> From<S> for ParamKey {
    fn from(s: S) -> Self {
        Self(s.into())
    }
}

/// What kind of value a parameter holds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ValueType {
    F64,
    Bool,
    Int,
    Enum(Vec<String>),
    String,
}

/// A concrete parameter value.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ParamValue {
    F64(f64),
    Bool(bool),
    Int(i64),
    Enum(String),
    String(String),
}

impl ParamValue {
    /// Extract as f64, returning None if not numeric.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::F64(v) => Some(*v),
            Self::Int(v) => Some(*v as f64),
            _ => None,
        }
    }
}

/// Units for display and conversion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ParamUnit {
    Length,
    Percent,
    Speed,
    FeedRate,
    FeedPerTooth,
    SpindleRpm,
    Angle,
    Count,
    None,
}

/// Definition of a parameter in the dictionary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamDef {
    pub key: ParamKey,
    pub display_name: String,
    pub description: String,
    pub unit: ParamUnit,
    pub value_type: ValueType,
    pub default: ParamValue,
    pub range: Option<(f64, f64)>,
}

/// State of a parameter at a tree node.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParamState {
    /// No local value — cascade from parent. (Nodes in this state have no
    /// ParamEntry; inherit is represented by absence.)
    Inherit,
    /// Local value — bulk edit CAN change it.
    Override,
    /// Local value — bulk edit SKIPS it.
    Locked,
}

/// A parameter value set at a node (Override or Locked only).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamEntry {
    pub value: ParamValue,
    pub state: ParamState,
}

/// Index into the node arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NodeId(pub u32);

/// What kind of tree node.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NodeKind {
    Folder,
    Operation(String),
}

/// Declares which dictionary params a toolpath type uses.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpTypeRegistration {
    pub op_type: String,
    pub display_name: String,
    pub required_params: Vec<ParamKey>,
    pub optional_params: Vec<ParamKey>,
}

/// Result of resolving a parameter through the cascade.
#[derive(Debug, Clone)]
pub struct ResolvedParam {
    pub value: ParamValue,
    pub source: NodeId,
    pub source_name: String,
    pub is_locked: bool,
    pub is_default: bool,
}

/// Result of a bulk edit operation.
#[derive(Debug, Clone)]
pub struct BulkEditResult {
    pub updated: Vec<NodeId>,
    pub skipped_locked: Vec<NodeId>,
}

/// Validation error for an operation node.
#[derive(Debug, Clone)]
pub enum ValidationError {
    MissingRequired { key: ParamKey, op_type: String },
    TypeMismatch { key: ParamKey, expected: ValueType, got: ParamValue },
    OutOfRange { key: ParamKey, value: f64, min: f64, max: f64 },
    UnknownParam { key: ParamKey },
}

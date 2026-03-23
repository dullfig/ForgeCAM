//! Persistent entity IDs and naming journal for ForgeCAM.
//!
//! This crate bridges between the kernel's transient arena indices (`FaceIdx`, `EdgeIdx`,
//! `VertexIdx`) and the stable [`EntityRef`] handles that annotations, MRSEV, and
//! the feature tree use to reference model entities.
//!
//! The central type is [`NamingJournal`], which processes [`ShapeEvolution`] records
//! from kernel operations and maintains a bidirectional mapping between stable
//! `EntityRef` handles and live arena indices.

pub mod naming;

pub use naming::{FeatureId, NamingJournal, NamingKey, ResolveResult};

//! # param-cascade
//!
//! Hierarchical machining parameter dictionary with folder inheritance and lock/protect.
//!
//! Every toolpath operation binds to canonical parameter keys from a shared dictionary.
//! Parameters cascade through a folder tree: operation → parent folder → global default.
//! Any parameter can be locked to prevent bulk edits from changing it.

mod types;
mod tree;
mod cascade;
mod bulk;
mod query;
mod registry;

pub use types::*;
pub use tree::{OpTree, OpNode};

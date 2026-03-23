//! Persistent entity reference — an opaque handle that survives kernel operations.
//!
//! `EntityRef` lives in the geometry crate so that every downstream crate (annotations,
//! toolpath, mrsev) can reference model entities without importing the kernel.
//!
//! The top 2 bits encode the entity kind (face/edge/vertex). The remaining 62 bits
//! are a monotonically increasing counter managed by the `NamingJournal` in
//! `parameter-registry`.

use serde::{Deserialize, Serialize};
use std::fmt;

/// What kind of topological entity an `EntityRef` points to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u8)]
pub enum EntityKind {
    Face = 0,
    Edge = 1,
    Vertex = 2,
}

/// Opaque, persistent handle to a topological entity.
///
/// Encodes [`EntityKind`] in the top 2 bits and a unique counter in the lower 62.
/// Two `EntityRef` values are equal iff they were allocated by the same journal
/// entry — meaning they identify the same logical entity across feature-tree replays.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct EntityRef(u64);

const KIND_SHIFT: u32 = 62;
const COUNTER_MASK: u64 = (1u64 << KIND_SHIFT) - 1;

impl EntityRef {
    /// Create a new `EntityRef` from a kind and counter value.
    ///
    /// # Panics
    /// Panics if `counter` exceeds 62-bit range.
    pub fn new(kind: EntityKind, counter: u64) -> Self {
        assert!(
            counter <= COUNTER_MASK,
            "EntityRef counter overflow: {counter}"
        );
        let bits = (kind as u64) << KIND_SHIFT | counter;
        Self(bits)
    }

    /// The kind of entity this reference points to.
    pub fn kind(self) -> EntityKind {
        match self.0 >> KIND_SHIFT {
            0 => EntityKind::Face,
            1 => EntityKind::Edge,
            2 => EntityKind::Vertex,
            _ => unreachable!("invalid EntityKind bits"),
        }
    }

    /// The raw counter value (lower 62 bits). Unique within a journal.
    pub fn counter(self) -> u64 {
        self.0 & COUNTER_MASK
    }

    /// The raw u64 representation.
    pub fn raw(self) -> u64 {
        self.0
    }

    /// Reconstruct from a raw u64 (e.g., after deserialization).
    pub fn from_raw(raw: u64) -> Self {
        Self(raw)
    }
}

impl fmt::Debug for EntityRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "EntityRef({:?}:{})", self.kind(), self.counter())
    }
}

impl fmt::Display for EntityRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind_char = match self.kind() {
            EntityKind::Face => 'F',
            EntityKind::Edge => 'E',
            EntityKind::Vertex => 'V',
        };
        write!(f, "{kind_char}{}", self.counter())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_kind_and_counter() {
        for kind in [EntityKind::Face, EntityKind::Edge, EntityKind::Vertex] {
            let er = EntityRef::new(kind, 42);
            assert_eq!(er.kind(), kind);
            assert_eq!(er.counter(), 42);
        }
    }

    #[test]
    fn raw_round_trip() {
        let er = EntityRef::new(EntityKind::Edge, 12345);
        let raw = er.raw();
        let er2 = EntityRef::from_raw(raw);
        assert_eq!(er, er2);
        assert_eq!(er2.kind(), EntityKind::Edge);
        assert_eq!(er2.counter(), 12345);
    }

    #[test]
    fn display_format() {
        assert_eq!(format!("{}", EntityRef::new(EntityKind::Face, 7)), "F7");
        assert_eq!(format!("{}", EntityRef::new(EntityKind::Edge, 99)), "E99");
        assert_eq!(format!("{}", EntityRef::new(EntityKind::Vertex, 0)), "V0");
    }

    #[test]
    fn different_kinds_not_equal() {
        let a = EntityRef::new(EntityKind::Face, 1);
        let b = EntityRef::new(EntityKind::Edge, 1);
        assert_ne!(a, b);
    }

    #[test]
    fn hash_works() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(EntityRef::new(EntityKind::Face, 5));
        set.insert(EntityRef::new(EntityKind::Face, 5));
        set.insert(EntityRef::new(EntityKind::Edge, 5));
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn max_counter() {
        let max = (1u64 << 62) - 1;
        let er = EntityRef::new(EntityKind::Vertex, max);
        assert_eq!(er.counter(), max);
        assert_eq!(er.kind(), EntityKind::Vertex);
    }

    #[test]
    #[should_panic(expected = "counter overflow")]
    fn counter_overflow_panics() {
        EntityRef::new(EntityKind::Face, 1u64 << 62);
    }

    #[test]
    fn serde_round_trip() {
        let er = EntityRef::new(EntityKind::Edge, 777);
        let json = serde_json::to_string(&er).unwrap();
        let er2: EntityRef = serde_json::from_str(&json).unwrap();
        assert_eq!(er, er2);
    }
}

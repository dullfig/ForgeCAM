use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::types::*;

/// A node in the operation tree (folder or leaf operation).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpNode {
    pub id: NodeId,
    pub name: String,
    pub kind: NodeKind,
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
    pub params: HashMap<ParamKey, ParamEntry>,
}

/// The parameter dictionary + operation tree.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpTree {
    pub(crate) dictionary: HashMap<ParamKey, ParamDef>,
    pub(crate) op_types: HashMap<String, OpTypeRegistration>,
    pub(crate) nodes: Vec<Option<OpNode>>, // arena with tombstones
    pub(crate) roots: Vec<NodeId>,
    pub(crate) next_id: u32,
}

impl OpTree {
    pub fn new() -> Self {
        Self {
            dictionary: HashMap::new(),
            op_types: HashMap::new(),
            nodes: Vec::new(),
            roots: Vec::new(),
            next_id: 0,
        }
    }

    // ── Dictionary management ─────────────────────────────────────

    pub fn register_param(&mut self, def: ParamDef) {
        self.dictionary.insert(def.key.clone(), def);
    }

    pub fn param_def(&self, key: &ParamKey) -> Option<&ParamDef> {
        self.dictionary.get(key)
    }

    pub fn all_param_defs(&self) -> impl Iterator<Item = &ParamDef> {
        self.dictionary.values()
    }

    // ── Op type registration ──────────────────────────────────────

    pub fn register_op_type(&mut self, reg: OpTypeRegistration) {
        self.op_types.insert(reg.op_type.clone(), reg);
    }

    pub fn op_type_params(&self, op_type: &str) -> Option<&OpTypeRegistration> {
        self.op_types.get(op_type)
    }

    pub fn all_op_types(&self) -> impl Iterator<Item = &OpTypeRegistration> {
        self.op_types.values()
    }

    pub fn op_types_using_param(&self, key: &ParamKey) -> Vec<&str> {
        self.op_types
            .values()
            .filter(|reg| {
                reg.required_params.contains(key) || reg.optional_params.contains(key)
            })
            .map(|reg| reg.op_type.as_str())
            .collect()
    }

    // ── Tree structure ────────────────────────────────────────────

    fn alloc_node(&mut self, node: OpNode) -> NodeId {
        let id = node.id;
        let idx = id.0 as usize;
        if idx >= self.nodes.len() {
            self.nodes.resize_with(idx + 1, || None);
        }
        self.nodes[idx] = Some(node);
        id
    }

    fn next_node_id(&mut self) -> NodeId {
        let id = NodeId(self.next_id);
        self.next_id += 1;
        id
    }

    pub fn add_folder(&mut self, name: &str, parent: Option<NodeId>) -> NodeId {
        let id = self.next_node_id();
        let node = OpNode {
            id,
            name: name.to_string(),
            kind: NodeKind::Folder,
            parent,
            children: Vec::new(),
            params: HashMap::new(),
        };
        self.alloc_node(node);
        if let Some(pid) = parent {
            if let Some(Some(p)) = self.nodes.get_mut(pid.0 as usize) {
                p.children.push(id);
            }
        } else {
            self.roots.push(id);
        }
        id
    }

    pub fn add_operation(&mut self, name: &str, op_type: &str, parent: Option<NodeId>) -> NodeId {
        let id = self.next_node_id();
        let node = OpNode {
            id,
            name: name.to_string(),
            kind: NodeKind::Operation(op_type.to_string()),
            parent,
            children: Vec::new(),
            params: HashMap::new(),
        };
        self.alloc_node(node);
        if let Some(pid) = parent {
            if let Some(Some(p)) = self.nodes.get_mut(pid.0 as usize) {
                p.children.push(id);
            }
        } else {
            self.roots.push(id);
        }
        id
    }

    pub fn remove_node(&mut self, id: NodeId) {
        // Collect children first to remove recursively
        let children: Vec<NodeId> = self
            .node(id)
            .map(|n| n.children.clone())
            .unwrap_or_default();
        for child in children {
            self.remove_node(child);
        }

        // Remove from parent's children list
        if let Some(parent_id) = self.node(id).and_then(|n| n.parent) {
            if let Some(Some(parent)) = self.nodes.get_mut(parent_id.0 as usize) {
                parent.children.retain(|c| *c != id);
            }
        } else {
            self.roots.retain(|r| *r != id);
        }

        // Tombstone the slot
        if let Some(slot) = self.nodes.get_mut(id.0 as usize) {
            *slot = None;
        }
    }

    pub fn move_node(&mut self, id: NodeId, new_parent: Option<NodeId>) {
        // Detach from old parent
        let old_parent = self.node(id).and_then(|n| n.parent);
        if let Some(old_pid) = old_parent {
            if let Some(Some(p)) = self.nodes.get_mut(old_pid.0 as usize) {
                p.children.retain(|c| *c != id);
            }
        } else {
            self.roots.retain(|r| *r != id);
        }

        // Attach to new parent
        if let Some(new_pid) = new_parent {
            if let Some(Some(p)) = self.nodes.get_mut(new_pid.0 as usize) {
                p.children.push(id);
            }
        } else {
            self.roots.push(id);
        }

        // Update node's parent pointer
        if let Some(Some(node)) = self.nodes.get_mut(id.0 as usize) {
            node.parent = new_parent;
        }
    }

    pub fn clone_subtree(&mut self, id: NodeId, new_parent: Option<NodeId>) -> NodeId {
        let (name, kind, params, children) = {
            let node = match self.node(id) {
                Some(n) => n,
                None => return NodeId(u32::MAX),
            };
            (
                node.name.clone(),
                node.kind.clone(),
                node.params.clone(),
                node.children.clone(),
            )
        };

        let new_id = self.next_node_id();
        let new_node = OpNode {
            id: new_id,
            name,
            kind,
            parent: new_parent,
            children: Vec::new(),
            params,
        };
        self.alloc_node(new_node);

        if let Some(pid) = new_parent {
            if let Some(Some(p)) = self.nodes.get_mut(pid.0 as usize) {
                p.children.push(new_id);
            }
        } else {
            self.roots.push(new_id);
        }

        // Recursively clone children
        for child_id in children {
            self.clone_subtree(child_id, Some(new_id));
        }

        new_id
    }

    pub fn node(&self, id: NodeId) -> Option<&OpNode> {
        self.nodes.get(id.0 as usize).and_then(|slot| slot.as_ref())
    }

    pub fn node_mut(&mut self, id: NodeId) -> Option<&mut OpNode> {
        self.nodes
            .get_mut(id.0 as usize)
            .and_then(|slot| slot.as_mut())
    }

    pub fn children(&self, id: NodeId) -> &[NodeId] {
        self.node(id)
            .map(|n| n.children.as_slice())
            .unwrap_or(&[])
    }

    pub fn roots(&self) -> &[NodeId] {
        &self.roots
    }

    /// Walk up from a node to the root, yielding each ancestor.
    pub(crate) fn ancestors(&self, id: NodeId) -> Vec<NodeId> {
        let mut result = Vec::new();
        let mut current = self.node(id).and_then(|n| n.parent);
        while let Some(pid) = current {
            result.push(pid);
            current = self.node(pid).and_then(|n| n.parent);
        }
        result
    }

    /// Collect all leaf operation nodes under a given node (recursive).
    pub fn descendant_ops(&self, id: NodeId) -> Vec<NodeId> {
        let mut result = Vec::new();
        self.collect_ops(id, &mut result);
        result
    }

    fn collect_ops(&self, id: NodeId, out: &mut Vec<NodeId>) {
        if let Some(node) = self.node(id) {
            match &node.kind {
                NodeKind::Operation(_) => out.push(id),
                NodeKind::Folder => {
                    for &child in &node.children {
                        self.collect_ops(child, out);
                    }
                }
            }
        }
    }

    /// All leaf operation nodes in the tree.
    pub fn all_ops(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter_map(|slot| {
                slot.as_ref().and_then(|n| match &n.kind {
                    NodeKind::Operation(_) => Some(n.id),
                    _ => None,
                })
            })
            .collect()
    }

    // ── Parameter access ──────────────────────────────────────────

    pub fn set_param(&mut self, node: NodeId, key: ParamKey, value: ParamValue) {
        if let Some(n) = self.node_mut(node) {
            n.params.insert(
                key,
                ParamEntry {
                    value,
                    state: ParamState::Override,
                },
            );
        }
    }

    pub fn set_param_locked(&mut self, node: NodeId, key: ParamKey, value: ParamValue) {
        if let Some(n) = self.node_mut(node) {
            n.params.insert(
                key,
                ParamEntry {
                    value,
                    state: ParamState::Locked,
                },
            );
        }
    }

    pub fn lock_param(&mut self, node: NodeId, key: &ParamKey) {
        if let Some(n) = self.node_mut(node) {
            if let Some(entry) = n.params.get_mut(key) {
                entry.state = ParamState::Locked;
            }
        }
    }

    pub fn unlock_param(&mut self, node: NodeId, key: &ParamKey) {
        if let Some(n) = self.node_mut(node) {
            if let Some(entry) = n.params.get_mut(key) {
                entry.state = ParamState::Override;
            }
        }
    }

    pub fn clear_param(&mut self, node: NodeId, key: &ParamKey) {
        if let Some(n) = self.node_mut(node) {
            n.params.remove(key);
        }
    }

    pub fn local_param(&self, node: NodeId, key: &ParamKey) -> Option<&ParamEntry> {
        self.node(node).and_then(|n| n.params.get(key))
    }
}

impl Default for OpTree {
    fn default() -> Self {
        Self::new()
    }
}

use std::collections::HashMap;

use crate::tree::OpTree;
use crate::types::*;

impl OpTree {
    /// Resolve a parameter for a node by walking up the tree.
    ///
    /// Resolution order: node → parent → grandparent → ... → global default.
    /// Returns None only if the key isn't in the dictionary at all.
    pub fn resolve(&self, node: NodeId, key: &ParamKey) -> Option<ResolvedParam> {
        // Check the node itself
        if let Some(n) = self.node(node) {
            if let Some(entry) = n.params.get(key) {
                return Some(ResolvedParam {
                    value: entry.value.clone(),
                    source: node,
                    source_name: n.name.clone(),
                    is_locked: entry.state == ParamState::Locked,
                    is_default: false,
                });
            }
        }

        // Walk up ancestors
        for ancestor_id in self.ancestors(node) {
            if let Some(ancestor) = self.node(ancestor_id) {
                if let Some(entry) = ancestor.params.get(key) {
                    return Some(ResolvedParam {
                        value: entry.value.clone(),
                        source: ancestor_id,
                        source_name: ancestor.name.clone(),
                        is_locked: entry.state == ParamState::Locked,
                        is_default: false,
                    });
                }
            }
        }

        // Fall back to dictionary default
        self.dictionary.get(key).map(|def| ResolvedParam {
            value: def.default.clone(),
            source: node,
            source_name: "default".to_string(),
            is_locked: false,
            is_default: true,
        })
    }

    /// Resolve ALL known params for a node (full effective parameter set).
    ///
    /// Returns every param from the dictionary with its resolved value.
    pub fn resolve_all(&self, node: NodeId) -> HashMap<ParamKey, ResolvedParam> {
        let keys: Vec<ParamKey> = self.dictionary.keys().cloned().collect();
        let mut result = HashMap::with_capacity(keys.len());
        for key in keys {
            if let Some(resolved) = self.resolve(node, &key) {
                result.insert(key, resolved);
            }
        }
        result
    }

    /// Resolve only the params relevant to an operation's registered type.
    ///
    /// If the op type isn't registered, falls back to resolve_all.
    pub fn resolve_for_op_type(&self, node: NodeId) -> HashMap<ParamKey, ResolvedParam> {
        let op_type = match self.node(node) {
            Some(n) => match &n.kind {
                NodeKind::Operation(t) => t.clone(),
                _ => return self.resolve_all(node),
            },
            None => return HashMap::new(),
        };

        let keys: Vec<ParamKey> = match self.op_types.get(&op_type) {
            Some(reg) => reg
                .required_params
                .iter()
                .chain(reg.optional_params.iter())
                .cloned()
                .collect(),
            None => return self.resolve_all(node),
        };

        let mut result = HashMap::with_capacity(keys.len());
        for key in keys {
            if let Some(resolved) = self.resolve(node, &key) {
                result.insert(key, resolved);
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tree_with_dict() -> OpTree {
        let mut tree = OpTree::new();
        tree.register_param(ParamDef {
            key: ParamKey::new("stepover"),
            display_name: "Stepover".into(),
            description: "Distance between adjacent passes".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.500),
            range: Some((0.001, 10.0)),
        });
        tree.register_param(ParamDef {
            key: ParamKey::new("side_stock"),
            display_name: "Side Stock to Leave".into(),
            description: "Material to leave on side walls".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.010),
            range: Some((0.0, 1.0)),
        });
        tree.register_param(ParamDef {
            key: ParamKey::new("bottom_stock"),
            display_name: "Bottom Stock to Leave".into(),
            description: "Material to leave on floor".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.010),
            range: Some((0.0, 1.0)),
        });
        tree
    }

    #[test]
    fn resolve_inherits_from_folder() {
        let mut tree = make_tree_with_dict();
        let folder = tree.add_folder("ROUGHING", None);
        tree.set_param(folder, "stepover".into(), ParamValue::F64(0.300));

        let op = tree.add_operation("OP1 pocket", "adaptive_pocket", Some(folder));

        let resolved = tree.resolve(op, &ParamKey::new("stepover")).unwrap();
        assert_eq!(resolved.value, ParamValue::F64(0.300));
        assert_eq!(resolved.source, folder);
        assert!(!resolved.is_default);
    }

    #[test]
    fn resolve_local_override_wins() {
        let mut tree = make_tree_with_dict();
        let folder = tree.add_folder("ROUGHING", None);
        tree.set_param(folder, "stepover".into(), ParamValue::F64(0.300));

        let op = tree.add_operation("OP1 pocket", "adaptive_pocket", Some(folder));
        tree.set_param(op, "stepover".into(), ParamValue::F64(0.200));

        let resolved = tree.resolve(op, &ParamKey::new("stepover")).unwrap();
        assert_eq!(resolved.value, ParamValue::F64(0.200));
        assert_eq!(resolved.source, op);
    }

    #[test]
    fn resolve_falls_through_to_default() {
        let mut tree = make_tree_with_dict();
        let op = tree.add_operation("OP1 pocket", "adaptive_pocket", None);

        let resolved = tree.resolve(op, &ParamKey::new("stepover")).unwrap();
        assert_eq!(resolved.value, ParamValue::F64(0.500));
        assert!(resolved.is_default);
    }

    #[test]
    fn resolve_nested_folders() {
        let mut tree = make_tree_with_dict();
        let root = tree.add_folder("ALL_OPS", None);
        tree.set_param(root, "side_stock".into(), ParamValue::F64(0.020));

        let sub = tree.add_folder("ROUGHING", Some(root));
        tree.set_param(sub, "stepover".into(), ParamValue::F64(0.400));
        // side_stock NOT set on sub — should inherit from root

        let op = tree.add_operation("OP1", "adaptive_pocket", Some(sub));

        // stepover from ROUGHING
        let r1 = tree.resolve(op, &ParamKey::new("stepover")).unwrap();
        assert_eq!(r1.value, ParamValue::F64(0.400));
        assert_eq!(r1.source, sub);

        // side_stock from ALL_OPS
        let r2 = tree.resolve(op, &ParamKey::new("side_stock")).unwrap();
        assert_eq!(r2.value, ParamValue::F64(0.020));
        assert_eq!(r2.source, root);
    }

    #[test]
    fn resolve_all_returns_every_dict_param() {
        let mut tree = make_tree_with_dict();
        let op = tree.add_operation("OP1", "adaptive_pocket", None);
        let all = tree.resolve_all(op);
        assert_eq!(all.len(), 3); // stepover, side_stock, bottom_stock
        assert!(all.values().all(|r| r.is_default));
    }
}

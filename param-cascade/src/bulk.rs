use tracing::debug_span;

use crate::tree::OpTree;
use crate::types::*;

impl OpTree {
    /// Set a param on multiple nodes. Locked nodes are skipped.
    pub fn bulk_set(
        &mut self,
        nodes: &[NodeId],
        key: ParamKey,
        value: ParamValue,
    ) -> BulkEditResult {
        let _span = debug_span!("bulk_set", param = %key, count = nodes.len()).entered();

        let mut updated = Vec::new();
        let mut skipped_locked = Vec::new();

        for &node_id in nodes {
            let is_locked = self
                .local_param(node_id, &key)
                .map(|e| e.state == ParamState::Locked)
                .unwrap_or(false);

            if is_locked {
                tracing::debug!(node = node_id.0, "skipping locked param");
                skipped_locked.push(node_id);
            } else {
                self.set_param(node_id, key.clone(), value.clone());
                updated.push(node_id);
            }
        }

        BulkEditResult {
            updated,
            skipped_locked,
        }
    }

    /// Set a param on a folder and propagate to all descendant operations.
    /// The folder itself gets the value set. Descendant ops that have the
    /// param locked are skipped. Descendant ops that have an unlocked
    /// local override get it cleared (revert to inherit from folder).
    /// Descendant ops inheriting already pick up the folder's new value.
    pub fn cascade_set(
        &mut self,
        folder: NodeId,
        key: ParamKey,
        value: ParamValue,
    ) -> BulkEditResult {
        let _span = debug_span!("cascade_set", param = %key, folder = folder.0).entered();

        // Set on the folder itself
        self.set_param(folder, key.clone(), value.clone());

        let ops = self.descendant_ops(folder);
        let mut updated = vec![folder];
        let mut skipped_locked = Vec::new();

        for op_id in ops {
            let local = self.local_param(op_id, &key);
            match local.map(|e| e.state) {
                Some(ParamState::Locked) => {
                    tracing::debug!(node = op_id.0, "skipping locked param in cascade");
                    skipped_locked.push(op_id);
                }
                Some(ParamState::Override) => {
                    // Clear the local override so it inherits from folder
                    self.clear_param(op_id, &key);
                    updated.push(op_id);
                }
                _ => {
                    // Already inheriting — no action needed, value cascades automatically
                    updated.push(op_id);
                }
            }
        }

        BulkEditResult {
            updated,
            skipped_locked,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tree() -> OpTree {
        let mut tree = OpTree::new();
        tree.register_param(ParamDef {
            key: ParamKey::new("stepover"),
            display_name: "Stepover".into(),
            description: "Distance between passes".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.500),
            range: None,
        });
        tree.register_param(ParamDef {
            key: ParamKey::new("side_stock"),
            display_name: "Side Stock".into(),
            description: "Stock to leave on sides".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.010),
            range: None,
        });
        tree
    }

    #[test]
    fn bulk_set_skips_locked() {
        let mut tree = make_tree();
        let folder = tree.add_folder("ROUGH", None);
        let op1 = tree.add_operation("OP1", "pocket", Some(folder));
        let op2 = tree.add_operation("OP2", "pocket", Some(folder));
        let op3 = tree.add_operation("OP3", "pocket", Some(folder));

        // Lock OP2's stepover
        tree.set_param_locked(op2, "stepover".into(), ParamValue::F64(0.250));

        let result = tree.bulk_set(
            &[op1, op2, op3],
            "stepover".into(),
            ParamValue::F64(0.100),
        );

        assert_eq!(result.updated.len(), 2);
        assert_eq!(result.skipped_locked.len(), 1);
        assert_eq!(result.skipped_locked[0], op2);

        // OP1 and OP3 got the new value
        assert_eq!(
            tree.local_param(op1, &ParamKey::new("stepover"))
                .unwrap()
                .value,
            ParamValue::F64(0.100)
        );
        // OP2 kept its locked value
        assert_eq!(
            tree.local_param(op2, &ParamKey::new("stepover"))
                .unwrap()
                .value,
            ParamValue::F64(0.250)
        );
    }

    #[test]
    fn cascade_set_clears_overrides() {
        let mut tree = make_tree();
        let folder = tree.add_folder("FINISH", None);
        let op1 = tree.add_operation("OP1", "contour", Some(folder));
        let op2 = tree.add_operation("OP2", "contour", Some(folder));

        // OP1 has a local override, OP2 inherits
        tree.set_param(op1, "side_stock".into(), ParamValue::F64(0.005));

        // Cascade set on folder
        let result = tree.cascade_set(folder, "side_stock".into(), ParamValue::F64(0.000));

        assert!(result.skipped_locked.is_empty());

        // OP1's local override should be cleared — now inherits from folder
        assert!(tree
            .local_param(op1, &ParamKey::new("side_stock"))
            .is_none());

        // Both ops should resolve to 0.000 from folder
        let r1 = tree.resolve(op1, &ParamKey::new("side_stock")).unwrap();
        assert_eq!(r1.value, ParamValue::F64(0.000));
        let r2 = tree.resolve(op2, &ParamKey::new("side_stock")).unwrap();
        assert_eq!(r2.value, ParamValue::F64(0.000));
    }

    #[test]
    fn cascade_set_respects_locks() {
        let mut tree = make_tree();
        let folder = tree.add_folder("FINISH", None);
        let op1 = tree.add_operation("OP1", "contour", Some(folder));
        let op2 = tree.add_operation("OP2", "contour", Some(folder));

        // Lock OP2's side_stock — customer callout
        tree.set_param_locked(op2, "side_stock".into(), ParamValue::F64(0.015));

        let result = tree.cascade_set(folder, "side_stock".into(), ParamValue::F64(0.000));

        assert_eq!(result.skipped_locked.len(), 1);
        assert_eq!(result.skipped_locked[0], op2);

        // OP1 inherits from folder → 0.000
        let r1 = tree.resolve(op1, &ParamKey::new("side_stock")).unwrap();
        assert_eq!(r1.value, ParamValue::F64(0.000));

        // OP2 still locked at 0.015
        let r2 = tree.resolve(op2, &ParamKey::new("side_stock")).unwrap();
        assert_eq!(r2.value, ParamValue::F64(0.015));
        assert!(r2.is_locked);
    }

    #[test]
    fn clone_rough_to_finish_workflow() {
        let mut tree = make_tree();

        // Build roughing program
        let rough = tree.add_folder("ROUGHING", None);
        tree.set_param(rough, "side_stock".into(), ParamValue::F64(0.010));
        tree.set_param(rough, "stepover".into(), ParamValue::F64(0.500));

        let _op1 = tree.add_operation("OP1 pocket", "adaptive_pocket", Some(rough));
        let _op2 = tree.add_operation("OP2 pocket", "adaptive_pocket", Some(rough));
        let op3 = tree.add_operation("OP3 contour", "contour", Some(rough));
        // Lock OP3's side_stock — customer spec
        tree.set_param_locked(op3, "side_stock".into(), ParamValue::F64(0.015));

        // Clone to finish
        let finish = tree.clone_subtree(rough, None);
        // Rename the folder
        if let Some(n) = tree.node_mut(finish) {
            n.name = "FINISH_CUTS".to_string();
        }

        // Change folder params for finishing
        tree.cascade_set(finish, "side_stock".into(), ParamValue::F64(0.000));
        tree.cascade_set(finish, "stepover".into(), ParamValue::F64(0.030));

        // Verify finish ops
        let finish_ops = tree.descendant_ops(finish);
        assert_eq!(finish_ops.len(), 3);

        for &op in &finish_ops {
            let r = tree.resolve(op, &ParamKey::new("stepover")).unwrap();
            assert_eq!(r.value, ParamValue::F64(0.030), "stepover should be 0.030");
        }

        // Check side_stock — the cloned OP3 should still be locked at 0.015
        let cloned_op3 = finish_ops[2]; // last op in clone order
        let r = tree.resolve(cloned_op3, &ParamKey::new("side_stock")).unwrap();
        assert_eq!(r.value, ParamValue::F64(0.015));
        assert!(r.is_locked);

        // Other finish ops should have side_stock = 0.000
        let r1 = tree
            .resolve(finish_ops[0], &ParamKey::new("side_stock"))
            .unwrap();
        assert_eq!(r1.value, ParamValue::F64(0.000));
    }
}

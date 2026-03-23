use crate::tree::OpTree;
use crate::types::*;

impl OpTree {
    /// Find all leaf operations where a resolved param matches a predicate.
    pub fn query<F>(&self, key: &ParamKey, pred: F) -> Vec<NodeId>
    where
        F: Fn(&ParamValue) -> bool,
    {
        self.all_ops()
            .into_iter()
            .filter(|&op| {
                self.resolve(op, key)
                    .map(|r| pred(&r.value))
                    .unwrap_or(false)
            })
            .collect()
    }

    /// Find all leaf operations where a resolved param equals a value.
    pub fn query_eq(&self, key: &ParamKey, value: &ParamValue) -> Vec<NodeId> {
        self.query(key, |v| v == value)
    }

    /// Find all leaf operations where a resolved numeric param is in [min, max].
    pub fn query_range(&self, key: &ParamKey, min: f64, max: f64) -> Vec<NodeId> {
        self.query(key, |v| {
            v.as_f64().map(|n| n >= min && n <= max).unwrap_or(false)
        })
    }

    /// Find all leaf operations where a resolved numeric param is > threshold.
    pub fn query_gt(&self, key: &ParamKey, threshold: f64) -> Vec<NodeId> {
        self.query(key, |v| {
            v.as_f64().map(|n| n > threshold).unwrap_or(false)
        })
    }

    /// Find all leaf operations where a resolved numeric param is < threshold.
    pub fn query_lt(&self, key: &ParamKey, threshold: f64) -> Vec<NodeId> {
        self.query(key, |v| {
            v.as_f64().map(|n| n < threshold).unwrap_or(false)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tree() -> OpTree {
        let mut tree = OpTree::new();
        tree.register_param(ParamDef {
            key: ParamKey::new("bottom_stock"),
            display_name: "Bottom Stock".into(),
            description: "Stock to leave on floor".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.010),
            range: None,
        });
        tree.register_param(ParamDef {
            key: ParamKey::new("cut_direction"),
            display_name: "Cut Direction".into(),
            description: "Climb or conventional".into(),
            unit: ParamUnit::None,
            value_type: ValueType::Enum(vec!["Climb".into(), "Conventional".into()]),
            default: ParamValue::Enum("Climb".into()),
            range: None,
        });
        tree
    }

    #[test]
    fn query_where_leaving_stock() {
        let mut tree = make_tree();
        let rough = tree.add_folder("ROUGHING", None);
        tree.set_param(rough, "bottom_stock".into(), ParamValue::F64(0.010));

        let op1 = tree.add_operation("OP1", "pocket", Some(rough));
        let op2 = tree.add_operation("OP2", "pocket", Some(rough));
        // OP2 overrides to zero
        tree.set_param(op2, "bottom_stock".into(), ParamValue::F64(0.0));

        let op3 = tree.add_operation("OP3", "pocket", None);
        tree.set_param(op3, "bottom_stock".into(), ParamValue::F64(0.005));

        // "Where am I leaving stock on the bottom?"
        let results = tree.query_gt(&ParamKey::new("bottom_stock"), 0.0);
        assert_eq!(results.len(), 2); // OP1 (inherits 0.010) and OP3 (0.005)
        assert!(results.contains(&op1));
        assert!(results.contains(&op3));
        assert!(!results.contains(&op2));
    }

    #[test]
    fn query_by_enum_value() {
        let mut tree = make_tree();
        let _op1 = tree.add_operation("OP1", "contour", None);
        let op2 = tree.add_operation("OP2", "contour", None);
        tree.set_param(
            op2,
            "cut_direction".into(),
            ParamValue::Enum("Conventional".into()),
        );

        // "Which ops are conventional?"
        let results =
            tree.query_eq(
                &ParamKey::new("cut_direction"),
                &ParamValue::Enum("Conventional".into()),
            );
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], op2);
    }

    #[test]
    fn query_range() {
        let mut tree = make_tree();
        let op1 = tree.add_operation("OP1", "pocket", None);
        tree.set_param(op1, "bottom_stock".into(), ParamValue::F64(0.005));
        let op2 = tree.add_operation("OP2", "pocket", None);
        tree.set_param(op2, "bottom_stock".into(), ParamValue::F64(0.015));
        let op3 = tree.add_operation("OP3", "pocket", None);
        // OP3 inherits default 0.010

        // Find ops leaving between 0.005 and 0.010 inclusive
        let results = tree.query_range(&ParamKey::new("bottom_stock"), 0.005, 0.010);
        assert_eq!(results.len(), 2); // OP1 (0.005) and OP3 (0.010 default)
        assert!(results.contains(&op1));
        assert!(results.contains(&op3));
    }
}

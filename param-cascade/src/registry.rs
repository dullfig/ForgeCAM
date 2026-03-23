use crate::tree::OpTree;
use crate::types::*;

impl OpTree {
    /// Validate an operation node against its registered op type.
    ///
    /// Checks:
    /// - All required params resolve to a value
    /// - All local param values match the dictionary's ValueType
    /// - All numeric values are within the dictionary's range
    /// - No unknown param keys are set locally
    pub fn validate_op(&self, node: NodeId) -> Vec<ValidationError> {
        let mut errors = Vec::new();

        let op_type = match self.node(node) {
            Some(n) => match &n.kind {
                NodeKind::Operation(t) => t.clone(),
                _ => return errors, // folders don't validate
            },
            None => return errors,
        };

        // Check required params resolve
        if let Some(reg) = self.op_types.get(&op_type) {
            for key in &reg.required_params {
                match self.resolve(node, key) {
                    None => {
                        errors.push(ValidationError::MissingRequired {
                            key: key.clone(),
                            op_type: op_type.clone(),
                        });
                    }
                    Some(resolved) => {
                        // Check type match
                        if let Some(def) = self.dictionary.get(key) {
                            if let Some(err) = check_type_match(key, &def.value_type, &resolved.value) {
                                errors.push(err);
                            }
                            // Check range
                            if let Some((min, max)) = def.range {
                                if let Some(v) = resolved.value.as_f64() {
                                    if v < min || v > max {
                                        errors.push(ValidationError::OutOfRange {
                                            key: key.clone(),
                                            value: v,
                                            min,
                                            max,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Check for unknown local params
        if let Some(n) = self.node(node) {
            for key in n.params.keys() {
                if !self.dictionary.contains_key(key) {
                    errors.push(ValidationError::UnknownParam { key: key.clone() });
                }
            }
        }

        errors
    }

    /// Validate all operations in the tree.
    pub fn validate_all(&self) -> Vec<(NodeId, Vec<ValidationError>)> {
        self.all_ops()
            .into_iter()
            .map(|op| {
                let errors = self.validate_op(op);
                (op, errors)
            })
            .filter(|(_, errors)| !errors.is_empty())
            .collect()
    }
}

fn check_type_match(key: &ParamKey, expected: &ValueType, got: &ParamValue) -> Option<ValidationError> {
    let mismatch = match (expected, got) {
        (ValueType::F64, ParamValue::F64(_)) => false,
        (ValueType::Bool, ParamValue::Bool(_)) => false,
        (ValueType::Int, ParamValue::Int(_)) => false,
        (ValueType::String, ParamValue::String(_)) => false,
        (ValueType::Enum(choices), ParamValue::Enum(v)) => !choices.contains(v),
        _ => true,
    };

    if mismatch {
        Some(ValidationError::TypeMismatch {
            key: key.clone(),
            expected: expected.clone(),
            got: got.clone(),
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_missing_required() {
        let mut tree = OpTree::new();
        tree.register_param(ParamDef {
            key: ParamKey::new("spindle_rpm"),
            display_name: "Spindle Speed".into(),
            description: "RPM".into(),
            unit: ParamUnit::SpindleRpm,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.0), // zero = unset sentinel
            range: Some((1.0, 100_000.0)),
        });

        tree.register_op_type(OpTypeRegistration {
            op_type: "drill".into(),
            display_name: "Drill".into(),
            required_params: vec![ParamKey::new("spindle_rpm")],
            optional_params: vec![],
        });

        let op = tree.add_operation("OP1 drill", "drill", None);
        // Don't set spindle_rpm — it'll resolve to default (0.0), which is out of range

        let errors = tree.validate_op(op);
        assert!(errors.iter().any(|e| matches!(e,
            ValidationError::OutOfRange { key, .. } if key.0 == "spindle_rpm"
        )));
    }

    #[test]
    fn validate_unknown_param() {
        let mut tree = OpTree::new();
        let op = tree.add_operation("OP1", "pocket", None);
        tree.set_param(op, "my_custom_thing".into(), ParamValue::F64(42.0));

        let errors = tree.validate_op(op);
        assert!(errors.iter().any(|e| matches!(e,
            ValidationError::UnknownParam { key } if key.0 == "my_custom_thing"
        )));
    }

    #[test]
    fn validate_enum_mismatch() {
        let mut tree = OpTree::new();
        tree.register_param(ParamDef {
            key: ParamKey::new("cut_direction"),
            display_name: "Cut Direction".into(),
            description: "Climb or conventional".into(),
            unit: ParamUnit::None,
            value_type: ValueType::Enum(vec!["Climb".into(), "Conventional".into()]),
            default: ParamValue::Enum("Climb".into()),
            range: None,
        });

        tree.register_op_type(OpTypeRegistration {
            op_type: "contour".into(),
            display_name: "Contour".into(),
            required_params: vec![ParamKey::new("cut_direction")],
            optional_params: vec![],
        });

        let op = tree.add_operation("OP1", "contour", None);
        tree.set_param(
            op,
            "cut_direction".into(),
            ParamValue::Enum("Zigzag".into()), // not a valid choice
        );

        let errors = tree.validate_op(op);
        assert!(errors.iter().any(|e| matches!(e,
            ValidationError::TypeMismatch { key, .. } if key.0 == "cut_direction"
        )));
    }

    #[test]
    fn validate_clean_op() {
        let mut tree = OpTree::new();
        tree.register_param(ParamDef {
            key: ParamKey::new("stepover"),
            display_name: "Stepover".into(),
            description: "Pass spacing".into(),
            unit: ParamUnit::Length,
            value_type: ValueType::F64,
            default: ParamValue::F64(0.500),
            range: Some((0.001, 10.0)),
        });

        tree.register_op_type(OpTypeRegistration {
            op_type: "pocket".into(),
            display_name: "Pocket".into(),
            required_params: vec![ParamKey::new("stepover")],
            optional_params: vec![],
        });

        let folder = tree.add_folder("ROUGH", None);
        tree.set_param(folder, "stepover".into(), ParamValue::F64(0.300));
        let op = tree.add_operation("OP1", "pocket", Some(folder));

        let errors = tree.validate_op(op);
        assert!(errors.is_empty());
    }
}

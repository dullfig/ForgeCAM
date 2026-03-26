//! STEP Part 21 file parser.
//!
//! Parses the ISO 10303-21 clear text encoding into an entity graph.
//! This is format-level parsing only — it doesn't interpret what entities
//! mean (that's the geometry builder and topology assembler's job).
//!
//! # Format overview
//!
//! ```text
//! ISO-10303-21;
//! HEADER;
//! FILE_DESCRIPTION((...), '2;1');
//! FILE_NAME('part.stp', '2024-01-15', ('author'), ('org'), ...);
//! FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MBD_MIM_LF'));
//! ENDSEC;
//! DATA;
//! #1 = CARTESIAN_POINT('', (0.0, 0.0, 0.0));
//! #2 = DIRECTION('', (0.0, 0.0, 1.0));
//! ...
//! ENDSEC;
//! END-ISO-10303-21;
//! ```

mod tokenizer;
mod entities;

pub use entities::{StepEntity, Parameter, StepFile, SchemaType, HeaderData};

use std::collections::HashMap;
use std::path::Path;

/// Parse a STEP file from a file path.
pub fn parse_step_file(path: &Path) -> Result<StepFile, ParseError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| ParseError::Io(e.to_string()))?;
    parse_step_string(&content)
}

/// Parse a STEP file from a string.
pub fn parse_step_string(input: &str) -> Result<StepFile, ParseError> {
    let tokens = tokenizer::tokenize(input)?;
    entities::build_step_file(&tokens)
}

/// Errors that can occur during STEP parsing.
#[derive(Debug, Clone)]
pub enum ParseError {
    Io(String),
    /// Unexpected token or structure.
    Syntax(String),
    /// Missing required section (HEADER, DATA).
    MissingSection(String),
    /// Entity reference (#id) not found.
    UnresolvedReference(u64),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::Io(e) => write!(f, "IO error: {e}"),
            ParseError::Syntax(msg) => write!(f, "syntax error: {msg}"),
            ParseError::MissingSection(s) => write!(f, "missing section: {s}"),
            ParseError::UnresolvedReference(id) => write!(f, "unresolved reference: #{id}"),
        }
    }
}

impl std::error::Error for ParseError {}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_STEP: &str = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''), '2;1');
FILE_NAME('test.stp', '2024-01-01', (''), (''), '', '', '');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('origin', (0.0, 0.0, 0.0));
#2 = DIRECTION('z_axis', (0.0, 0.0, 1.0));
#3 = DIRECTION('x_axis', (1.0, 0.0, 0.0));
#4 = AXIS2_PLACEMENT_3D('placement', #1, #2, #3);
#5 = PLANE('base_plane', #4);
ENDSEC;
END-ISO-10303-21;
"#;

    #[test]
    fn parse_minimal_step() {
        let result = parse_step_string(MINIMAL_STEP);
        let file = result.expect("should parse minimal STEP file");

        assert_eq!(file.schema, SchemaType::AP203);
        assert_eq!(file.entities.len(), 5);

        // Check entity #1 is a CARTESIAN_POINT
        let e1 = &file.entities[&1];
        assert_eq!(e1.type_name, "CARTESIAN_POINT");

        // Check entity #5 is a PLANE referencing #4
        let e5 = &file.entities[&5];
        assert_eq!(e5.type_name, "PLANE");
        match &e5.params[1] {
            Parameter::Ref(id) => assert_eq!(*id, 4),
            other => panic!("expected Ref(4), got {:?}", other),
        }
    }

    #[test]
    fn detect_ap242_schema() {
        let input = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''), '2;1');
FILE_NAME('part.stp', '2024-01-01', (''), (''), '', '', '');
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MBD_MIM_LF'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;
"#;
        let file = parse_step_string(input).expect("should parse");
        assert_eq!(file.schema, SchemaType::AP242);
    }

    #[test]
    fn detect_ap214_schema() {
        let input = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''), '2;1');
FILE_NAME('part.stp', '2024-01-01', (''), (''), '', '', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;
"#;
        let file = parse_step_string(input).expect("should parse");
        assert_eq!(file.schema, SchemaType::AP214);
    }

    const COMPLEX_PARAMS: &str = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''), '2;1');
FILE_NAME('test.stp', '2024-01-01', (''), (''), '', '', '');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('', (1.5, -2.3E1, 4.0E-2));
#2 = B_SPLINE_CURVE_WITH_KNOTS('', 3, (#10, #11, #12, #13), .UNSPECIFIED., .F., .F., (4, 4), (0.0, 1.0), .UNSPECIFIED.);
#3 = ORIENTED_EDGE('', *, *, #50, .T.);
ENDSEC;
END-ISO-10303-21;
"#;

    #[test]
    fn parse_complex_parameters() {
        let file = parse_step_string(COMPLEX_PARAMS).expect("should parse");

        // Scientific notation
        let e1 = &file.entities[&1];
        match &e1.params[1] {
            Parameter::List(items) => {
                assert_eq!(items.len(), 3);
                match &items[1] {
                    Parameter::Float(v) => assert!((*v - (-23.0)).abs() < 1e-10),
                    other => panic!("expected Float, got {:?}", other),
                }
            }
            other => panic!("expected List, got {:?}", other),
        }

        // ORIENTED_EDGE('', *, *, #50, .T.)
        // params: [String(""), Unset, Unset, Ref(50), Enum("T")]
        let e3 = &file.entities[&3];
        match &e3.params[4] {
            Parameter::Enum(s) => assert_eq!(s, "T"),
            other => panic!("expected Enum('T'), got {:?}", other),
        }

        // Unset values (*)
        match &e3.params[1] {
            Parameter::Unset => {}
            other => panic!("expected Unset, got {:?}", other),
        }
    }

    #[test]
    fn parse_real_step_file() {
        // Real STEP file from the server — DRS Leonardo part, 9KB.
        let path = std::path::Path::new(
            r"Q:\001 -- PRINTS\DRS Leonardo\7016080-001 REV 1\7016080-001_1_2.stp"
        );
        if !path.exists() {
            eprintln!("skipping real file test — Q: drive not available");
            return;
        }
        let file = parse_step_file(path).expect("should parse real STEP file");

        eprintln!("Schema: {:?} ({})", file.schema, file.header.schema_name);
        eprintln!("Entities: {}", file.entities.len());
        eprintln!("File: {}", file.header.file_name);
        eprintln!("Date: {}", file.header.date);

        // Count entity types.
        let mut type_counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
        for e in file.entities.values() {
            *type_counts.entry(&e.type_name).or_insert(0) += 1;
        }
        let mut sorted: Vec<_> = type_counts.iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(a.1));
        eprintln!("\nTop entity types:");
        for (name, count) in sorted.iter().take(20) {
            eprintln!("  {count:>5}  {name}");
        }

        assert!(file.entities.len() > 0, "should have entities");
    }
}

//! STEP entity types and file structure.

use std::collections::HashMap;

use super::tokenizer::Token;
use super::ParseError;

/// A parsed STEP file.
#[derive(Debug, Clone)]
pub struct StepFile {
    /// Header metadata.
    pub header: HeaderData,
    /// Detected application protocol.
    pub schema: SchemaType,
    /// All data entities, keyed by their #id.
    pub entities: HashMap<u64, StepEntity>,
}

/// Header section data.
#[derive(Debug, Clone, Default)]
pub struct HeaderData {
    pub description: Vec<String>,
    pub file_name: String,
    pub date: String,
    pub schema_name: String,
}

/// Which application protocol the file uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SchemaType {
    AP203,
    AP214,
    AP242,
    Unknown,
}

/// A single STEP entity: `#id = TYPE_NAME(params);`
#[derive(Debug, Clone)]
pub struct StepEntity {
    /// Entity id (the number after #).
    pub id: u64,
    /// Entity type name (e.g., "CARTESIAN_POINT", "ADVANCED_FACE").
    pub type_name: String,
    /// Parameters in order.
    pub params: Vec<Parameter>,
}

/// A parameter value in a STEP entity.
#[derive(Debug, Clone)]
pub enum Parameter {
    /// Integer value.
    Int(i64),
    /// Floating-point value.
    Float(f64),
    /// String value (without quotes).
    String(String),
    /// Reference to another entity (#id).
    Ref(u64),
    /// Enumeration value (.NAME.).
    Enum(String),
    /// List of parameters (parenthesized).
    List(Vec<Parameter>),
    /// Unset/derived value (*).
    Unset,
    /// Omitted value ($).
    Omitted,
    /// Typed parameter: TYPE_NAME(value, ...) — e.g., POSITIVE_LENGTH_MEASURE(1.0).
    Typed(String, Vec<Parameter>),
}

impl Parameter {
    /// Try to extract as a float (converts Int to Float).
    pub fn as_float(&self) -> Option<f64> {
        match self {
            Parameter::Float(v) => Some(*v),
            Parameter::Int(v) => Some(*v as f64),
            _ => None,
        }
    }

    /// Try to extract as a reference id.
    pub fn as_ref(&self) -> Option<u64> {
        match self {
            Parameter::Ref(id) => Some(*id),
            _ => None,
        }
    }

    /// Try to extract as a list.
    pub fn as_list(&self) -> Option<&[Parameter]> {
        match self {
            Parameter::List(items) => Some(items),
            _ => None,
        }
    }

    /// Try to extract as a string.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Parameter::String(s) => Some(s),
            _ => None,
        }
    }

    /// Try to extract as an enum name.
    pub fn as_enum(&self) -> Option<&str> {
        match self {
            Parameter::Enum(s) => Some(s),
            _ => None,
        }
    }

    /// Extract a list of floats (e.g., coordinate tuple).
    pub fn as_float_list(&self) -> Option<Vec<f64>> {
        let items = self.as_list()?;
        items.iter().map(|p| p.as_float()).collect()
    }

    /// Extract a list of references.
    pub fn as_ref_list(&self) -> Option<Vec<u64>> {
        let items = self.as_list()?;
        items.iter().map(|p| p.as_ref()).collect()
    }
}

/// Build a StepFile from a token stream.
pub fn build_step_file(tokens: &[Token]) -> Result<StepFile, ParseError> {
    let mut pos = 0;

    // Skip ISO-10303-21;
    skip_until_keyword(tokens, &mut pos, "HEADER")?;

    // Parse HEADER section.
    let (header, schema) = parse_header(tokens, &mut pos)?;

    // Skip to DATA section.
    skip_until_keyword(tokens, &mut pos, "DATA")?;

    // Parse DATA entities.
    let entities = parse_data_section(tokens, &mut pos)?;

    Ok(StepFile {
        header,
        schema,
        entities,
    })
}

fn skip_until_keyword(tokens: &[Token], pos: &mut usize, keyword: &str) -> Result<(), ParseError> {
    while *pos < tokens.len() {
        if let Token::Keyword(kw) = &tokens[*pos] {
            if kw.eq_ignore_ascii_case(keyword) {
                *pos += 1; // consume keyword
                // Skip the semicolon after it.
                if *pos < tokens.len() && matches!(&tokens[*pos], Token::Semi) {
                    *pos += 1;
                }
                return Ok(());
            }
        }
        *pos += 1;
    }
    Err(ParseError::MissingSection(keyword.to_string()))
}

fn parse_header(tokens: &[Token], pos: &mut usize) -> Result<(HeaderData, SchemaType), ParseError> {
    let mut header = HeaderData::default();
    let mut schema = SchemaType::Unknown;

    // Parse header entities until ENDSEC.
    while *pos < tokens.len() {
        match &tokens[*pos] {
            Token::Keyword(kw) if kw.eq_ignore_ascii_case("ENDSEC") => {
                *pos += 1;
                if *pos < tokens.len() && matches!(&tokens[*pos], Token::Semi) {
                    *pos += 1;
                }
                break;
            }
            Token::Keyword(kw) => {
                let name = kw.to_uppercase();
                *pos += 1;
                // Parse the parameter list.
                let params = parse_param_list(tokens, pos)?;
                // Skip trailing semicolon.
                if *pos < tokens.len() && matches!(&tokens[*pos], Token::Semi) {
                    *pos += 1;
                }

                match name.as_str() {
                    "FILE_SCHEMA" => {
                        // FILE_SCHEMA(('SCHEMA_NAME'))
                        if let Some(Parameter::List(items)) = params.first() {
                            if let Some(Parameter::String(s)) = items.first() {
                                header.schema_name = s.clone();
                                schema = detect_schema(s);
                            }
                        }
                    }
                    "FILE_NAME" => {
                        if let Some(Parameter::String(s)) = params.first() {
                            header.file_name = s.clone();
                        }
                        if let Some(Parameter::String(s)) = params.get(1) {
                            header.date = s.clone();
                        }
                    }
                    "FILE_DESCRIPTION" => {
                        if let Some(Parameter::List(items)) = params.first() {
                            for item in items {
                                if let Parameter::String(s) = item {
                                    header.description.push(s.clone());
                                }
                            }
                        }
                    }
                    _ => {} // ignore unknown header entities
                }
            }
            _ => { *pos += 1; }
        }
    }

    Ok((header, schema))
}

fn detect_schema(schema_str: &str) -> SchemaType {
    let upper = schema_str.to_uppercase();
    if upper.contains("AP242") || upper.contains("MBD") {
        SchemaType::AP242
    } else if upper.contains("AUTOMOTIVE_DESIGN") || upper.contains("AP214") {
        SchemaType::AP214
    } else if upper.contains("CONFIG_CONTROL") || upper.contains("AP203") {
        SchemaType::AP203
    } else {
        SchemaType::Unknown
    }
}

fn parse_data_section(tokens: &[Token], pos: &mut usize) -> Result<HashMap<u64, StepEntity>, ParseError> {
    let mut entities = HashMap::new();

    while *pos < tokens.len() {
        match &tokens[*pos] {
            Token::Keyword(kw) if kw.eq_ignore_ascii_case("ENDSEC") => {
                *pos += 1;
                break;
            }
            Token::EntityId(id) => {
                let id = *id;
                *pos += 1;

                // Expect '='
                if *pos < tokens.len() && matches!(&tokens[*pos], Token::Equals) {
                    *pos += 1;
                }

                if *pos >= tokens.len() {
                    break;
                }

                match &tokens[*pos] {
                    Token::Keyword(name) => {
                        // Simple entity: #id = TYPE_NAME(params);
                        let type_name = name.to_uppercase();
                        *pos += 1;
                        let params = parse_param_list(tokens, pos)?;
                        if *pos < tokens.len() && matches!(&tokens[*pos], Token::Semi) {
                            *pos += 1;
                        }
                        entities.insert(id, StepEntity { id, type_name, params });
                    }
                    Token::LParen => {
                        // Complex entity: #id = (TYPE_A(...) TYPE_B(...) ...);
                        // Multiple entity types sharing one ID.
                        // We store the first type as the primary and merge all params.
                        *pos += 1; // skip outer '('
                        let mut first_type = String::new();
                        let mut all_params = Vec::new();

                        while *pos < tokens.len() && !matches!(&tokens[*pos], Token::RParen) {
                            if let Token::Keyword(name) = &tokens[*pos] {
                                let sub_type = name.to_uppercase();
                                *pos += 1;
                                let sub_params = parse_param_list(tokens, pos)?;
                                if first_type.is_empty() {
                                    first_type = sub_type;
                                    all_params = sub_params;
                                }
                                // Additional types stored as typed params.
                                // (The geometry builder can inspect them if needed.)
                            } else {
                                *pos += 1; // skip unexpected tokens
                            }
                        }
                        if *pos < tokens.len() && matches!(&tokens[*pos], Token::RParen) {
                            *pos += 1;
                        }
                        if *pos < tokens.len() && matches!(&tokens[*pos], Token::Semi) {
                            *pos += 1;
                        }
                        if !first_type.is_empty() {
                            entities.insert(id, StepEntity { id, type_name: first_type, params: all_params });
                        }
                    }
                    other => {
                        return Err(ParseError::Syntax(
                            format!("expected entity type name after #{id} =, got {:?}", other),
                        ));
                    }
                }
            }
            _ => { *pos += 1; }
        }
    }

    Ok(entities)
}

/// Parse a parenthesized parameter list: `(param, param, ...)`
/// Returns the inner parameters (not including the parens).
fn parse_param_list(tokens: &[Token], pos: &mut usize) -> Result<Vec<Parameter>, ParseError> {
    if *pos >= tokens.len() || !matches!(&tokens[*pos], Token::LParen) {
        return Ok(Vec::new());
    }
    *pos += 1; // consume '('

    let mut params = Vec::new();

    loop {
        if *pos >= tokens.len() {
            return Err(ParseError::Syntax("unexpected end of file in parameter list".into()));
        }
        match &tokens[*pos] {
            Token::RParen => {
                *pos += 1;
                break;
            }
            Token::Comma => {
                *pos += 1;
                continue;
            }
            _ => {
                let param = parse_one_param(tokens, pos)?;
                params.push(param);
            }
        }
    }

    Ok(params)
}

/// Parse a single parameter value.
fn parse_one_param(tokens: &[Token], pos: &mut usize) -> Result<Parameter, ParseError> {
    if *pos >= tokens.len() {
        return Err(ParseError::Syntax("unexpected end of file".into()));
    }
    match &tokens[*pos] {
        Token::Integer(v) => {
            let v = *v;
            *pos += 1;
            Ok(Parameter::Int(v))
        }
        Token::Float(v) => {
            let v = *v;
            *pos += 1;
            Ok(Parameter::Float(v))
        }
        Token::StepString(s) => {
            let s = s.clone();
            *pos += 1;
            Ok(Parameter::String(s))
        }
        Token::EntityId(id) => {
            let id = *id;
            *pos += 1;
            Ok(Parameter::Ref(id))
        }
        Token::Enum(e) => {
            let e = e.clone();
            *pos += 1;
            Ok(Parameter::Enum(e))
        }
        Token::Star => {
            *pos += 1;
            Ok(Parameter::Unset)
        }
        Token::Dollar => {
            *pos += 1;
            Ok(Parameter::Omitted)
        }
        Token::LParen => {
            // Nested list.
            let inner = parse_param_list(tokens, pos)?;
            Ok(Parameter::List(inner))
        }
        Token::Keyword(name) => {
            // Typed parameter: TYPE_NAME(value) — e.g., POSITIVE_LENGTH_MEASURE(1.0)
            let name = name.clone();
            *pos += 1;
            if *pos < tokens.len() && matches!(&tokens[*pos], Token::LParen) {
                let inner = parse_param_list(tokens, pos)?;
                // Wrap as a typed parameter: store as a list with the type name.
                Ok(Parameter::Typed(name, inner))
            } else {
                // Bare keyword (shouldn't happen in well-formed STEP, but be lenient).
                Ok(Parameter::String(name))
            }
        }
        other => Err(ParseError::Syntax(format!("unexpected token in parameter: {:?}", other))),
    }
}

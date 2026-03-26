//! Tokenizer for STEP Part 21 files.
//!
//! Splits input text into tokens: keywords, entity IDs (#N), strings,
//! numbers, punctuation, enums (.NAME.), and special values (*, $).

use super::ParseError;

/// A token from a STEP Part 21 file.
#[derive(Debug, Clone)]
pub enum Token {
    /// A keyword or entity type name (e.g., "HEADER", "CARTESIAN_POINT").
    Keyword(String),
    /// Entity reference: #123
    EntityId(u64),
    /// Integer literal.
    Integer(i64),
    /// Floating-point literal (including scientific notation).
    Float(f64),
    /// String literal (content without quotes).
    StepString(String),
    /// Enumeration: .NAME. (content without dots).
    Enum(String),
    /// '('
    LParen,
    /// ')'
    RParen,
    /// ','
    Comma,
    /// ';'
    Semi,
    /// '='
    Equals,
    /// '*' (unset/derived value)
    Star,
    /// '$' (omitted value)
    Dollar,
}

/// Tokenize a STEP Part 21 string into a token stream.
pub fn tokenize(input: &str) -> Result<Vec<Token>, ParseError> {
    let mut tokens = Vec::new();
    let bytes = input.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    while i < len {
        let b = bytes[i];
        match b {
            // Whitespace — skip
            b' ' | b'\t' | b'\r' | b'\n' => { i += 1; }

            // Comment: /* ... */
            b'/' if i + 1 < len && bytes[i + 1] == b'*' => {
                i += 2;
                while i + 1 < len && !(bytes[i] == b'*' && bytes[i + 1] == b'/') {
                    i += 1;
                }
                i += 2; // skip */
            }

            // Punctuation
            b'(' => { tokens.push(Token::LParen); i += 1; }
            b')' => { tokens.push(Token::RParen); i += 1; }
            b',' => { tokens.push(Token::Comma); i += 1; }
            b';' => { tokens.push(Token::Semi); i += 1; }
            b'=' => { tokens.push(Token::Equals); i += 1; }
            b'*' => { tokens.push(Token::Star); i += 1; }
            b'$' => { tokens.push(Token::Dollar); i += 1; }

            // Entity ID: #123
            b'#' => {
                i += 1;
                let start = i;
                while i < len && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let num_str = &input[start..i];
                let id: u64 = num_str.parse().map_err(|_| {
                    ParseError::Syntax(format!("invalid entity ID: #{num_str}"))
                })?;
                tokens.push(Token::EntityId(id));
            }

            // String: 'content' (escaped '' becomes ')
            b'\'' => {
                i += 1;
                let start = i;
                loop {
                    if i >= len {
                        break;
                    }
                    if bytes[i] == b'\'' {
                        // Check for escaped quote ('')
                        if i + 1 < len && bytes[i + 1] == b'\'' {
                            i += 2; // skip both quotes, keep scanning
                        } else {
                            break; // end of string
                        }
                    } else {
                        i += 1;
                    }
                }
                let content = input[start..i].replace("''", "'");
                if i < len { i += 1; } // skip closing quote
                tokens.push(Token::StepString(content));
            }

            // Enum: .NAME.
            b'.' => {
                i += 1;
                let start = i;
                while i < len && bytes[i] != b'.' {
                    i += 1;
                }
                let name = input[start..i].to_string();
                i += 1; // skip closing dot
                tokens.push(Token::Enum(name));
            }

            // Number (integer or float, possibly negative, possibly scientific notation)
            b'-' | b'+' | b'0'..=b'9'
                if b.is_ascii_digit()
                    || ((b == b'-' || b == b'+') && i + 1 < len && bytes[i + 1].is_ascii_digit())
            => {
                let start = i;
                if b == b'-' || b == b'+' {
                    i += 1;
                }
                while i < len && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let mut is_float = false;
                // Decimal point
                if i < len && bytes[i] == b'.' {
                    // Check it's not an enum dot — peek ahead for a digit or E (scientific)
                    if i + 1 < len && (bytes[i + 1].is_ascii_digit() || bytes[i + 1] == b'E' || bytes[i + 1] == b'e') {
                        is_float = true;
                        i += 1; // skip '.'
                        while i < len && bytes[i].is_ascii_digit() {
                            i += 1;
                        }
                    }
                }
                // Scientific notation: E or e
                if i < len && (bytes[i] == b'E' || bytes[i] == b'e') {
                    is_float = true;
                    i += 1;
                    if i < len && (bytes[i] == b'+' || bytes[i] == b'-') {
                        i += 1;
                    }
                    while i < len && bytes[i].is_ascii_digit() {
                        i += 1;
                    }
                }

                let num_str = &input[start..i];
                if is_float {
                    let v: f64 = num_str.parse().map_err(|_| {
                        ParseError::Syntax(format!("invalid float: {num_str}"))
                    })?;
                    tokens.push(Token::Float(v));
                } else {
                    let v: i64 = num_str.parse().map_err(|_| {
                        ParseError::Syntax(format!("invalid integer: {num_str}"))
                    })?;
                    tokens.push(Token::Integer(v));
                }
            }

            // Keyword or entity type name: starts with letter or underscore
            _ if b.is_ascii_alphabetic() || b == b'_' => {
                let start = i;
                while i < len && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_' || bytes[i] == b'-') {
                    i += 1;
                }
                let word = input[start..i].to_string();
                tokens.push(Token::Keyword(word));
            }

            // Unknown character — skip (be lenient with real-world files)
            _ => { i += 1; }
        }
    }

    Ok(tokens)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenize_entity_line() {
        let input = "#1 = CARTESIAN_POINT('origin', (0.0, 0.0, 0.0));";
        let tokens = tokenize(input).unwrap();

        assert!(matches!(&tokens[0], Token::EntityId(1)));
        assert!(matches!(&tokens[1], Token::Equals));
        assert!(matches!(&tokens[2], Token::Keyword(k) if k == "CARTESIAN_POINT"));
        assert!(matches!(&tokens[3], Token::LParen));
        assert!(matches!(&tokens[4], Token::StepString(s) if s == "origin"));
    }

    #[test]
    fn tokenize_scientific_notation() {
        let input = "(1.5, -2.3E1, 4.0E-2)";
        let tokens = tokenize(input).unwrap();

        // 1.5
        match &tokens[1] {
            Token::Float(v) => assert!((*v - 1.5).abs() < 1e-10),
            other => panic!("expected Float(1.5), got {:?}", other),
        }
        // -2.3E1 = -23.0
        match &tokens[3] {
            Token::Float(v) => assert!((*v - (-23.0)).abs() < 1e-10),
            other => panic!("expected Float(-23.0), got {:?}", other),
        }
        // 4.0E-2 = 0.04
        match &tokens[5] {
            Token::Float(v) => assert!((*v - 0.04).abs() < 1e-10),
            other => panic!("expected Float(0.04), got {:?}", other),
        }
    }

    #[test]
    fn tokenize_enum_and_star() {
        let input = ".T., .UNSPECIFIED., *, $";
        let tokens = tokenize(input).unwrap();

        assert!(matches!(&tokens[0], Token::Enum(s) if s == "T"));
        assert!(matches!(&tokens[2], Token::Enum(s) if s == "UNSPECIFIED"));
        assert!(matches!(&tokens[4], Token::Star));
        assert!(matches!(&tokens[6], Token::Dollar));
    }

    #[test]
    fn tokenize_comment() {
        let input = "#1 = /* this is a comment */ POINT('', (0.0));";
        let tokens = tokenize(input).unwrap();

        assert!(matches!(&tokens[0], Token::EntityId(1)));
        assert!(matches!(&tokens[1], Token::Equals));
        assert!(matches!(&tokens[2], Token::Keyword(k) if k == "POINT"));
    }

    #[test]
    fn tokenize_escaped_string() {
        let input = "'it''s a test'";
        let tokens = tokenize(input).unwrap();
        match &tokens[0] {
            Token::StepString(s) => assert_eq!(s, "it's a test"),
            other => panic!("expected string, got {:?}", other),
        }
    }
}

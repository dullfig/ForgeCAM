use rustkernel_topology::topo::FaceIdx;

use crate::boolean::face_classifier::FacePosition;

/// Boolean operation type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    Fuse,   // union: A ∪ B
    Cut,    // difference: A − B
    Common, // intersection: A ∩ B
}

/// The set of faces to keep in the boolean result.
pub struct SelectedFaces {
    /// Faces from solid A to keep (unchanged normals).
    pub keep_from_a: Vec<FaceIdx>,
    /// Faces from solid B to keep (unchanged normals).
    pub keep_from_b: Vec<FaceIdx>,
    /// Faces from solid B to keep with flipped normals (used in Cut).
    pub flip_from_b: Vec<FaceIdx>,
}

/// Given classified faces from both solids, select which faces to keep
/// based on the boolean operation type.
pub fn select_faces(
    op: BooleanOp,
    faces_a: &[(FaceIdx, FacePosition)],
    faces_b: &[(FaceIdx, FacePosition)],
) -> SelectedFaces {
    let mut result = SelectedFaces {
        keep_from_a: Vec::new(),
        keep_from_b: Vec::new(),
        flip_from_b: Vec::new(),
    };

    for &(face, pos) in faces_a {
        let keep = match op {
            BooleanOp::Fuse => matches!(pos, FacePosition::Outside | FacePosition::OnSame),
            BooleanOp::Cut => matches!(pos, FacePosition::Outside | FacePosition::OnOpposite),
            BooleanOp::Common => matches!(pos, FacePosition::Inside | FacePosition::OnSame),
        };
        if keep {
            result.keep_from_a.push(face);
        }
    }

    for &(face, pos) in faces_b {
        match op {
            BooleanOp::Fuse => {
                if matches!(pos, FacePosition::Outside) {
                    result.keep_from_b.push(face);
                }
            }
            BooleanOp::Cut => {
                if matches!(pos, FacePosition::Inside) {
                    result.flip_from_b.push(face);
                }
            }
            BooleanOp::Common => {
                if matches!(pos, FacePosition::Inside) {
                    result.keep_from_b.push(face);
                }
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_topology::arena::Idx;

    fn face(n: u32) -> FaceIdx {
        Idx::from_raw(n)
    }

    #[test]
    fn test_fuse_selection() {
        let faces_a = vec![
            (face(0), FacePosition::Outside),
            (face(1), FacePosition::Inside),
            (face(2), FacePosition::OnSame),
        ];
        let faces_b = vec![
            (face(3), FacePosition::Outside),
            (face(4), FacePosition::Inside),
        ];

        let sel = select_faces(BooleanOp::Fuse, &faces_a, &faces_b);
        assert_eq!(sel.keep_from_a.len(), 2); // Outside + OnSame
        assert_eq!(sel.keep_from_b.len(), 1); // Outside
        assert_eq!(sel.flip_from_b.len(), 0);
    }

    #[test]
    fn test_cut_selection() {
        let faces_a = vec![
            (face(0), FacePosition::Outside),
            (face(1), FacePosition::Inside),
            (face(2), FacePosition::OnOpposite),
        ];
        let faces_b = vec![
            (face(3), FacePosition::Outside),
            (face(4), FacePosition::Inside),
        ];

        let sel = select_faces(BooleanOp::Cut, &faces_a, &faces_b);
        assert_eq!(sel.keep_from_a.len(), 2); // Outside + OnOpposite
        assert_eq!(sel.keep_from_b.len(), 0);
        assert_eq!(sel.flip_from_b.len(), 1); // Inside (flipped)
    }

    #[test]
    fn test_common_selection() {
        let faces_a = vec![
            (face(0), FacePosition::Outside),
            (face(1), FacePosition::Inside),
            (face(2), FacePosition::OnSame),
        ];
        let faces_b = vec![
            (face(3), FacePosition::Outside),
            (face(4), FacePosition::Inside),
        ];

        let sel = select_faces(BooleanOp::Common, &faces_a, &faces_b);
        assert_eq!(sel.keep_from_a.len(), 2); // Inside + OnSame
        assert_eq!(sel.keep_from_b.len(), 1); // Inside
        assert_eq!(sel.flip_from_b.len(), 0);
    }
}

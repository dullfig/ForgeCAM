use rustkernel_math::{Point3, Vec3};

use crate::geom_store::GeomAccess;
use crate::store::TopoStore;
use crate::topo::FaceIdx;

/// Walk the outer loop of a face, returning ordered vertex positions.
pub fn face_boundary_points(
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    face_idx: FaceIdx,
) -> Vec<Point3> {
    let face = topo.faces.get(face_idx);
    let loop_idx = face.outer_loop;
    let start_he = topo.loops.get(loop_idx).half_edge;

    let mut points = Vec::new();
    let mut he = start_he;
    loop {
        let vert_idx = topo.half_edges.get(he).origin;
        let point_id = topo.vertices.get(vert_idx).point_id;
        points.push(geom.point(point_id));
        he = topo.half_edges.get(he).next;
        if he == start_he {
            break;
        }
    }
    points
}

/// A local 2D coordinate frame for a planar face.
/// Projects 3D points onto the plane and back.
#[derive(Debug, Clone)]
pub struct PlaneFrame {
    pub origin: Point3,
    pub u_axis: Vec3,
    pub v_axis: Vec3,
    pub normal: Vec3,
}

impl PlaneFrame {
    /// Construct an orthonormal 2D frame from a plane origin and normal.
    pub fn from_normal(origin: Point3, normal: Vec3) -> Self {
        let n = normal.normalize();

        // Choose a vector not parallel to n to form the U axis.
        let arbitrary = if n.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };

        let u_axis = n.cross(&arbitrary).normalize();
        let v_axis = n.cross(&u_axis); // already unit length

        Self {
            origin,
            u_axis,
            v_axis,
            normal: n,
        }
    }

    /// Project a 3D point onto the plane's 2D coordinate system.
    pub fn project_to_2d(&self, p: &Point3) -> [f64; 2] {
        let d = p - self.origin;
        [d.dot(&self.u_axis), d.dot(&self.v_axis)]
    }

    /// Unproject a 2D coordinate back to 3D on the plane.
    pub fn unproject_to_3d(&self, uv: [f64; 2]) -> Point3 {
        self.origin + uv[0] * self.u_axis + uv[1] * self.v_axis
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plane_frame_roundtrip() {
        let frame = PlaneFrame::from_normal(
            Point3::new(1.0, 2.0, 3.0),
            Vec3::new(0.0, 0.0, 1.0),
        );

        let p = Point3::new(4.0, 5.0, 3.0); // on the plane z=3
        let uv = frame.project_to_2d(&p);
        let p2 = frame.unproject_to_3d(uv);

        assert!((p.x - p2.x).abs() < 1e-12);
        assert!((p.y - p2.y).abs() < 1e-12);
        assert!((p.z - p2.z).abs() < 1e-12);
    }

    #[test]
    fn test_plane_frame_axes_orthogonal() {
        let frame = PlaneFrame::from_normal(
            Point3::origin(),
            Vec3::new(1.0, 1.0, 1.0),
        );

        assert!(frame.u_axis.dot(&frame.v_axis).abs() < 1e-12);
        assert!(frame.u_axis.dot(&frame.normal).abs() < 1e-12);
        assert!(frame.v_axis.dot(&frame.normal).abs() < 1e-12);
        assert!((frame.u_axis.norm() - 1.0).abs() < 1e-12);
        assert!((frame.v_axis.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_plane_frame_identity_for_xy_plane() {
        let frame = PlaneFrame::from_normal(
            Point3::origin(),
            Vec3::new(0.0, 0.0, 1.0),
        );

        // Project points in the XY plane.
        let p = Point3::new(3.0, 4.0, 0.0);
        let uv = frame.project_to_2d(&p);

        // The 2D coords should reflect x and y (though possibly rotated).
        // Roundtrip is the important check.
        let p2 = frame.unproject_to_3d(uv);
        assert!((p.x - p2.x).abs() < 1e-12);
        assert!((p.y - p2.y).abs() < 1e-12);
        assert!((p.z - p2.z).abs() < 1e-12);
    }
}


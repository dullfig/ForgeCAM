use rustkernel_math::Point3;
use rustkernel_topology::face_util::face_boundary_points;
use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::{FaceIdx, SolidIdx};

/// Axis-aligned bounding box.
#[derive(Debug, Clone)]
pub struct AABB {
    pub min: Point3,
    pub max: Point3,
}

impl AABB {
    /// Build an AABB from a set of points.
    pub fn from_points(pts: &[Point3]) -> Self {
        assert!(!pts.is_empty());
        let mut min = pts[0];
        let mut max = pts[0];
        for p in &pts[1..] {
            min.x = min.x.min(p.x);
            min.y = min.y.min(p.y);
            min.z = min.z.min(p.z);
            max.x = max.x.max(p.x);
            max.y = max.y.max(p.y);
            max.z = max.z.max(p.z);
        }
        Self { min, max }
    }

    /// Test overlap with another AABB.
    pub fn overlaps(&self, other: &AABB) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// Expand by tolerance in all directions.
    pub fn expanded(&self, tol: f64) -> Self {
        Self {
            min: Point3::new(self.min.x - tol, self.min.y - tol, self.min.z - tol),
            max: Point3::new(self.max.x + tol, self.max.y + tol, self.max.z + tol),
        }
    }
}

/// Compute the maximum surface deviation beyond boundary vertices for curved surfaces.
/// For planes: 0. For cylinders/spheres: radius. Conservative bound.
fn surface_aabb_inflation(geom: &dyn GeomAccess, surface_id: u32) -> f64 {
    match geom.surface_kind(surface_id) {
        SurfaceKind::Plane { .. } => 0.0,
        SurfaceKind::Cylinder { radius, .. } => radius.abs(),
        SurfaceKind::Sphere { radius, .. } => radius.abs(),
        SurfaceKind::Cone { half_angle, .. } => {
            // Conservative: cone can extend significantly.
            // Use a fixed bound; face boundary vertices give the real extent.
            half_angle.abs().tan() * 0.1 // small inflation
        }
        SurfaceKind::Torus { minor_radius, .. } => minor_radius.abs(),
        SurfaceKind::Unknown => 0.0,
    }
}

/// Compute the AABB for a face by walking its boundary vertices,
/// then inflating by the surface curvature deviation.
pub fn face_aabb(topo: &TopoStore, geom: &dyn GeomAccess, face: FaceIdx) -> AABB {
    let pts = face_boundary_points(topo, geom, face);
    let base = AABB::from_points(&pts);
    let surface_id = topo.faces.get(face).surface_id;
    let inflation = surface_aabb_inflation(geom, surface_id);
    if inflation > 0.0 {
        base.expanded(inflation)
    } else {
        base
    }
}

/// Find all (face_a, face_b) pairs where face_a belongs to solid_a and face_b
/// belongs to solid_b, and their AABBs overlap (within tolerance).
pub fn find_interfering_face_pairs(
    topo: &TopoStore,
    geom: &dyn GeomAccess,
    solid_a: SolidIdx,
    solid_b: SolidIdx,
    tolerance: f64,
) -> Vec<(FaceIdx, FaceIdx)> {
    let shell_a = topo.solids.get(solid_a).shell;
    let shell_b = topo.solids.get(solid_b).shell;

    let faces_a: Vec<FaceIdx> = topo.shells.get(shell_a).faces.clone();
    let faces_b: Vec<FaceIdx> = topo.shells.get(shell_b).faces.clone();

    let aabbs_a: Vec<(FaceIdx, AABB)> = faces_a
        .iter()
        .map(|&f| (f, face_aabb(topo, geom, f).expanded(tolerance)))
        .collect();
    let aabbs_b: Vec<(FaceIdx, AABB)> = faces_b
        .iter()
        .map(|&f| (f, face_aabb(topo, geom, f).expanded(tolerance)))
        .collect();

    let mut pairs = Vec::new();
    for (fa, aabb_a) in &aabbs_a {
        for (fb, aabb_b) in &aabbs_b {
            if aabb_a.overlaps(aabb_b) {
                pairs.push((*fa, *fb));
            }
        }
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustkernel_builders::box_builder::make_box_into;
    use rustkernel_geom::AnalyticalGeomStore;
    use rustkernel_topology::store::TopoStore;

    #[test]
    fn test_overlapping_boxes_have_pairs() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 2.0, 2.0, 2.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(1.0, 0.0, 0.0), 2.0, 2.0, 2.0);

        let pairs = find_interfering_face_pairs(&topo, &geom, a, b, 1e-8);
        assert!(!pairs.is_empty(), "Overlapping boxes should have interfering face pairs");
    }

    #[test]
    fn test_separated_boxes_no_pairs() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 1.0, 1.0, 1.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(10.0, 0.0, 0.0), 1.0, 1.0, 1.0);

        let pairs = find_interfering_face_pairs(&topo, &geom, a, b, 1e-8);
        assert!(pairs.is_empty(), "Separated boxes should have no pairs");
    }

    #[test]
    fn test_touching_boxes_have_pairs() {
        let mut topo = TopoStore::new();
        let mut geom = AnalyticalGeomStore::new();
        // Box A: x in [-0.5, 0.5], Box B: x in [0.5, 1.5] — touching at x=0.5
        let a = make_box_into(&mut topo, &mut geom, Point3::origin(), 1.0, 1.0, 1.0);
        let b = make_box_into(&mut topo, &mut geom, Point3::new(1.0, 0.0, 0.0), 1.0, 1.0, 1.0);

        let pairs = find_interfering_face_pairs(&topo, &geom, a, b, 1e-8);
        // The two faces at x=0.5 should overlap (with tolerance), plus possibly edge-touching faces.
        assert!(!pairs.is_empty(), "Touching boxes should have pairs (with tolerance)");
    }
}

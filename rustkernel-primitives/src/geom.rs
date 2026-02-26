use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::geom_store::{GeomAccess, SurfaceKind};

/// A line segment between two points.
pub struct LineSegment {
    pub start: Point3,
    pub end: Point3,
}

/// A plane defined by an origin point and a normal vector.
pub struct Plane {
    pub origin: Point3,
    pub normal: Vec3,
}

/// Concrete geometry store holding analytical geometry (planes, line segments).
pub struct AnalyticalGeomStore {
    pub points: Vec<Point3>,
    pub curves: Vec<LineSegment>,
    pub surfaces: Vec<Plane>,
}

impl AnalyticalGeomStore {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            curves: Vec::new(),
            surfaces: Vec::new(),
        }
    }

    pub fn add_point(&mut self, p: Point3) -> u32 {
        let id = self.points.len() as u32;
        self.points.push(p);
        id
    }

    pub fn add_curve(&mut self, seg: LineSegment) -> u32 {
        let id = self.curves.len() as u32;
        self.curves.push(seg);
        id
    }

    pub fn add_surface(&mut self, plane: Plane) -> u32 {
        let id = self.surfaces.len() as u32;
        self.surfaces.push(plane);
        id
    }
}

impl Default for AnalyticalGeomStore {
    fn default() -> Self {
        Self::new()
    }
}

impl GeomAccess for AnalyticalGeomStore {
    fn point(&self, point_id: u32) -> Point3 {
        self.points[point_id as usize]
    }

    fn curve_eval(&self, curve_id: u32, t: f64) -> Point3 {
        let seg = &self.curves[curve_id as usize];
        Point3::from(seg.start.coords.lerp(&seg.end.coords, t))
    }

    fn surface_eval(&self, surface_id: u32, u: f64, v: f64) -> Point3 {
        let plane = &self.surfaces[surface_id as usize];
        // For a plane, we need a local frame. For M1 box faces this is unused
        // in tessellation (we read vertex points directly), but we provide a basic impl.
        let _ = (u, v);
        plane.origin
    }

    fn surface_normal(&self, surface_id: u32, _u: f64, _v: f64) -> Vec3 {
        self.surfaces[surface_id as usize].normal
    }

    fn surface_kind(&self, surface_id: u32) -> SurfaceKind {
        let plane = &self.surfaces[surface_id as usize];
        SurfaceKind::Plane {
            origin: plane.origin,
            normal: plane.normal,
        }
    }
}

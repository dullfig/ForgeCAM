use rustkernel_math::Point3;
use rustkernel_topology::topo::SolidIdx;
use rustkernel_builders::box_builder::make_box_into;
use rustkernel_builders::cylinder_builder::make_cylinder_into;
use rustkernel_builders::sphere_builder::make_sphere_into;
use rustkernel_builders::cone_builder::make_cone_into;
use rustkernel_builders::torus_builder::make_torus_into;

use crate::Kernel;

impl Kernel {
    /// Build an axis-aligned box centered at the origin.
    pub fn make_box(&mut self, dx: f64, dy: f64, dz: f64) -> SolidIdx {
        let s = make_box_into(&mut self.topo, &mut self.geom, Point3::origin(), dx, dy, dz);
        self.record_primitive_evolution(s);
        s
    }

    /// Build an axis-aligned box centered at `center`.
    pub fn make_box_at(&mut self, center: [f64; 3], dx: f64, dy: f64, dz: f64) -> SolidIdx {
        let c = Point3::new(center[0], center[1], center[2]);
        let s = make_box_into(&mut self.topo, &mut self.geom, c, dx, dy, dz);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a cylinder centered at the origin with +Z axis.
    pub fn make_cylinder(&mut self, radius: f64, height: f64) -> SolidIdx {
        let s = make_cylinder_into(&mut self.topo, &mut self.geom, Point3::origin(), radius, height, 32);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a cylinder at `center` with +Z axis.
    pub fn make_cylinder_at(&mut self, center: [f64; 3], radius: f64, height: f64) -> SolidIdx {
        let c = Point3::new(center[0], center[1], center[2]);
        let s = make_cylinder_into(&mut self.topo, &mut self.geom, c, radius, height, 32);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a sphere centered at the origin.
    pub fn make_sphere(&mut self, radius: f64) -> SolidIdx {
        let s = make_sphere_into(&mut self.topo, &mut self.geom, Point3::origin(), radius, 16, 8);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a sphere at `center`.
    pub fn make_sphere_at(&mut self, center: [f64; 3], radius: f64) -> SolidIdx {
        let c = Point3::new(center[0], center[1], center[2]);
        let s = make_sphere_into(&mut self.topo, &mut self.geom, c, radius, 16, 8);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a cone (or frustum) at the origin with +Z axis.
    /// `r1` is the bottom radius, `r2` is the top radius (0 for pointed cone).
    pub fn make_cone(&mut self, r1: f64, r2: f64, height: f64) -> SolidIdx {
        let s = make_cone_into(&mut self.topo, &mut self.geom, Point3::origin(), r1, r2, height, 32);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a cone at `center` with +Z axis.
    pub fn make_cone_at(&mut self, center: [f64; 3], r1: f64, r2: f64, height: f64) -> SolidIdx {
        let c = Point3::new(center[0], center[1], center[2]);
        let s = make_cone_into(&mut self.topo, &mut self.geom, c, r1, r2, height, 32);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a torus centered at the origin with +Z axis.
    pub fn make_torus(&mut self, major_r: f64, minor_r: f64) -> SolidIdx {
        let s = make_torus_into(&mut self.topo, &mut self.geom, Point3::origin(), major_r, minor_r, 16, 8);
        self.record_primitive_evolution(s);
        s
    }

    /// Build a torus at `center` with +Z axis.
    pub fn make_torus_at(&mut self, center: [f64; 3], major_r: f64, minor_r: f64) -> SolidIdx {
        let c = Point3::new(center[0], center[1], center[2]);
        let s = make_torus_into(&mut self.topo, &mut self.geom, c, major_r, minor_r, 16, 8);
        self.record_primitive_evolution(s);
        s
    }
}

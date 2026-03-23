use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::topo::SolidIdx;
use rustkernel_geom::CurveDef;
use rustkernel_builders::extrude_builder::make_extrude_into;
use rustkernel_builders::revolve_builder::make_revolve_into;
use rustkernel_builders::nurbs_extrude_builder::make_nurbs_extrude_solid_into;
use rustkernel_builders::nurbs_revolve_builder::make_nurbs_revolve_solid_into;
use curvo::prelude::{Interpolation, NurbsCurve3D, NurbsSurface3D};

use crate::Kernel;

impl Kernel {
    // --- Sketch + Extrude ---

    /// Create a new 2D sketch on a workplane defined by an origin and normal.
    pub fn create_sketch(&self, origin: Point3, normal: Vec3) -> rustkernel_sketch::Sketch {
        rustkernel_sketch::Sketch::new(origin, normal)
    }

    /// Extrude a closed 3D profile polygon along a direction.
    pub fn extrude(&mut self, profile: &[Point3], direction: Vec3, height: f64) -> SolidIdx {
        let s = make_extrude_into(&mut self.topo, &mut self.geom, profile, direction, height);
        self.record_primitive_evolution(s);
        s
    }

    /// Revolve a closed 3D profile polygon around an axis.
    pub fn revolve(&mut self, profile: &[Point3], axis_origin: Point3, axis_dir: Vec3, angle: f64) -> SolidIdx {
        let s = make_revolve_into(&mut self.topo, &mut self.geom, profile, axis_origin, axis_dir, angle, 32);
        self.record_primitive_evolution(s);
        s
    }

    // --- NURBS ---

    /// Create a NURBS curve interpolating through 3D points.
    pub fn interpolate_curve(&mut self, points: &[Point3], degree: usize) -> u32 {
        let pts: Vec<Point3> = points.to_vec();
        let curve = NurbsCurve3D::<f64>::interpolate(&pts, degree)
            .expect("NURBS interpolation failed");
        self.geom.add_nurbs_curve(curve)
    }

    /// Loft a NURBS surface through multiple NURBS curves.
    pub fn loft(&mut self, curve_ids: &[u32], degree_v: usize) -> u32 {
        let curves: Vec<NurbsCurve3D<f64>> = curve_ids
            .iter()
            .map(|&id| {
                match &self.geom.curves[id as usize] {
                    CurveDef::Nurbs(nc) => nc.clone(),
                    _ => panic!("loft requires NURBS curves, curve {} is not NURBS", id),
                }
            })
            .collect();
        let surface = NurbsSurface3D::<f64>::try_loft(&curves, Some(degree_v))
            .expect("NURBS loft failed");
        self.geom.add_nurbs_surface(surface)
    }

    /// Extrude a NURBS curve into a NURBS surface along direction * length.
    pub fn nurbs_extrude(&mut self, curve_id: u32, direction: Vec3, length: f64) -> u32 {
        let curve = match &self.geom.curves[curve_id as usize] {
            CurveDef::Nurbs(nc) => nc.clone(),
            _ => panic!("nurbs_extrude requires a NURBS curve, curve {} is not NURBS", curve_id),
        };
        let translation = direction.normalize() * length;
        let surface = NurbsSurface3D::<f64>::extrude(&curve, &translation);
        self.geom.add_nurbs_surface(surface)
    }

    /// Store a pre-built NURBS curve, returning its curve ID.
    pub fn add_nurbs_curve(&mut self, curve: NurbsCurve3D<f64>) -> u32 {
        self.geom.add_nurbs_curve(curve)
    }

    /// Store a pre-built NURBS surface, returning its surface ID.
    pub fn add_nurbs_surface(&mut self, surface: NurbsSurface3D<f64>) -> u32 {
        self.geom.add_nurbs_surface(surface)
    }

    // --- NURBS Solid Builders ---

    /// Extrude a NURBS curve into a solid along direction * height.
    /// The curve must be closed (start ≈ end). Returns a SolidIdx.
    pub fn make_nurbs_extrude_solid(
        &mut self,
        curve_id: u32,
        direction: Vec3,
        height: f64,
    ) -> SolidIdx {
        let curve = match &self.geom.curves[curve_id as usize] {
            CurveDef::Nurbs(nc) => nc.clone(),
            _ => panic!(
                "make_nurbs_extrude_solid requires a NURBS curve, curve {} is not NURBS",
                curve_id
            ),
        };
        let s = make_nurbs_extrude_solid_into(
            &mut self.topo,
            &mut self.geom,
            &curve,
            direction,
            height,
            32,
        );
        self.record_primitive_evolution(s);
        s
    }

    /// Revolve a NURBS curve into a solid around an axis.
    /// The curve must be closed (start ≈ end). Returns a SolidIdx.
    pub fn make_nurbs_revolve_solid(
        &mut self,
        curve_id: u32,
        axis_origin: Point3,
        axis_dir: Vec3,
        angle: f64,
    ) -> SolidIdx {
        let curve = match &self.geom.curves[curve_id as usize] {
            CurveDef::Nurbs(nc) => nc.clone(),
            _ => panic!(
                "make_nurbs_revolve_solid requires a NURBS curve, curve {} is not NURBS",
                curve_id
            ),
        };
        let s = make_nurbs_revolve_solid_into(
            &mut self.topo,
            &mut self.geom,
            &curve,
            axis_origin,
            axis_dir,
            angle,
            32,
            32,
        );
        self.record_primitive_evolution(s);
        s
    }
}

use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::topo::{SolidIdx, LoopIdx};
use tracing::info_span;

use crate::Kernel;

impl Kernel {
    /// Union of two solids.
    pub fn fuse(
        &mut self,
        a: SolidIdx,
        b: SolidIdx,
    ) -> Result<SolidIdx, rustkernel_boolean::ops::BooleanError> {
        let _span = info_span!("kernel.fuse", a = a.raw(), b = b.raw()).entered();
        let br = rustkernel_boolean::ops::boolean_op(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            a,
            b,
            rustkernel_boolean::face_selector::BooleanOp::Fuse,
        )?;
        self.last_evolution = Some(br.evolution);
        Ok(br.solid)
    }

    /// Subtraction: A minus B.
    pub fn cut(
        &mut self,
        a: SolidIdx,
        b: SolidIdx,
    ) -> Result<SolidIdx, rustkernel_boolean::ops::BooleanError> {
        let _span = info_span!("kernel.cut", a = a.raw(), b = b.raw()).entered();
        let br = rustkernel_boolean::ops::boolean_op(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            a,
            b,
            rustkernel_boolean::face_selector::BooleanOp::Cut,
        )?;
        self.last_evolution = Some(br.evolution);
        Ok(br.solid)
    }

    /// Intersection of two solids.
    pub fn common(
        &mut self,
        a: SolidIdx,
        b: SolidIdx,
    ) -> Result<SolidIdx, rustkernel_boolean::ops::BooleanError> {
        let _span = info_span!("kernel.common", a = a.raw(), b = b.raw()).entered();
        let br = rustkernel_boolean::ops::boolean_op(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            a,
            b,
            rustkernel_boolean::face_selector::BooleanOp::Common,
        )?;
        self.last_evolution = Some(br.evolution);
        Ok(br.solid)
    }

    /// Union of multiple solids (sequential fold).
    pub fn fuse_many(
        &mut self,
        solids: &[SolidIdx],
    ) -> Result<SolidIdx, rustkernel_boolean::ops::BooleanError> {
        let _span = info_span!("kernel.fuse_many", count = solids.len()).entered();
        if solids.is_empty() {
            return Err(rustkernel_boolean::ops::BooleanError::DegenerateInput(
                "empty list".into(),
            ));
        }
        let mut result = solids[0];
        for &s in &solids[1..] {
            result = self.fuse(result, s)?;
        }
        // last_evolution holds the evolution from the final fuse step.
        Ok(result)
    }

    /// Intersect a solid with a plane, returning section wire loops.
    pub fn section(
        &mut self,
        solid: SolidIdx,
        plane_origin: [f64; 3],
        plane_normal: [f64; 3],
    ) -> Result<Vec<LoopIdx>, rustkernel_boolean::ops::BooleanError> {
        rustkernel_boolean::section::section_solid(
            &mut self.topo,
            &mut self.geom,
            &self.pipeline,
            solid,
            Point3::new(plane_origin[0], plane_origin[1], plane_origin[2]),
            Vec3::new(plane_normal[0], plane_normal[1], plane_normal[2]),
        )
    }
}

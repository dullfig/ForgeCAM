use rustkernel_math::Point3;
use rustkernel_topology::store::TopoStore;
use rustkernel_topology::topo::*;
use rustkernel_geom::AnalyticalGeomStore;
use tracing::info_span;

use crate::Kernel;

/// Mass properties of a solid body computed from its tessellated mesh.
#[derive(Debug, Clone)]
pub struct MassProperties {
    pub volume: f64,
    pub surface_area: f64,
    pub center_of_mass: Point3,
    /// Inertia tensor about center of mass, row-major \[\[Ixx,Ixy,Ixz\],\[Iyx,Iyy,Iyz\],\[Izx,Izy,Izz\]\]
    pub inertia_tensor: [[f64; 3]; 3],
    pub mass: f64,
}

/// Compute mass properties of a solid via signed-tetrahedra decomposition on its tessellated mesh.
pub fn compute_mass_properties(
    topo: &mut TopoStore,
    geom: &AnalyticalGeomStore,
    solid: SolidIdx,
    density: f64,
) -> MassProperties {
    use rustkernel_topology::tessellate::tessellate_shell;

    let _span = info_span!("compute_mass_properties").entered();

    for &sh in &topo.solids.get(solid).shells.clone() {
        tessellate_shell(topo, sh, geom);
    }

    let mut vol6: f64 = 0.0;
    let mut area: f64 = 0.0;
    let mut cx: f64 = 0.0;
    let mut cy: f64 = 0.0;
    let mut cz: f64 = 0.0;
    // Second moments (diagonal and off-diagonal)
    let mut sum_xx: f64 = 0.0;
    let mut sum_yy: f64 = 0.0;
    let mut sum_zz: f64 = 0.0;
    let mut sum_xy: f64 = 0.0;
    let mut sum_xz: f64 = 0.0;
    let mut sum_yz: f64 = 0.0;

    let faces: Vec<FaceIdx> = topo.solid_faces(solid);
    for face_idx in faces {
        if let Some(mesh) = &topo.faces.get(face_idx).mesh_cache {
            for tri in 0..mesh.triangle_count() {
                let i0 = mesh.indices[tri * 3] as usize;
                let i1 = mesh.indices[tri * 3 + 1] as usize;
                let i2 = mesh.indices[tri * 3 + 2] as usize;
                let a = mesh.positions[i0];
                let b = mesh.positions[i1];
                let c = mesh.positions[i2];

                // Signed volume of tetrahedron (origin, a, b, c) * 6
                let cross = (b - Point3::origin()).cross(&(c - Point3::origin()));
                let det = (a - Point3::origin()).dot(&cross);

                vol6 += det;

                // Surface area
                let ab = b - a;
                let ac = c - a;
                area += 0.5 * ab.cross(&ac).norm();

                // First moments
                cx += det * (a.x + b.x + c.x);
                cy += det * (a.y + b.y + c.y);
                cz += det * (a.z + b.z + c.z);

                // Second moments (diagonal)
                sum_xx += det * (a.x * a.x + b.x * b.x + c.x * c.x + a.x * b.x + a.x * c.x + b.x * c.x);
                sum_yy += det * (a.y * a.y + b.y * b.y + c.y * c.y + a.y * b.y + a.y * c.y + b.y * c.y);
                sum_zz += det * (a.z * a.z + b.z * b.z + c.z * c.z + a.z * b.z + a.z * c.z + b.z * c.z);

                // Second moments (off-diagonal)
                sum_xy += det * (2.0 * a.x * a.y + 2.0 * b.x * b.y + 2.0 * c.x * c.y
                    + a.x * b.y + a.y * b.x + a.x * c.y + a.y * c.x + b.x * c.y + b.y * c.x);
                sum_xz += det * (2.0 * a.x * a.z + 2.0 * b.x * b.z + 2.0 * c.x * c.z
                    + a.x * b.z + a.z * b.x + a.x * c.z + a.z * c.x + b.x * c.z + b.z * c.x);
                sum_yz += det * (2.0 * a.y * a.z + 2.0 * b.y * b.z + 2.0 * c.y * c.z
                    + a.y * b.z + a.z * b.y + a.y * c.z + a.z * c.y + b.y * c.z + b.z * c.y);
            }
        }
    }

    let volume = vol6 / 6.0;
    let com = if vol6.abs() > 1e-30 {
        Point3::new(cx / (4.0 * vol6), cy / (4.0 * vol6), cz / (4.0 * vol6))
    } else {
        Point3::origin()
    };
    let mass = density * volume;

    // Integrals of x^2, y^2, z^2, xy, xz, yz over the solid
    let xx = sum_xx / 60.0;
    let yy = sum_yy / 60.0;
    let zz = sum_zz / 60.0;
    let xy = sum_xy / 120.0;
    let xz = sum_xz / 120.0;
    let yz = sum_yz / 120.0;

    // Inertia about origin
    let ixx_o = density * (yy + zz);
    let iyy_o = density * (xx + zz);
    let izz_o = density * (xx + yy);
    let ixy_o = -density * xy;
    let ixz_o = -density * xz;
    let iyz_o = -density * yz;

    // Parallel axis theorem: shift to center of mass
    let d = com;
    let ixx = ixx_o - mass * (d.y * d.y + d.z * d.z);
    let iyy = iyy_o - mass * (d.x * d.x + d.z * d.z);
    let izz = izz_o - mass * (d.x * d.x + d.y * d.y);
    let ixy = ixy_o + mass * d.x * d.y;
    let ixz = ixz_o + mass * d.x * d.z;
    let iyz = iyz_o + mass * d.y * d.z;

    MassProperties {
        volume,
        surface_area: area,
        center_of_mass: com,
        inertia_tensor: [
            [ixx, ixy, ixz],
            [ixy, iyy, iyz],
            [ixz, iyz, izz],
        ],
        mass,
    }
}

impl Kernel {
    /// Compute mass properties (volume, surface area, center of mass, inertia tensor) for a solid.
    pub fn mass_properties(&mut self, solid: SolidIdx, density: f64) -> MassProperties {
        let _span = info_span!("kernel.mass_properties").entered();
        compute_mass_properties(&mut self.topo, &self.geom, solid, density)
    }
}

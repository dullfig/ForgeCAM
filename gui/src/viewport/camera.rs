//! Orbit/pan/zoom camera controller.
//!
//! Mastercam-style navigation:
//! - Scroll wheel: zoom (dolly toward/away from target)
//! - Middle mouse drag: orbit (rotate around target)
//! - Shift + middle mouse drag: pan (translate target in screen plane)

use nalgebra::{Matrix4, Point3, Vector3};

/// Orbital camera that looks at a target point.
pub struct Camera {
    /// Point the camera orbits around / looks at.
    pub target: Point3<f32>,
    /// Distance from target to eye.
    pub distance: f32,
    /// Azimuth angle in radians (rotation around Y axis).
    pub azimuth: f32,
    /// Elevation angle in radians (rotation up from XZ plane).
    /// Clamped to (-89°, 89°) to avoid gimbal lock at poles.
    pub elevation: f32,
    /// Vertical field of view in radians.
    pub fov: f32,
    /// If true, use orthographic projection instead of perspective.
    pub orthographic: bool,
}

impl Camera {
    pub fn new(target: Point3<f32>, distance: f32, azimuth: f32, elevation: f32) -> Self {
        Self {
            target,
            distance,
            azimuth,
            elevation: elevation.clamp(-1.55, 1.55),
            fov: std::f32::consts::FRAC_PI_4,
            orthographic: false,
        }
    }

    /// Compute the eye position from spherical coordinates.
    pub fn eye(&self) -> Point3<f32> {
        let cos_elev = self.elevation.cos();
        let x = self.distance * cos_elev * self.azimuth.sin();
        let y = self.distance * self.elevation.sin();
        let z = self.distance * cos_elev * self.azimuth.cos();
        Point3::new(
            self.target.x + x,
            self.target.y + y,
            self.target.z + z,
        )
    }

    /// Build the view matrix (look-at from eye to target).
    pub fn view_matrix(&self) -> Matrix4<f32> {
        let eye = self.eye();
        let up = Vector3::new(0.0, 1.0, 0.0);
        Matrix4::look_at_rh(&eye, &self.target, &up)
    }

    /// Build the projection matrix (perspective or orthographic).
    pub fn projection_matrix(&self, aspect: f32) -> Matrix4<f32> {
        if self.orthographic {
            // Orthographic half-extents scale with distance so zoom feels the same.
            let half_h = self.distance * (self.fov * 0.5).tan();
            let half_w = half_h * aspect;
            Matrix4::new_orthographic(-half_w, half_w, -half_h, half_h, 0.01, 1000.0)
        } else {
            Matrix4::new_perspective(aspect, self.fov, 0.01, 1000.0)
        }
    }

    /// Build the combined MVP matrix (no model transform — identity).
    pub fn mvp(&self, aspect: f32) -> Matrix4<f32> {
        self.projection_matrix(aspect) * self.view_matrix()
    }

    /// Orbit: rotate azimuth and elevation by delta angles.
    pub fn orbit(&mut self, delta_azimuth: f32, delta_elevation: f32) {
        self.azimuth += delta_azimuth;
        self.elevation = (self.elevation + delta_elevation).clamp(-1.55, 1.55);
    }

    /// Pan: translate target in the screen plane.
    pub fn pan(&mut self, delta_x: f32, delta_y: f32) {
        let view = self.view_matrix();
        // Screen-right and screen-up in world space are rows of the view matrix.
        let right = Vector3::new(view[(0, 0)], view[(0, 1)], view[(0, 2)]);
        let up = Vector3::new(view[(1, 0)], view[(1, 1)], view[(1, 2)]);

        // Scale pan by distance so it feels consistent at any zoom level.
        let scale = self.distance * 0.002;
        self.target += right * (-delta_x * scale) + up * (delta_y * scale);
    }

    /// Zoom: adjust distance (dolly toward/away from target).
    /// `delta` is raw scroll pixels (typically ±120 per notch on Windows).
    pub fn zoom(&mut self, delta: f32) {
        // Multiplicative zoom so it feels consistent at any distance.
        // Scale down heavily — scroll deltas are large pixel values.
        let factor = 1.0 - delta * 0.001;
        self.distance = (self.distance * factor).clamp(0.01, 10000.0);
    }
}

impl Default for Camera {
    fn default() -> Self {
        Self::new(
            Point3::new(1.5, 0.5, 0.5),  // centered on typical demo scene
            8.0,                           // distance
            0.4,                           // azimuth (slight angle)
            0.35,                          // elevation (looking slightly down)
        )
    }
}

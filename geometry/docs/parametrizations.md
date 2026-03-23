# Surface & Curve Parametrizations

These parametrizations are load-bearing. The kernel's inverse UV mapping, face
tessellation, and boolean pipeline all depend on these exact conventions.

## Curves

All analytical curves use **t in [0, 1]** as the parameter domain.

| Curve | t = 0 | t = 1 | Mapping |
|-------|-------|-------|---------|
| `LineSegment` | `start` | `end` | Linear interpolation |
| `CircleCurve` | `center + radius * ref_dir` | Same (full loop) | `t * 2pi` |
| `CircularArcCurve` | At `start_angle` | At `end_angle` | `start + t * (end - start)` |
| `EllipseCurve` | `center + semi_major * major_dir` | Same (full loop) | `t * 2pi` |
| NURBS | Knot domain start | Knot domain end | Native knot space |

### Circle/Arc geometry

Position at angle theta:
```
P = center + radius * (cos(theta) * ref_dir + sin(theta) * (axis x ref_dir))
```

`axis` is the unit normal to the circle plane. `ref_dir` is the unit vector defining
angle = 0. They must be perpendicular.

### Ellipse geometry

```
P = center + semi_major * cos(theta) * major_dir
           + semi_minor * sin(theta) * (axis x major_dir)
```

## Surfaces

### Plane

- **eval(u, v)**: `origin + u * basis_u + v * basis_v`
  where `(basis_u, basis_v) = orthonormal_basis(normal)`
- **normal**: constant everywhere
- **domain**: unbounded
- **inverse_uv**: project point onto plane, decompose into basis

### Cylinder

- **u** = angle `[0, 2pi]`
- **v** = axial distance (unbounded)
- **eval(u, v)**: `origin + v*axis + |radius| * (cos(u)*u_dir + sin(u)*v_dir)`
  where `(u_dir, v_dir) = orthonormal_basis(axis)`
- **normal(u, v)**: `sign(radius) * (cos(u)*u_dir + sin(u)*v_dir)`
- **inverse_uv**: project onto axis for v, atan2 of lateral component for u

### Sphere

- **u** = azimuth `[0, 2pi]`
- **v** = elevation `[-pi/2, pi/2]`
- **Polar axis**: hardcoded to z-axis `(0, 0, 1)`
- **eval(u, v)**: `center + |radius| * (cos(v)*cos(u), cos(v)*sin(u), sin(v))`
- **normal(u, v)**: `sign(radius) * (cos(v)*cos(u), cos(v)*sin(u), sin(v))`
- **inverse_uv**: `u = atan2(dy, dx)`, `v = asin(dz / |d|)`

### Cone

- **u** = angle `[0, 2pi]`
- **v** = distance along axis from apex `[0, inf)`
- **eval(u, v)**: `apex + v*axis + v*tan(|alpha|) * (cos(u)*u_dir + sin(u)*v_dir)`
- **normal(u, v)**: `sign(half_angle) * (cos(|alpha|)*radial - sin(|alpha|)*axis)`
- **inverse_uv**: project onto axis for v, atan2 of lateral for u

### Torus

- **u** = major angle `[0, 2pi]` (around the axis of revolution)
- **v** = minor angle `[0, 2pi]` (around the tube cross-section)
- **eval(u, v)**:
  ```
  major_dir = cos(u)*u_dir + sin(u)*v_dir
  P = center + (R + r*cos(v)) * major_dir + r*sin(v) * axis
  ```
  where `R = major_radius`, `r = |minor_radius|`
- **normal(u, v)**: `sign(minor_radius) * (cos(v)*major_dir + sin(v)*axis)`

### NURBS (FlippableNurbs)

- **u, v** in native knot domain (from `knots_domain()`)
- **eval**: delegates to curvo `point_at(u, v)`
- **normal**: delegates to curvo `normal_at(u, v)`, flipped if `flipped == true`
- **inverse_uv**: 8x8 grid search + 5 Gauss-Newton iterations (finite-difference Jacobian, step 1e-6)

## Negative radius convention

This is the most important convention in the crate. **Do not change it.**

| Surface | Positive | Negative |
|---------|----------|----------|
| Cylinder | Normal points outward (away from axis) | Normal points inward (toward axis) |
| Sphere | Normal points outward (away from center) | Normal points inward (toward center) |
| Torus | Normal points outward (away from tube center) | Normal points inward |
| Cone | `half_angle > 0`: normal away from axis | `half_angle < 0`: normal toward axis |

The kernel uses this to distinguish between the inside and outside of a solid's
faces. A face on the outside of a cylinder has positive radius; a hole through a
solid has negative radius on the cylinder surface.

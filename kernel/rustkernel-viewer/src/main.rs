mod convert;

use convert::{
    build_edge_mesh, compute_silhouette_lines, extract_shell_edges, extract_silhouette_candidates,
    merge_shell_meshes, SilhouetteCandidate,
};
use rustkernel_kernel::Kernel;
use rustkernel_math::{Point3, Vec3};
use rustkernel_topology::tessellate::TessellationOptions;
use rustkernel_topology::topo::SolidIdx;
use three_d::*;

// ---------------------------------------------------------------------------
// Quality presets
// ---------------------------------------------------------------------------

const QUALITY_NAMES: &[&str] = &["Low", "Medium", "High", "Ultra"];

fn quality_options(level: usize) -> TessellationOptions {
    match level {
        0 => TessellationOptions {
            min_segments: 16,
            angular_tolerance: std::f64::consts::PI / 16.0,
            chordal_tolerance: 0.01,
        },
        1 => TessellationOptions {
            min_segments: 24,
            angular_tolerance: std::f64::consts::PI / 32.0,
            chordal_tolerance: 0.002,
        },
        2 => TessellationOptions {
            min_segments: 32,
            angular_tolerance: std::f64::consts::PI / 64.0,
            chordal_tolerance: 0.0005,
        },
        _ => TessellationOptions {
            min_segments: 64,
            angular_tolerance: std::f64::consts::PI / 128.0,
            chordal_tolerance: 0.0001,
        },
    }
}

// ---------------------------------------------------------------------------
// Display modes
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum DisplayMode {
    Shaded,
    ShadedWithEdges,
    HiddenLine,
    Wireframe,
}

const MODE_NAMES: &[&str] = &["Shaded", "Shaded+Edges", "Hidden Line", "Wireframe"];

impl DisplayMode {
    fn index(self) -> usize {
        match self {
            Self::Shaded => 0,
            Self::ShadedWithEdges => 1,
            Self::HiddenLine => 2,
            Self::Wireframe => 3,
        }
    }
    fn next(self) -> Self {
        match self {
            Self::Shaded => Self::ShadedWithEdges,
            Self::ShadedWithEdges => Self::HiddenLine,
            Self::HiddenLine => Self::Wireframe,
            Self::Wireframe => Self::Shaded,
        }
    }
}

// ---------------------------------------------------------------------------
// Scene definitions — each builds geometry via Kernel and returns a SolidIdx
// ---------------------------------------------------------------------------

struct Scene {
    name: &'static str,
    build: fn(&mut Kernel) -> SolidIdx,
}

fn scene_box(k: &mut Kernel) -> SolidIdx {
    k.make_box(2.0, 1.5, 1.0)
}

fn scene_cylinder(k: &mut Kernel) -> SolidIdx {
    k.make_cylinder(1.0, 2.0)
}

fn scene_sphere(k: &mut Kernel) -> SolidIdx {
    k.make_sphere(1.2)
}

fn scene_cone(k: &mut Kernel) -> SolidIdx {
    k.make_cone(1.0, 0.3, 2.0)
}

fn scene_pointed_cone(k: &mut Kernel) -> SolidIdx {
    k.make_cone(1.2, 0.0, 2.0)
}

fn scene_torus(k: &mut Kernel) -> SolidIdx {
    k.make_torus(1.5, 0.4)
}

fn scene_chamfered_box(k: &mut Kernel) -> SolidIdx {
    let solid = k.make_box(2.0, 2.0, 1.5);
    let edges = collect_edges(k, solid);
    // Chamfer one edge to avoid SharedVertex errors on adjacent edges
    let edge = &edges[0..1];
    match k.chamfer_edges(solid, edge, 0.2) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Chamfer failed: {:?}, showing unchamfered", e);
            solid
        }
    }
}

fn scene_filleted_box(k: &mut Kernel) -> SolidIdx {
    let solid = k.make_box(2.0, 2.0, 1.5);
    let edges = collect_edges(k, solid);
    // Fillet one edge at a time to avoid SharedVertex errors on adjacent edges
    let edge = &edges[0..1];
    match k.fillet_edges(solid, edge, 0.3) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Fillet failed: {:?}, showing unfilleted", e);
            solid
        }
    }
}

fn scene_boolean_cut(k: &mut Kernel) -> SolidIdx {
    let block = k.make_box(2.0, 2.0, 2.0);
    let tool = k.make_box_at([1.0, 0.0, 0.0], 2.0, 2.0, 2.0);
    match k.cut(block, tool) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Boolean cut failed: {:?}", e);
            block
        }
    }
}

fn scene_boolean_fuse(k: &mut Kernel) -> SolidIdx {
    let block = k.make_box(2.0, 1.5, 1.0);
    let block2 = k.make_box_at([1.0, 0.0, 0.0], 1.5, 1.5, 1.0);
    match k.fuse(block, block2) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Boolean fuse failed: {:?}", e);
            block
        }
    }
}

fn scene_boolean_common(k: &mut Kernel) -> SolidIdx {
    let block = k.make_box(2.0, 2.0, 2.0);
    let block2 = k.make_box_at([1.0, 0.0, 0.0], 2.0, 2.0, 2.0);
    match k.common(block, block2) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Boolean common failed: {:?}", e);
            block
        }
    }
}

fn scene_interior_pocket(k: &mut Kernel) -> SolidIdx {
    // Rectangular pocket punched through a plate — matches test_cut_single_interior_pocket
    // Block: 4x4x1, pocket: 0.8x0.8 at boundary corner, punches through Z
    let plate = k.make_box_at([-2.0, -2.0, -0.5], 4.0, 4.0, 1.0);
    let pocket = k.make_box_at([-0.4, -0.4, -0.6], 0.8, 0.8, 1.5);
    match k.cut(plate, pocket) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Interior pocket cut failed: {:?}", e);
            plate
        }
    }
}

fn scene_extrude_l(k: &mut Kernel) -> SolidIdx {
    let profile = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 0.8, 0.0),
        Point3::new(0.8, 0.8, 0.0),
        Point3::new(0.8, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    k.extrude(&profile, Vec3::new(0.0, 0.0, 1.0), 1.0)
}

fn scene_revolve(k: &mut Kernel) -> SolidIdx {
    // Profile in XZ meridional plane (y=0): x = radius, z = height
    let profile = vec![
        Point3::new(0.5, 0.0, 0.0),
        Point3::new(1.2, 0.0, 0.0),
        Point3::new(1.2, 0.0, 1.0),
        Point3::new(0.5, 0.0, 1.0),
    ];
    k.revolve(
        &profile,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::TAU,
    )
}

fn scene_dissimilar_fillet(k: &mut Kernel) -> SolidIdx {
    use rustkernel_topology::topo::EdgeIdx;
    use std::collections::HashSet;
    let solid = k.make_box(2.0, 2.0, 2.0);
    let target = Point3::new(0.0, 0.0, 0.0);

    // Collect all unique edges, find those incident on the corner vertex near target.
    let shell_idx = k.topo().solids.get(solid).outer_shell();
    let faces = k.topo().shells.get(shell_idx).faces.clone();
    let mut seen = HashSet::new();
    let mut corner_edges: Vec<EdgeIdx> = Vec::new();
    for fi in &faces {
        let start = k.topo().loops.get(k.topo().faces.get(*fi).outer_loop).half_edge;
        let mut cur = start;
        loop {
            let he = k.topo().half_edges.get(cur);
            if seen.insert(he.edge.raw()) {
                let pa = k.geom().points[k.topo().vertices.get(he.origin).point_id as usize];
                let pb_vid = if let Some(twin) = he.twin {
                    k.topo().half_edges.get(twin).origin
                } else { he.origin };
                let pb = k.geom().points[k.topo().vertices.get(pb_vid).point_id as usize];
                if (pa - target).norm() < 1e-6 || (pb - target).norm() < 1e-6 {
                    corner_edges.push(he.edge);
                }
            }
            cur = he.next;
            if cur == start { break; }
        }
    }

    if corner_edges.len() >= 3 {
        match k.euler_fillet_edges_with_radii(solid, &corner_edges[..3], &[0.2, 0.3, 0.4]) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Dissimilar fillet failed: {:?}", e);
                solid
            }
        }
    } else {
        eprintln!("Could not find 3 corner edges, found {}", corner_edges.len());
        solid
    }
}

fn scene_multi_boolean(k: &mut Kernel) -> SolidIdx {
    // Fuse two offset boxes into an L-shaped solid — demonstrates multi-step boolean
    let block_a = k.make_box(2.0, 1.5, 1.0);
    let block_b = k.make_box_at([1.0, 0.0, 0.0], 1.5, 1.5, 1.0);
    let block_c = k.make_box_at([0.0, 1.0, 0.0], 1.5, 1.5, 1.0);
    let mut result = match k.fuse(block_a, block_b) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Fuse A+B failed: {:?}", e);
            return block_a;
        }
    };
    result = match k.fuse(result, block_c) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Fuse (A+B)+C failed: {:?}", e);
            result
        }
    };
    result
}

const SCENES: &[Scene] = &[
    Scene { name: "Box", build: scene_box },
    Scene { name: "Cylinder", build: scene_cylinder },
    Scene { name: "Sphere", build: scene_sphere },
    Scene { name: "Cone (frustum)", build: scene_cone },
    Scene { name: "Cone (pointed)", build: scene_pointed_cone },
    Scene { name: "Torus", build: scene_torus },
    Scene { name: "Chamfered box", build: scene_chamfered_box },
    Scene { name: "Filleted box", build: scene_filleted_box },
    Scene { name: "Extrude (L-shape)", build: scene_extrude_l },
    Scene { name: "Revolve (ring)", build: scene_revolve },
    Scene { name: "Boolean cut (box - box)", build: scene_boolean_cut },
    Scene { name: "Boolean fuse (box + box)", build: scene_boolean_fuse },
    Scene { name: "Boolean common (box & box)", build: scene_boolean_common },
    Scene { name: "Cut: interior pocket", build: scene_interior_pocket },
    Scene { name: "Multi-fuse (3 boxes)", build: scene_multi_boolean },
    Scene { name: "Dissimilar fillet (3 radii)", build: scene_dissimilar_fillet },
];

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Collect all unique edge indices from a solid.
fn collect_edges(
    k: &Kernel,
    solid: SolidIdx,
) -> Vec<rustkernel_topology::topo::EdgeIdx> {
    let shell_idx = k.topo().solids.get(solid).outer_shell();
    let faces = k.topo().shells.get(shell_idx).faces.clone();
    let mut edges = Vec::new();
    for fi in &faces {
        let face = k.topo().faces.get(*fi);
        let start = k.topo().loops.get(face.outer_loop).half_edge;
        let mut cur = start;
        loop {
            let he = k.topo().half_edges.get(cur);
            if !edges.contains(&he.edge) {
                edges.push(he.edge);
            }
            cur = he.next;
            if cur == start { break; }
        }
    }
    edges
}

/// Scene data: model mesh + edge line mesh + silhouette candidates.
struct SceneData {
    name: String,
    model_positions: Vec<Vector3<f32>>,
    model_normals: Vec<Vector3<f32>>,
    model_indices: Vec<u32>,
    edge_positions: Vec<Vector3<f32>>,
    edge_normals: Vec<Vector3<f32>>,
    edge_indices: Vec<u32>,
    silhouette_candidates: Vec<SilhouetteCandidate>,
}

/// Build a scene: create kernel, build geometry, tessellate, extract model + edges.
fn build_scene(scene_idx: usize, quality: usize) -> SceneData {
    let scene = &SCENES[scene_idx];
    let mut kernel = Kernel::new();
    let solid = (scene.build)(&mut kernel);

    let opts = quality_options(quality);
    kernel.tessellate_solid_with_options(solid, &opts);

    let shell_idx = kernel.topo().solids.get(solid).outer_shell();

    // Model mesh
    let merged = merge_shell_meshes(kernel.topo(), shell_idx);
    let model_positions = merged.positions.iter().map(|p| Vector3::new(p[0], p[1], p[2])).collect();
    let model_normals = merged.normals.iter().map(|n| Vector3::new(n[0], n[1], n[2])).collect();

    // Silhouette candidates (from tessellation mesh edges)
    let silhouette_candidates = extract_silhouette_candidates(kernel.topo(), shell_idx);

    // Edge mesh (sharp B-Rep edges only)
    let edge_lines = extract_shell_edges(kernel.topo(), kernel.geom(), shell_idx);
    let edge_merged = build_edge_mesh(&edge_lines);
    let edge_positions = edge_merged.positions.iter().map(|p| Vector3::new(p[0], p[1], p[2])).collect();
    let edge_normals = edge_merged.normals.iter().map(|n| Vector3::new(n[0], n[1], n[2])).collect();

    println!(
        "[{}] {}: {} verts, {} tris, {} edges, {} silhouette candidates",
        QUALITY_NAMES[quality], scene.name,
        merged.positions.len(), merged.indices.len() / 3,
        edge_lines.len(), silhouette_candidates.len(),
    );

    SceneData {
        name: scene.name.to_string(),
        model_positions,
        model_normals,
        model_indices: merged.indices,
        edge_positions,
        edge_normals,
        edge_indices: edge_merged.indices,
        silhouette_candidates,
    }
}

/// Create GPU model from scene data.
fn make_model(ctx: &Context, data: &SceneData) -> Gm<Mesh, PhysicalMaterial> {
    Gm::new(
        Mesh::new(ctx, &CpuMesh {
            positions: Positions::F32(data.model_positions.clone()),
            normals: Some(data.model_normals.clone()),
            indices: Indices::U32(data.model_indices.clone()),
            ..Default::default()
        }),
        PhysicalMaterial::new_opaque(
            ctx,
            &CpuMaterial {
                albedo: Srgba::new(120, 160, 220, 255),
                roughness: 0.7,
                metallic: 0.0,
                ..Default::default()
            },
        ),
    )
}

/// Create GPU model for hidden-line mode (flat white, writes depth only).
fn make_hidden_model(ctx: &Context, data: &SceneData) -> Gm<Mesh, ColorMaterial> {
    Gm::new(
        Mesh::new(ctx, &CpuMesh {
            positions: Positions::F32(data.model_positions.clone()),
            normals: Some(data.model_normals.clone()),
            indices: Indices::U32(data.model_indices.clone()),
            ..Default::default()
        }),
        ColorMaterial {
            color: Srgba::new(200, 200, 205, 255),
            ..Default::default()
        },
    )
}

/// Create GPU edge overlay mesh (black, no backface culling).
fn make_edge_model_from_parts(
    ctx: &Context,
    positions: &[Vector3<f32>],
    normals: &[Vector3<f32>],
    indices: &[u32],
) -> Gm<Mesh, ColorMaterial> {
    Gm::new(
        Mesh::new(ctx, &CpuMesh {
            positions: Positions::F32(positions.to_vec()),
            normals: Some(normals.to_vec()),
            indices: Indices::U32(indices.to_vec()),
            ..Default::default()
        }),
        ColorMaterial {
            color: Srgba::new(20, 20, 20, 255),
            render_states: RenderStates {
                cull: Cull::None,
                ..Default::default()
            },
            ..Default::default()
        },
    )
}

/// Create GPU edge overlay mesh from SceneData.
fn make_edge_model(ctx: &Context, data: &SceneData) -> Gm<Mesh, ColorMaterial> {
    make_edge_model_from_parts(ctx, &data.edge_positions, &data.edge_normals, &data.edge_indices)
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // --list: print scenes and exit
    if args.iter().any(|a| a == "--list") {
        println!("Available scenes:");
        for (i, s) in SCENES.iter().enumerate() {
            println!("  {:2}: {}", i, s.name);
        }
        return;
    }

    // Optional scene index as first arg (default 0)
    let scene_idx: usize = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0)
        .min(SCENES.len() - 1);

    let mut quality: usize = 1; // Start at Medium
    let mut display_mode = DisplayMode::Shaded;

    println!("Scene {}: {}", scene_idx, SCENES[scene_idx].name);
    println!("Quality: {}", QUALITY_NAMES[quality]);
    println!("MMB: orbit | Shift+MMB: pan | Scroll: zoom | Left/Right: scene | Q: quality | D: display mode");

    let window = Window::new(WindowSettings {
        title: format!("ForgeCAM Viewer — {}", SCENES[scene_idx].name),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();

    // Gradient background quad (light grey top → darker grey bottom)
    let bg_quad = Gm::new(
        Mesh::new(
            &context,
            &CpuMesh {
                positions: Positions::F32(vec![
                    Vector3::new(-10.0, -1.0, 0.0),
                    Vector3::new(10.0, -1.0, 0.0),
                    Vector3::new(10.0, 1.0, 0.0),
                    Vector3::new(-10.0, 1.0, 0.0),
                ]),
                indices: Indices::U32(vec![0, 1, 2, 0, 2, 3]),
                colors: Some(vec![
                    Srgba::new(140, 140, 145, 255),
                    Srgba::new(140, 140, 145, 255),
                    Srgba::new(220, 220, 225, 255),
                    Srgba::new(220, 220, 225, 255),
                ]),
                ..Default::default()
            },
        ),
        ColorMaterial {
            render_states: RenderStates {
                depth_test: DepthTest::Always,
                write_mask: WriteMask::COLOR,
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let bg_camera = Camera::new_orthographic(
        window.viewport(),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        2.0,
        0.1,
        10.0,
    );

    // Build initial scene
    let data = build_scene(scene_idx, quality);
    let mut model = make_model(&context, &data);
    let mut hidden_model = make_hidden_model(&context, &data);
    let mut edge_model = make_edge_model(&context, &data);
    let mut silhouette_candidates = data.silhouette_candidates;
    let mut silhouette_model: Option<Gm<Mesh, ColorMaterial>> = None;
    let mut last_eye_pos: Vector3<f32> = vec3(0.0, 0.0, 0.0);

    let mut camera = Camera::new_perspective(
        window.viewport(),
        vec3(4.0, 4.0, 3.0),
        vec3(0.0, 0.0, 0.5),
        vec3(0.0, 0.0, 1.0),
        degrees(45.0),
        0.1,
        100.0,
    );

    let orbit_target = vec3(0.0, 0.0, 0.5);
    let min_dist = 0.5_f32;
    let max_dist = 50.0_f32;

    // Camera-following headlamp + fill light + ambient
    let mut headlamp = DirectionalLight::new(&context, 1.8, Srgba::WHITE, vec3(0.0, 0.0, -1.0));
    let mut fill = DirectionalLight::new(&context, 0.6, Srgba::new(200, 210, 230, 255), vec3(0.0, 0.0, -1.0));
    let ambient = AmbientLight::new(&context, 0.35, Srgba::WHITE);

    let mut current = scene_idx;

    window.render_loop(move |mut frame_input| {
        camera.set_viewport(frame_input.viewport);

        // Update headlamp to follow camera view direction
        let view_dir = camera.view_direction();
        let right = camera.right_direction();
        let up_dir = right.cross(view_dir);
        headlamp.direction = view_dir + 0.15 * right - 0.1 * up_dir;
        fill.direction = -view_dir + 0.5 * right + 0.3 * up_dir;

        let mut needs_rebuild = false;

        for event in frame_input.events.iter_mut() {
            match event {
                Event::MouseMotion {
                    button,
                    delta,
                    modifiers,
                    handled,
                    ..
                } => {
                    if *button == Some(MouseButton::Middle) && !*handled {
                        if modifiers.shift {
                            let speed = 0.005
                                * (camera.position() - orbit_target).magnitude();
                            let right = camera.right_direction();
                            let up = camera.up();
                            let pan = right * (-delta.0 * speed)
                                + up * (delta.1 * speed);
                            camera.translate(pan);
                        } else {
                            camera.rotate_around_with_fixed_up(
                                orbit_target,
                                0.01 * delta.0,
                                0.01 * delta.1,
                            );
                        }
                        *handled = true;
                    }
                }
                Event::MouseWheel {
                    delta,
                    handled,
                    ..
                } => {
                    if !*handled {
                        let dist = (camera.position() - orbit_target).magnitude();
                        let speed = 0.01 * dist + 0.001;
                        camera.zoom_towards(
                            orbit_target,
                            speed * delta.1,
                            min_dist,
                            max_dist,
                        );
                        *handled = true;
                    }
                }
                Event::KeyPress { kind, handled, .. } => {
                    if !*handled {
                        match kind {
                            Key::ArrowRight => {
                                current = (current + 1) % SCENES.len();
                                needs_rebuild = true;
                                *handled = true;
                            }
                            Key::ArrowLeft => {
                                current = if current == 0 {
                                    SCENES.len() - 1
                                } else {
                                    current - 1
                                };
                                needs_rebuild = true;
                                *handled = true;
                            }
                            Key::Q => {
                                quality = (quality + 1) % QUALITY_NAMES.len();
                                needs_rebuild = true;
                                *handled = true;
                            }
                            Key::D => {
                                display_mode = display_mode.next();
                                println!("Display: {}", MODE_NAMES[display_mode.index()]);
                                *handled = true;
                            }
                            _ => {}
                        }
                    }
                }
                _ => {}
            }
        }

        if needs_rebuild {
            let data = build_scene(current, quality);
            model = make_model(&context, &data);
            hidden_model = make_hidden_model(&context, &data);
            edge_model = make_edge_model(&context, &data);
            silhouette_candidates = data.silhouette_candidates;
            silhouette_model = None;
            last_eye_pos = vec3(0.0, 0.0, 0.0); // force silhouette rebuild
        }

        // Rebuild silhouette edges when camera moves (for HiddenLine mode).
        let eye = camera.position();
        if display_mode == DisplayMode::HiddenLine && (eye - last_eye_pos).magnitude() > 1e-4 {
            last_eye_pos = eye;
            let sil_lines = compute_silhouette_lines(
                &silhouette_candidates,
                [eye.x, eye.y, eye.z],
            );
            if !sil_lines.is_empty() {
                let sil_mesh = build_edge_mesh(&sil_lines);
                let positions: Vec<Vector3<f32>> = sil_mesh.positions.iter().map(|p| Vector3::new(p[0], p[1], p[2])).collect();
                let normals: Vec<Vector3<f32>> = sil_mesh.normals.iter().map(|n| Vector3::new(n[0], n[1], n[2])).collect();
                silhouette_model = Some(make_edge_model_from_parts(
                    &context, &positions, &normals, &sil_mesh.indices,
                ));
            } else {
                silhouette_model = None;
            }
        }

        model.set_transformation(Mat4::identity());
        hidden_model.set_transformation(Mat4::identity());
        edge_model.set_transformation(Mat4::identity());
        if let Some(ref mut sm) = silhouette_model {
            sm.set_transformation(Mat4::identity());
        }

        let binding = frame_input.screen();
        let screen = binding
            .clear(ClearState::color_and_depth(0.0, 0.0, 0.0, 1.0, 1.0))
            .render(&bg_camera, &bg_quad, &[]);

        match display_mode {
            DisplayMode::Shaded => {
                screen.render(&camera, &model, &[&headlamp, &fill, &ambient]);
            }
            DisplayMode::ShadedWithEdges => {
                screen
                    .render(&camera, &model, &[&headlamp, &fill, &ambient])
                    .render(&camera, &edge_model, &[]);
            }
            DisplayMode::HiddenLine => {
                let s = screen
                    .render(&camera, &hidden_model, &[])
                    .render(&camera, &edge_model, &[]);
                if let Some(ref sm) = silhouette_model {
                    s.render(&camera, sm, &[]);
                }
            }
            DisplayMode::Wireframe => {
                screen.render(&camera, &edge_model, &[]);
            }
        }

        FrameOutput::default()
    });
}

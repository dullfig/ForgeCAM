mod convert;

use convert::merge_shell_meshes;
use rustkernel_primitives::box_builder::make_box;
use rustkernel_topology::tessellate::tessellate_shell;
use three_d::*;

fn main() {
    let window = Window::new(WindowSettings {
        title: "RustKernel — Box on Screen".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();

    // Build a box and tessellate it.
    let (mut topo, geom, solid_idx) = make_box(2.0, 1.5, 1.0);
    let shell_idx = topo.solids.get(solid_idx).shell;
    tessellate_shell(&mut topo, shell_idx, &geom);

    // Merge face meshes into a single GPU-ready mesh.
    let merged = merge_shell_meshes(&topo, shell_idx);

    // Build three-d CpuMesh from raw data.
    let positions: Vec<Vector3<f32>> = merged
        .positions
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();
    let normals: Vec<Vector3<f32>> = merged
        .normals
        .iter()
        .map(|n| Vector3::new(n[0], n[1], n[2]))
        .collect();
    let indices_flat: Vec<u32> = merged.indices.clone();

    let cpu_mesh = CpuMesh {
        positions: Positions::F32(positions),
        normals: Some(normals),
        indices: Indices::U32(indices_flat),
        ..Default::default()
    };

    let mut model = Gm::new(
        Mesh::new(&context, &cpu_mesh),
        PhysicalMaterial::new_opaque(
            &context,
            &CpuMaterial {
                albedo: Srgba::new(100, 149, 237, 255), // cornflower blue
                roughness: 0.5,
                metallic: 0.0,
                ..Default::default()
            },
        ),
    );

    let mut camera = Camera::new_perspective(
        window.viewport(),
        vec3(3.0, 3.0, 3.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 0.0, 1.0),
        degrees(45.0),
        0.1,
        100.0,
    );

    let mut orbit = OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 10.0);

    let light0 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, vec3(-1.0, -1.0, -1.0));
    let light1 = DirectionalLight::new(&context, 1.0, Srgba::WHITE, vec3(1.0, 1.0, 1.0));
    let ambient = AmbientLight::new(&context, 0.3, Srgba::WHITE);

    window.render_loop(move |mut frame_input| {
        camera.set_viewport(frame_input.viewport);
        orbit.handle_events(&mut camera, &mut frame_input.events);

        model.set_transformation(Mat4::identity());

        frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.15, 0.15, 0.15, 1.0, 1.0))
            .render(&camera, &model, &[&light0, &light1, &ambient]);

        FrameOutput::default()
    });
}

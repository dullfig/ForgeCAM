slint::include_modules!();

mod viewport;

use slint::wgpu_28::wgpu;
use slint::wgpu_28::{WGPUConfiguration, WGPUSettings};
use viewport::Vertex;

// ── Gradient uniforms ────────────────────────────────────────────────

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct GradientUniforms {
    color_top: [f32; 4],
    color_bot: [f32; 4],
}

impl GradientUniforms {
    fn from_hex(top: u32, bot: u32) -> Self {
        Self {
            color_top: [
                ((top >> 16) & 0xFF) as f32 / 255.0,
                ((top >> 8) & 0xFF) as f32 / 255.0,
                (top & 0xFF) as f32 / 255.0,
                1.0,
            ],
            color_bot: [
                ((bot >> 16) & 0xFF) as f32 / 255.0,
                ((bot >> 8) & 0xFF) as f32 / 255.0,
                (bot & 0xFF) as f32 / 255.0,
                1.0,
            ],
        }
    }
}

// Vertex layout constant for wgpu pipeline setup.
const VERTEX_LAYOUT: wgpu::VertexBufferLayout<'static> = wgpu::VertexBufferLayout {
    array_stride: std::mem::size_of::<Vertex>() as u64,
    step_mode: wgpu::VertexStepMode::Vertex,
    attributes: &[
        wgpu::VertexAttribute {
            offset: 0,
            shader_location: 0,
            format: wgpu::VertexFormat::Float32x3,
        },
        wgpu::VertexAttribute {
            offset: 12,
            shader_location: 1,
            format: wgpu::VertexFormat::Float32x3,
        },
    ],
};

// ── Scene uniforms (MVP + light) ─────────────────────────────────────

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct SceneUniforms {
    mvp: [[f32; 4]; 4],
    model: [[f32; 4]; 4],
    light_dir: [f32; 4],    // normalized, world space
    base_color: [f32; 4],   // object color
}

fn mat4_to_array(m: &nalgebra::Matrix4<f32>) -> [[f32; 4]; 4] {
    // nalgebra is column-major, WGSL expects column-major — direct copy
    let s = m.as_slice();
    [
        [s[0], s[1], s[2], s[3]],
        [s[4], s[5], s[6], s[7]],
        [s[8], s[9], s[10], s[11]],
        [s[12], s[13], s[14], s[15]],
    ]
}

fn build_scene_uniforms(time: f32, aspect: f32) -> SceneUniforms {
    use nalgebra::{Matrix4, Point3, Vector3};

    // Slow rotation
    let angle = time * 0.4;

    let model = Matrix4::from_euler_angles(0.3, angle, 0.1);

    // Camera positioned for real-world dimensions (kernel uses inches).
    let eye = Point3::new(3.0, 4.0, 8.0);
    let target = Point3::new(1.5, 0.5, 0.5);
    let up = Vector3::new(0.0, 1.0, 0.0);
    let view = Matrix4::look_at_rh(&eye, &target, &up);

    let fov = std::f32::consts::FRAC_PI_4;
    let near = 0.1;
    let far = 100.0;
    let proj = Matrix4::new_perspective(aspect, fov, near, far);

    let mvp = proj * view * model;

    // Light from upper-left-front
    let light = Vector3::new(0.5, 0.8, 0.6).normalize();

    SceneUniforms {
        mvp: mat4_to_array(&mvp),
        model: mat4_to_array(&model),
        light_dir: [light.x, light.y, light.z, 0.0],
        base_color: [0.35, 0.55, 0.75, 1.0], // steel blue
    }
}

// ── Viewport renderer ────────────────────────────────────────────────

struct ViewportRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    // Gradient
    gradient_pipeline: wgpu::RenderPipeline,
    gradient_uniform_buf: wgpu::Buffer,
    gradient_bind_group: wgpu::BindGroup,
    // Model (was "cube" — now renders kernel solids)
    model_pipeline: wgpu::RenderPipeline,
    model_vertex_buf: wgpu::Buffer,
    model_index_buf: wgpu::Buffer,
    model_index_count: u32,
    scene_uniform_buf: wgpu::Buffer,
    scene_bind_group: wgpu::BindGroup,
    // Time
    start_time: std::time::Instant,
}

impl ViewportRenderer {
    fn new(device: &wgpu::Device, queue: &wgpu::Queue, mesh_vertices: &[Vertex], mesh_indices: &[u32]) -> Self {
        let color_format = wgpu::TextureFormat::Rgba8UnormSrgb;
        let depth_format = wgpu::TextureFormat::Depth32Float;

        // ── Gradient pipeline ────────────────────────────────────
        let grad_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("gradient_shader"),
            source: wgpu::ShaderSource::Wgsl(GRADIENT_SHADER.into()),
        });

        let grad_uniforms =
            GradientUniforms::from_hex(DEFAULT_GRADIENT_TOP, DEFAULT_GRADIENT_BOT);
        let gradient_uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gradient_uniforms"),
            size: std::mem::size_of::<GradientUniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&gradient_uniform_buf, 0, bytemuck::bytes_of(&grad_uniforms));

        let grad_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("gradient_bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let gradient_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("gradient_bg"),
            layout: &grad_bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: gradient_uniform_buf.as_entire_binding(),
            }],
        });

        let grad_pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("gradient_pl"),
            bind_group_layouts: &[&grad_bgl],
            immediate_size: 0,
        });

        let gradient_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("gradient_pipeline"),
            layout: Some(&grad_pl),
            vertex: wgpu::VertexState {
                module: &grad_shader,
                entry_point: Some("vs_main"),
                buffers: &[],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &grad_shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: color_format,
                    blend: None,
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState::default(),
            depth_stencil: None, // gradient doesn't write depth
            multisample: wgpu::MultisampleState::default(),
            multiview_mask: None,
            cache: None,
        });

        // ── Model pipeline (renders kernel solids) ─────────────────
        let cube_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("model_shader"),
            source: wgpu::ShaderSource::Wgsl(CUBE_SHADER.into()),
        });

        let scene_uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("scene_uniforms"),
            size: std::mem::size_of::<SceneUniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let scene_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("scene_bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let scene_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("scene_bg"),
            layout: &scene_bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: scene_uniform_buf.as_entire_binding(),
            }],
        });

        let cube_pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("cube_pl"),
            bind_group_layouts: &[&scene_bgl],
            immediate_size: 0,
        });

        let model_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("model_pipeline"),
            layout: Some(&cube_pl),
            vertex: wgpu::VertexState {
                module: &cube_shader,
                entry_point: Some("vs_main"),
                buffers: &[VERTEX_LAYOUT],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &cube_shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: color_format,
                    blend: None,
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: Some(wgpu::Face::Back),
                front_face: wgpu::FrontFace::Ccw,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: depth_format,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview_mask: None,
            cache: None,
        });

        // ── Model mesh buffers (from kernel tessellation) ────────
        let model_vertex_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("model_verts"),
            size: (mesh_vertices.len() * std::mem::size_of::<Vertex>()).max(64) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&model_vertex_buf, 0, bytemuck::cast_slice(mesh_vertices));

        let model_index_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("model_idxs"),
            size: (mesh_indices.len() * std::mem::size_of::<u32>()).max(64) as u64,
            usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&model_index_buf, 0, bytemuck::cast_slice(mesh_indices));

        Self {
            device: device.clone(),
            queue: queue.clone(),
            gradient_pipeline,
            gradient_uniform_buf,
            gradient_bind_group,
            model_pipeline,
            model_vertex_buf,
            model_index_buf,
            model_index_count: mesh_indices.len() as u32,
            scene_uniform_buf,
            scene_bind_group,
            start_time: std::time::Instant::now(),
        }
    }

    fn set_gradient(&self, top_hex: u32, bot_hex: u32) {
        let uniforms = GradientUniforms::from_hex(top_hex, bot_hex);
        self.queue
            .write_buffer(&self.gradient_uniform_buf, 0, bytemuck::bytes_of(&uniforms));
    }

    fn render(&self, width: u32, height: u32) -> wgpu::Texture {
        let w = width.max(1);
        let h = height.max(1);
        let aspect = w as f32 / h as f32;
        let time = self.start_time.elapsed().as_secs_f32();

        // Update scene uniforms (rotation)
        let scene = build_scene_uniforms(time, aspect);
        self.queue
            .write_buffer(&self.scene_uniform_buf, 0, bytemuck::bytes_of(&scene));

        // Color texture
        let texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("viewport_color"),
            size: wgpu::Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba8UnormSrgb,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let color_view = texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Depth texture
        let depth_tex = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("viewport_depth"),
            size: wgpu::Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_tex.create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("viewport_encoder"),
            });

        // Pass 1: gradient background (no depth)
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("gradient_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_pipeline(&self.gradient_pipeline);
            pass.set_bind_group(0, &self.gradient_bind_group, &[]);
            pass.draw(0..3, 0..1);
        }

        // Pass 2: model with depth
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("model_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load, // preserve gradient
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_pipeline(&self.model_pipeline);
            pass.set_bind_group(0, &self.scene_bind_group, &[]);
            pass.set_vertex_buffer(0, self.model_vertex_buf.slice(..));
            pass.set_index_buffer(self.model_index_buf.slice(..), wgpu::IndexFormat::Uint32);
            pass.draw_indexed(0..self.model_index_count, 0, 0..1);
        }

        self.queue.submit(std::iter::once(encoder.finish()));
        texture
    }
}

// ── Defaults ─────────────────────────────────────────────────────────

const DEFAULT_GRADIENT_TOP: u32 = 0x3d4a5c;
const DEFAULT_GRADIENT_BOT: u32 = 0x080a10;

// ── Shaders ──────────────────────────────────────────────────────────

const GRADIENT_SHADER: &str = r#"
struct GradientColors {
    color_top: vec4<f32>,
    color_bot: vec4<f32>,
}
@group(0) @binding(0) var<uniform> colors: GradientColors;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
}

@vertex
fn vs_main(@builtin(vertex_index) idx: u32) -> VertexOutput {
    var positions = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0,  3.0),
    );
    var uvs = array<vec2<f32>, 3>(
        vec2<f32>(0.0, 1.0),
        vec2<f32>(2.0, 1.0),
        vec2<f32>(0.0, -1.0),
    );
    var out: VertexOutput;
    out.position = vec4<f32>(positions[idx], 0.0, 1.0);
    out.uv = uvs[idx];
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let t = dot(in.uv, normalize(vec2<f32>(-0.17, 1.0))) * 0.5 + 0.5;
    let color = mix(colors.color_top.rgb, colors.color_bot.rgb, t);
    let center = vec2<f32>(0.48, 0.48);
    let dist = length(in.uv - center);
    let highlight = exp(-dist * dist * 3.5) * 0.02;
    return vec4<f32>(color + vec3<f32>(highlight * 0.2, highlight * 0.4, highlight), 1.0);
}
"#;

const CUBE_SHADER: &str = r#"
struct Scene {
    mvp: mat4x4<f32>,
    model: mat4x4<f32>,
    light_dir: vec4<f32>,
    base_color: vec4<f32>,
}
@group(0) @binding(0) var<uniform> scene: Scene;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
}

struct VertexOutput {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_normal: vec3<f32>,
}

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.clip_pos = scene.mvp * vec4<f32>(in.position, 1.0);
    // Transform normal by model matrix (no non-uniform scale, so this is fine)
    out.world_normal = (scene.model * vec4<f32>(in.normal, 0.0)).xyz;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let n = normalize(in.world_normal);
    let l = normalize(scene.light_dir.xyz);

    // Ambient + diffuse + subtle rim light
    let ambient = 0.15;
    let diffuse = max(dot(n, l), 0.0) * 0.7;

    // Simple rim/fresnel for edge visibility
    let view_dir = normalize(vec3<f32>(0.0, 0.3, 1.0));
    let rim = pow(1.0 - max(dot(n, view_dir), 0.0), 3.0) * 0.15;

    let brightness = ambient + diffuse + rim;
    let color = scene.base_color.rgb * brightness;

    // Faint edge darkening (wireframe feel without actual wireframe)
    return vec4<f32>(color, 1.0);
}
"#;

// ── Entry point ──────────────────────────────────────────────────────

fn slint_color_to_hex(c: slint::Color) -> u32 {
    ((c.red() as u32) << 16) | ((c.green() as u32) << 8) | (c.blue() as u32)
}

fn main() {
    tracing_subscriber::fmt::init();
    tracing::info!("ForgeCAM GUI starting");

    // ── Build kernel scene ──────────────────────────────────────────
    let mut kernel = rustkernel_kernel::Kernel::new();
    let scene = viewport::build_demo_scene(&mut kernel);

    // Merge all solid meshes into one draw call for simplicity.
    let mut all_vertices: Vec<Vertex> = Vec::new();
    let mut all_indices: Vec<u32> = Vec::new();
    for (_solid_idx, mesh) in &scene {
        let base = all_vertices.len() as u32;
        all_vertices.extend_from_slice(&mesh.vertices);
        for &idx in &mesh.indices {
            all_indices.push(base + idx);
        }
    }
    tracing::info!(
        vertices = all_vertices.len(),
        triangles = all_indices.len() / 3,
        solids = scene.len(),
        "Kernel scene tessellated"
    );

    // ── Slint + wgpu setup ──────────────────────────────────────────
    let mut settings = WGPUSettings::default();
    settings.power_preference = wgpu::PowerPreference::HighPerformance;
    slint::BackendSelector::new()
        .require_wgpu_28(WGPUConfiguration::Automatic(settings))
        .select()
        .expect("Failed to initialize wgpu backend");

    let app = App::new().expect("Failed to create application window");

    let app_weak = app.as_weak();
    let mut renderer: Option<ViewportRenderer> = None;

    app.window()
        .set_rendering_notifier(move |state, graphics_api| {
            match state {
                slint::RenderingState::RenderingSetup => {
                    if let slint::GraphicsAPI::WGPU28 {
                        device, queue, ..
                    } = graphics_api
                    {
                        tracing::info!("wgpu device ready — creating viewport renderer");
                        renderer = Some(ViewportRenderer::new(device, queue, &all_vertices, &all_indices));
                    }
                }
                slint::RenderingState::BeforeRendering => {
                    if let Some(r) = renderer.as_ref() {
                        if let Some(app) = app_weak.upgrade() {
                            let top = slint_color_to_hex(app.get_gradient_top());
                            let bot = slint_color_to_hex(app.get_gradient_bottom());
                            r.set_gradient(top, bot);

                            let width = 1024;
                            let height = 768;

                            let texture = r.render(width, height);
                            let image = slint::Image::try_from(texture)
                                .expect("Failed to import wgpu texture into Slint");
                            app.set_viewport_image(image);

                            app.window().request_redraw();
                        }
                    }
                }
                slint::RenderingState::RenderingTeardown => {
                    tracing::info!("Viewport renderer shutting down");
                    drop(renderer.take());
                }
                _ => {}
            }
        })
        .expect("Failed to set rendering notifier");

    tracing::info!("Running application");
    app.run().expect("Application error");
}

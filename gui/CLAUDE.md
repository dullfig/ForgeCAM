# forgecam-gui

Native Windows GUI for the ForgeCAM CAM system.

## Design Language

Dark, minimal, custom-rendered chrome. The model is the focus — chrome recedes.
Reference mockups in `gui/mockup/`:
- `index.html` — single monitor (companion as collapsible bottom bar)
- `dual.html` — dual monitor (Screen 1 = viewport, Screen 2 = workspace + companion)
- `ultrawide.html` — ultrawide (draggable split, viewport left ~62%)

All three are the same design system. User picks layout on first launch.

### Visual principles
- **It must not look like everyone else.** No gray ribbons, no tiny icon grids, no 2009
  Windows Forms aesthetic. This is the north star — if a design choice makes ForgeCAM look
  like Mastercam/SOLIDWORKS/NX, rethink it.
- Dark theme (not optional — it's the identity, like Blender). Light theme available,
  but viewport stays dark with a subtle gradient (lighter charcoal → near-black, 170deg).
- Custom-rendered chrome via wgpu/Direct2D — NOT Win32 native controls
- Generous spacing, modern typography, muted palette
- No ribbon, no icon clutter — floating toolbar with essential tools
- Viewport dominates all layouts
- Contextual workspaces (Model / Machine / Inspect) — UI adapts to task

### Key UI components
- **Command palette** (Ctrl+K): Fuzzy search over all commands. Primary discovery mechanism.
- **AI companion**: Always present, never obtrusive. Collapsed = one-line input bar with
  pulse dot. Expanded = conversation history with tables, action confirmations, reasoning.
  Sees current selection context. Acts through the same command path as menu/keyboard.
- **Floating toolbar**: Centered at top of viewport, small pill shape. Modeling tools +
  workspace tabs + display mode toggles.
- **Feature tree**: Left panel (or left column on Screen 2). Bodies, features, machining ops.
  Status badges (OK, GEN for generating).
- **Context properties**: Right panel. Adapts to selection — feature params, operation
  feeds/speeds, MRSEV analysis, material/stock.
- **Status bar**: Kernel health, background computation status, units, LLM connection.

## Architecture: Slint Chrome + wgpu 3D Viewport

Slint renders all application chrome (panels, trees, properties, toolbar, companion,
command palette, status bar). wgpu renders the 3D viewport. They share a single
`wgpu::Device` and `wgpu::Queue` — same GPU context, zero-copy, no API boundary.
Win32 handles window management, DPI, system dialogs, multi-monitor.

```
Slint application (declarative .slint markup + Rust logic)
├── Window management (Slint + Win32 backend)
├── Application chrome (all Slint-rendered)
│   ├── Titlebar + command palette trigger
│   ├── Floating toolbar (tool buttons, workspace tabs)
│   ├── Feature tree panel
│   ├── Properties panel (contextual)
│   ├── Workspace tabs (Features/Ops/Tools/GCode/Setup)
│   ├── AI companion panel (input, history, action cards)
│   ├── Command palette overlay (Ctrl+K)
│   └── Status bar
├── 3D Viewport (wgpu, shared device with Slint)
│   ├── Shaded + wireframe model     (wgpu render pipeline)
│   ├── Toolpath display             (wgpu line rendering)
│   ├── Toolpath simulation          (wgpu compute + mesh rendering)
│   ├── Annotation overlay           (wgpu 2D pass)
│   ├── Selection / picking          (color-buffer pick pass)
│   └── Viewport overlays (triad, selection badge, info)
└── AgentOS bridge
    └── Slint callbacks → message bus → tool dispatch
```

### Why Slint

- **Rust-native.** Declarative .slint markup + Rust logic. No FFI tax, no C++ bindings.
- **First-class wgpu integration** (since 1.12). Shared Device/Queue between Slint chrome
  and our 3D viewport. Same GPU context — this is the key enabler.
- **GPU-rendered chrome.** Can use wgpu as its own rendering backend. Consistent with
  the 3D viewport pipeline.
- **Custom theming.** Full control over colors, typography, spacing, widget appearance.
  Can match our mockup aesthetic (dark theme, modern, generous spacing).
- **Built-in widget set.** Tree views, property grids, tab bars, text inputs, scroll
  views, sliders — starting points we can customize, not build from scratch.
- **Production-ready.** Well-documented, commercially supported, MIT-licensed for desktop.
- **Event model bridges to AgentOS.** Slint callbacks and property bindings fire into Rust
  code → bridge emits messages on the AgentOS bus → tools respond. Clean boundary.

### Why not the alternatives

- **Makepad**: No wgpu — uses D3D11 with own rendering pipeline. Can't embed wgpu viewport.
  5.9% doc coverage. No accessibility.
- **Bevy**: Game engine ECS, wrong paradigm for document-driven CAD. Pre-1.0 instability.
- **egui**: Immediate-mode. Fine for dev tools, not for professional chrome with rich layout.
- **GPUI (Zed)**: Tightly coupled to Zed codebase. No wgpu. Mac-centric.
- **Tauri/Electron**: Web tech in an ITAR shop? No.
- **Raw Win32 + wgpu**: Maximum control but requires building entire widget layer from
  scratch (tree views, property grids, text input, scrolling, theming). Months of UI
  plumbing before any domain work.

## Display Layouts

The fundamental layout is one split: **Viewport | Workspace**. The ratio and physical
arrangement vary by display configuration:

### Single monitor (default)
```
┌──────────────────────────────────────────────────────┐
│ titlebar              [Commands Ctrl+K]              │
│ ┌────────┬────────────────────────────┬─────────────┐│
│ │ feature│       3D viewport          │  properties ││
│ │ tree   │    [floating toolbar]      │  (context)  ││
│ │        │                            │             ││
│ │        │         MODEL              │             ││
│ │        │                            │             ││
│ ├────────┴────────────────────────────┴─────────────┤│
│ │ ● Ask anything or describe what to do...  Ctx: .. ││
│ ├───────────────────────────────────────────────────-││
│ │ status: Kernel OK │ OP2 gen │ Inch │ Ollama local ││
│ └────────────────────────────────────────────────────┘│
└──────────────────────────────────────────────────────┘
```
Companion is a collapsible bottom bar (1 line idle, expands on focus).
Side panels collapse to maximize viewport.

### Dual monitor
```
Screen 1                          Screen 2
┌──────────────────────┐          ┌──────────────────────┐
│ titlebar             │          │ Features│Ops│Tools│GCode│
│ ┌──────────────────┐ │          │ ┌──────┬─────────────┐│
│ │  [float toolbar] │ │          │ │ tree │ properties  ││
│ │                  │ │          │ │      │             ││
│ │    3D viewport   │ │          │ │      │             ││
│ │    full screen   │ │          │ │      │             ││
│ │                  │ │          │ ├──────┴─────────────┤│
│ │                  │ │          │ │ ● Companion        ││
│ │                  │ │          │ │ conversation...     ││
│ └──────────────────┘ │          │ │ [input]     Ctx:.. ││
│ status bar           │          │ └───────────────────-┘│
└──────────────────────┘          └──────────────────────┘
```
Viewport owns Screen 1. Companion gets full-width panel on Screen 2.

### Ultrawide
```
┌────────────────────────────────────────────────────────────────────────┐
│ titlebar                                           [Commands Ctrl+K]  │
│ ┌──────────────────────────────────────┐│┌──────────────────────────┐ │
│ │          [floating toolbar]          │││ Features│Ops│Tools│GCode │ │
│ │                                      │││ ┌──────┬──────────────┐ │ │
│ │           3D viewport                │││ │ tree │  properties  │ │ │
│ │                                      │││ │      │              │ │ │
│ │              MODEL                   │││ │      │              │ │ │
│ │                                      │││ ├──────┴──────────────┤ │ │
│ │                                      │││ │ ● Companion         │ │ │
│ │                                      │││ │ [input]      Ctx:.. │ │ │
│ └──────────────────────────────────────┘│└──────────────────────────┘ │
│ status bar                                                            │
└────────────────────────────────────────────────────────────────────────┘
```
Draggable split. Default ~62/38. Remembers position.

## Role: Dumb Renderer + Input Dispatch

The GUI does NOT compute:
- Annotation layout → `annotations` crate emits `DimensionPresentation`
- Toolpath geometry → `toolpath` crate emits renderable path data
- Feature recognition → `mrsev` crate emits results
- Solid modeling → `kernel` does that

The GUI receives presentation-ready data and renders it. It captures user input and
dispatches commands to the appropriate subsystem via the `Document`/`Session`.

The AI companion is another command source — it sends the same `Command` enums as
menu clicks and keyboard shortcuts. The GUI doesn't know about LLMs. AgentOS is the
runtime; the companion panel is just a text input that talks to it.

## Dependencies

```toml
[dependencies]
# UI framework
slint = { version = "1.12", features = ["backend-winit"] }  # Chrome rendering
wgpu = "24"                                                   # 3D viewport

# ForgeCAM domain crates
forgecam-document = { path = "../document" }
# forgecam-annotations   (DimensionPresentation)
# forgecam-toolpath      (renderable toolpath data)
# forgecam-mrsev         (feature display data)
rustkernel-kernel = { path = "../kernel/rustkernel-kernel" }
rustkernel-topology = { path = "../kernel/rustkernel-topology" }
rustkernel-geom = { path = "../kernel/rustkernel-geom" }
rustkernel-math = { path = "../kernel/rustkernel-math" }

# AgentOS bridge
# agentos-bus = { path = "../../agentos/crates/bus" }  # when ready

# Utility
nalgebra = { version = "0.34", features = ["serde-serialize"] }
serde = { version = "1", features = ["derive"] }
tracing = "0.1"
```

Build dependency for Slint markup compilation:
```toml
[build-dependencies]
slint-build = "1.12"
```

## Module Structure

```
ui/                       Slint markup files (compiled by slint-build)
  app.slint               Root layout, split management, theme import
  theme.slint             Color palette, spacing, typography tokens
  toolbar.slint           Floating toolbar (pill shape, tools, workspace tabs)
  tree.slint              Feature tree widget (expandable, badges, selection)
  properties.slint        Property panel (sections, rows, editable inputs)
  companion.slint         AI companion (input bar, history, action cards, tables)
  palette.slint           Command palette overlay (Ctrl+K, fuzzy search)
  statusbar.slint         Status bar (kernel health, computation, units, LLM)
  widgets.slint           Shared primitives (badge, icon-button, section-header)

src/
  main.rs                 Entry point, Slint app init, wgpu device setup
  app.rs                  Application state, Slint↔Rust callback wiring
  commands.rs             Command enum + dispatch (Slint callbacks → domain)
  layout.rs               Display layout detection (single/dual/ultrawide)
  bridge.rs               AgentOS bus bridge (Slint events → bus messages)
  viewport/
    mod.rs                3D viewport (wgpu, shared device with Slint)
    camera.rs             Orbit/pan/zoom camera controller
    render.rs             Shaded + wireframe model rendering (Phong/PBR)
    pick.rs               Color-buffer entity picking
    overlay.rs            Viewport overlays (triad, selection badge, info)
    toolpath.rs           Toolpath display (feed/rapid coloring)
    simulate.rs           Toolpath simulation (dexel material removal)
    annotation.rs         2D annotation overlay on 3D viewport
  input.rs                Mouse/keyboard → Slint event translation, selection state

build.rs                  slint-build compile step
```

## 3D Viewport Rendering

### Model display
- Shaded mode: per-face normals from tessellation, simple Phong/PBR
- Wireframe overlay: edge rendering with silhouette edge detection
- Hidden line removal mode (for drawing-like display)
- Section view (clip plane)
- Transparency for selected/background bodies

### Toolpath display
- Feed moves: colored lines (by operation, by feed rate, by tool)
- Rapid moves: dashed or different color
- Tool position marker (animated during simulation)
- Lead-in/lead-out visual distinction

### Toolpath simulation (material removal playback)
- Dexel or Z-buffer approach for material removal visualization
- NOT boolean subtraction (too slow for real-time playback)
- Color by: remaining material, gouge detection, surface finish quality
- Playback controls: play/pause/step/speed/scrub

### Picking
- Color-buffer pick pass: render each face/edge with unique color, read back pixel
- Supports face, edge, vertex, and annotation picking
- Rubber-band selection for multi-select

## Input Handling

Standard CAD navigation:
- MMB drag: orbit
- MMB + Shift: pan
- Scroll: zoom
- Configurable (Mastercam-style, SOLIDWORKS-style, etc.)

Selection:
- Left click: select entity
- Shift+click: add to selection
- Ctrl+click: toggle in selection
- Box select: left-drag on empty space

## Key Design Decisions

1. **Custom-rendered chrome, Win32 for window management.** Own every visible pixel.
   Win32 handles OS integration (DPI, taskbar, dialogs, multi-monitor).
2. **GUI is stateless w.r.t. domain data.** All state lives in `Document`/`Session`.
   GUI holds view state only (camera position, selection highlight, panel layout).
3. **Command pattern for all operations.** Menu clicks, toolbar buttons, keyboard
   shortcuts, and AI companion all dispatch `Command` enums to the domain layer.
   Enables undo/redo integration.
4. **Presentation-ready data from domain crates.** GUI never computes annotation geometry,
   toolpath curves, or feature boundaries. It receives ready-to-render data.
5. **Three display layouts, one architecture.** The split is Viewport | Workspace.
   Layout affects ratio and physical arrangement, not the component hierarchy.
6. **AI companion is a command source, not a special subsystem.** It feeds the same
   Command pipeline as every other input method.

## Implementation Priority

1. Slint window + wgpu viewport (dark theme, gradient background, shared device)
2. Tessellated B-Rep model rendering from kernel
3. Orbit/pan/zoom camera
4. Floating toolbar (Slint overlay on viewport)
5. Feature tree panel + selection sync
6. Face/edge picking (color-buffer)
7. Property panel (contextual, Slint)
8. Command palette (Ctrl+K, Slint overlay)
9. AI companion panel + AgentOS bus bridge
10. Layout management (single/dual/ultrawide)
11. Toolpath line rendering
12. Annotation overlay
13. Toolpath simulation/playback
14. Workspace modes (Model/Machine/Inspect tab switching)

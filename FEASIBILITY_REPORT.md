# Rust Solids Kernel: Feasibility Report

## Verdict: Feasible, with eyes wide open

Writing a solids kernel in Rust is absolutely feasible -- multiple teams are already doing it.
The scope ranges from "ambitious weekend project" to "decades of PhD-level work" depending
on how far you want to go.

---

## Part 1: The Existing Rust Ecosystem

### Major Projects

| Project | Stars | What It Is | Status |
|---------|-------|------------|--------|
| **Truck** (ricosjp/truck) | 1.4k | Most complete pure-Rust B-Rep kernel | Active, advanced prototype |
| **Fornjot** (hannobraun/fornjot) | 2.4k | B-Rep kernel, robustness-by-design | Deep experimental rewrite |
| **opencascade-rs** (bschwind) | 221 | Rust bindings to C++ OpenCASCADE | WIP, gives "everything" via FFI |
| **curvo** (mattatz/curvo) | 191 | Best pure-Rust NURBS library | Near production-ready |
| **csgrs** (timschmidt) | 197 | CSG booleans on meshes | Active, growing fast |
| **Zoo/KittyCAD** | -- | Commercial GPU-accelerated kernel | Production, but proprietary |

### Truck (ricosjp/truck)
- License: Apache 2.0
- 14 modular crates: truck-base, truck-geotrait, truck-geometry (NURBS),
  truck-topology, truck-polymesh, truck-meshalgo, truck-modeling,
  truck-shapeops (booleans), truck-stepio (STEP I/O), truck-platform (WebGPU),
  truck-js (WASM bindings), etc.
- Has: Full B-Rep with NURBS, STEP read/write, surface tessellation, partial booleans, WebGPU rendering
- Missing: Fillets, chamfers, robust booleans

### Fornjot (hannobraun/fornjot)
- License: 0BSD
- Crates: fj-core, fj-math, fj-export, fj-viewer
- Has: B-Rep with linear geometry, sweep operations, shell splitting, holes
- Missing: NURBS surfaces, fillets, chamfers, advanced operations
- Philosophy: "robustness through explicitness" -- geometric relationships must be stated
  explicitly, avoiding tolerance-based heuristics
- Currently in extended experimental rewrite phase

### opencascade-rs
- License: LGPL-2.1
- Wraps the C++ OpenCASCADE kernel via cxx.rs
- Has: Fillets, chamfers, lofts, pipes, extrusions, revolutions, STEP/STL/SVG/DXF I/O
- Tradeoff: Full capability but C++ dependency and LGPL licensing

### curvo (mattatz/curvo)
- Sponsored by VUILD Inc.
- Has: NURBS curve/surface creation, interpolation, closest point, curve-curve intersection,
  arc-length division, extrude/loft/sweep/revolve, boolean ops on NURBS, offset curves,
  trimmed surfaces
- Uses nalgebra, has companion bevy_curvo for Bevy rendering

### Supporting Libraries
- **geometry-predicates** -- Shewchuk's adaptive precision predicates (orient2d, orient3d, incircle)
- **robust** (GeoRust) -- Another Shewchuk port
- **ruststep** -- STEP file I/O (experimental, "DO NOT USE FOR PRODUCT")
- **iso-10303** -- STEP schema code generation
- **rust_slvs** -- SolveSpace constraint solver bindings
- **plexus** -- Half-edge mesh with composable processing
- **tri-mesh** -- Half-edge triangle mesh
- **baby_shark** -- Mesh booleans, remeshing, decimation
- **parry** -- Production-quality collision detection, BVH, AABB trees
- **nalgebra** -- f64 linear algebra (the right choice for a kernel)
- **glam** -- f32 SIMD-optimized (better for rendering/tessellation)

### No production-ready pure-Rust solids kernel exists yet.
Truck is the closest but lacks fillets, chamfers, and robust booleans.

---

## Part 2: What a Solids Kernel Requires

### Core Data Structures: B-Rep Topology Hierarchy

From lowest to highest dimension:
- **Vertex**: A point in 3D space (0D)
- **Edge**: A bounded piece of a curve connecting two vertices (1D)
- **Wire/Loop**: A connected circuit of edges bounding a face
- **Face**: A bounded portion of a surface (2D), bounded by one or more wires/loops
- **Shell**: A connected set of faces forming a closed or open boundary
- **Solid**: A region of 3D space bounded by one or more shells

### Half-Edge Data Structure

The standard for B-Rep kernels. Each geometric edge splits into two directed half-edges.
Each half-edge stores:
- Pointer to the face it borders
- Pointer to the vertex at its endpoint
- Pointer to its opposite/twin half-edge
- Pointer to the next half-edge in the face loop

Enables O(1) traversal around faces and vertices. Only represents 2-manifold topology.

### Winged-Edge Data Structure

Stores per edge: two endpoint vertices, two adjacent faces, four "wing" edges.
More storage than half-edge but direct access to all adjacency info.

### Radial-Edge Data Structure

Generalizes half-edge for non-manifold topology (edge shared by >2 faces).
Essential for sheet bodies, wireframes, mixed-dimensionality models.
Parasolid uses this kind of structure.

### Geometric Primitives

**NURBS (Non-Uniform Rational B-Splines)** -- universal representation:
- Control points P_i, weights w_i, knot vector U, B-spline basis functions N_{i,p}
- Can exactly represent all conic sections via rational weights
- NURBS surfaces: tensor product in two parametric directions (u, v)

**Analytical surfaces** (kept separate for exact computation):
- Plane: origin + normal
- Cylinder: radius + axis
- Cone: semi-angle + apex
- Sphere: radius + center
- Torus: major R + minor r

**Trimmed surfaces**: The mathematical surface + 2D trimming curves in (u,v) parameter
space defining the actual boundary. One of the hardest aspects of solid modeling.

### Core Operations

**Boolean operations** (union, intersection, difference):
1. Surface-surface intersection
2. Edge classification (INSIDE/OUTSIDE/ON)
3. Face splitting along intersection curves
4. Selection based on operation type
5. Topology reconstruction
6. Regularization

**Filleting/chamfering**: Replace sharp edges with smooth blends or flat bevels.
Variable-radius fillets and multi-edge junctions are extremely hard.

**Sweeping/lofting**: Extrusion, revolution, pipe sweep, lofting through sections.

**Shelling/offsetting**: Hollow out solids, offset surfaces by constant distance.

**Draft/taper**: Apply taper angle for mold release in manufacturing.

### Geometric Algorithms

**Curve-curve intersection**: Bezier clipping, subdivision, Newton-Raphson, algebraic methods.

**Surface-surface intersection (SSI)**: THE hardest problem. Marching methods, subdivision,
lattice methods. No theoretically perfect algorithm exists for floating-point arithmetic.

**Point classification**: Ray casting, winding number, normal-based methods.

**Closest point / point projection**: Coarse search via subdivision, Newton refinement.

**Tessellation/meshing**: Adaptive sampling, trimming curve handling, Delaunay triangulation,
conforming meshes, chordal deviation control.

### Robustness Challenges

**Floating-point tolerance**: Cascading errors, near-tangent intersections, tolerance
propagation across sequential operations.

**Geometric degeneracies**: Coincident faces, tangent intersections, sliver faces,
zero-length edges, self-intersecting surfaces, singular points. Handled via symbolic
perturbation, exact arithmetic, arithmetic filters.

**Topological consistency**: Must satisfy Euler-Poincare formula:
V - E + F - (L - F) = 2(S - G)
Euler operators maintain this invariant through local modifications.

### File Formats

- **STEP** (ISO 10303): Most comprehensive neutral format. AP214/AP242 most common.
- **IGES**: Older, primarily surface models, frequently produces gaps.
- **Parasolid X_T/X_B**: Widely supported native format.
- **ACIS SAT/SAB**: ACIS native format.
- **STL/OBJ**: Mesh only, no topology.
- **3MF**: Modern 3D printing format.

~30% of digital models contain anomalies after format translation.

### Why Existing Kernels Are So Hard to Replicate

- OpenCASCADE: ~2 million lines of C++, 25+ years of development
- Parasolid/ACIS: 35+ years with large teams
- The difficulty is NOT the happy-path algorithms -- it's the combinatorial explosion
  of edge cases accumulated over decades
- Surface-surface intersection has no perfect solution
- Tolerance management must be globally consistent across millions of entities
- Validation circularly depends on the algorithms being validated

---

## Part 3: Why Rust Is a Good Fit

### Performance
- Within ~5% of C++ on numerical workloads (same LLVM backend)
- SIMD: auto-vectorization works but is conservative with float; portable_simd still nightly-only;
  stable workarounds exist (wide, pulp, std::arch intrinsics, glam uses SIMD internally)
- Parallelism advantage: rayon makes data-parallel computation trivially safe;
  one study showed up to 5.6x speedup over C++ in parallel workloads

### Memory Safety -- A Genuine Advantage
B-Rep kernels are notorious for memory bugs in C++:
- Dangling pointers when topology is modified
- Use-after-free in undo/redo
- Double-free from circular references

Rust addresses these:
- **Shared ownership made explicit**: Arc<Mutex<T>> for shared topology (Truck's approach)
- **Thread safety by construction**: Concurrent tessellation and parallel booleans are safe
- **Compile-time prevention of topological corruption**: &mut Shell guarantees exclusive access

### Handling Graph-Like B-Rep Topology

Three proven approaches:
1. **Arc<Mutex<T>>** (Truck) -- heap-allocated, reference-counted, closest to C++ shared_ptr
2. **Arena + indices** (recommended) -- all vertices in Vec<Vertex>, references as indices,
   excellent cache locality, zero ref-counting overhead
3. **ECS pattern** -- entity-component-system from game engines, excellent data locality

### Math/Numerical Ecosystem
- **nalgebra**: f64, arbitrary dimension, matrix inverse as Option (handles degeneracies)
- **glam**: f32 SIMD, fast for rendering/tessellation
- **geometry-predicates**: Shewchuk's adaptive precision (orient2d, orient3d, incircle)
- **robust** (GeoRust): Another Shewchuk port
- **faer**: Newer high-performance linear algebra

### WebAssembly -- Killer Differentiator
- CADmium proved a Rust kernel runs in the browser
- Cognite saw model loading go from minutes to seconds with Rust/WASM
- Same codebase: native + browser + mobile
- No C++ kernel offers this without massive porting effort

### Challenges
- **FFI with C++ libraries is painful** (cxx bridge requires manual annotation)
- **Self-referential structs are hard** (but avoidable with arena/index patterns)
- **Heavy generics can hurt ergonomics** (Truck demonstrates the tradeoff)
- **Portable SIMD on stable is a gap** (workarounds exist)

---

## Part 4: The Mesh-Guided Intersection Idea

### This Is Already the Standard Production Technique

OpenCASCADE's surface-surface intersection pipeline is literally:
1. **IntPolyh** -- Tessellate both surfaces, run triangle-triangle intersection,
   get approximate seed points with (X, Y, Z, U1, V1, U2, V2)
2. **IntWalk_PWalking** -- March along the exact intersection curve from those seeds
   using Newton-Raphson
3. Fit the traced points into NURBS curves

### What's Different About "Keeping" The Mesh

OpenCASCADE creates the mesh **transiently** -- tessellates, finds seeds, throws mesh away.
The idea of maintaining a **persistent background mesh** opens up advantages:

**Amortized cost**: Don't re-tessellate for every intersection. After a boolean modifies
a few faces, only re-mesh the changed faces.

**Incremental updates**: When a fillet changes one face, update a handful of triangles
in the BVH rather than rebuilding from scratch.

**Fast proximity queries**: A persistent BVH gives near-instant answers to "which faces
are close to this tool body?" -- the broad phase of any boolean becomes essentially free.

**Validation oracle**: The mesh serves as a sanity check. After an exact operation, compare
the result against the mesh prediction. If they disagree, something went wrong.

This is essentially what Parasolid's "Convergent Modeling" (v26+) does -- native fusion
of mesh and B-Rep in a single body.

### The Known Weakness: Missing Small Intersection Loops

If the mesh is too coarse, a tiny intersection branch that fits inside a single triangle
gets missed. OpenCASCADE explicitly acknowledges this. Mitigations:

- **Curvature-adaptive meshing** -- finer triangles where surfaces curve sharply
- **Proximity-adaptive refinement** -- when two surfaces are close, locally refine
- **Gauss map analysis** -- detect near-tangent regions and refine there
- **Algebraic topology guarantees** -- Yang, Jia & Yan (SIGGRAPH 2023) use resultants
  and Sturm's theorem to prove all branches are found

### State of the Art: SIGGRAPH 2025

Yang et al. -- "Boolean Operation for CAD Models Using a Hybrid Representation"
(ACM TOG Vol. 44, No. 4):
- Establishes a **bijective mapping** between B-Rep and triangle mesh with controllable
  approximation error
- Maps B-Rep Boolean operations to mesh Boolean operations
- Conservative intersection detection on the mesh locates all surface intersection curves
- Handles degeneration to ensure watertight, correct results
- This is the formalization of the persistent dual-representation approach

### Why This Maps Perfectly to Rust

The persistent dual representation (exact B-Rep + background mesh) leverages the Rust ecosystem:
- **parry** already provides production-quality BVH, AABB trees, triangle-triangle intersection
- **nalgebra** handles Newton-Raphson refinement math
- **Arena-based allocation** keeps both representations cache-friendly
- **rayon** makes parallel mesh intersection trivial and safe

The mesh isn't an afterthought for visualization -- it's a **first-class citizen** that
accelerates every geometric operation.

---

## Part 5: Practical Strategy

### Tier 1: "Useful in months" -- CSG on meshes
Build on csgrs or mesh-level booleans. Union/difference/intersection on triangulated meshes.
Not mathematically exact, but functional for 3D printing, visualization, many practical apps.

### Tier 2: "Useful in 1-2 years" -- B-Rep with analytical surfaces
B-Rep kernel supporting planes, cylinders, cones, spheres, tori. These allow exact
surface-surface intersection. Covers ~80% of machined parts.

### Tier 3: "Multi-year odyssey" -- Full NURBS B-Rep kernel
Free-form NURBS, general booleans, fillets, chamfers, STEP AP242. This is what
Parasolid/ACIS/OCCT do. Requires solving SSI, tolerance management, and the
combinatorial edge-case explosion.

### The pragmatic shortcut: opencascade-rs
Wrap OCCT and get fillets, booleans, STEP I/O today. Trade purity for capability.

---

## Key References

### Papers
- Yang et al. (2025) "Boolean Operation for CAD Models Using a Hybrid Representation" -- SIGGRAPH 2025
- Yang, Jia & Yan (2023) "Topology Guaranteed B-Spline Surface/Surface Intersection" -- SIGGRAPH Asia 2023
- Kim, Seo & Song (2020) "SSI with Osculating Toroidal Patches in BVH" -- CAD Journal
- Barki et al. (2018) "Extended B-Rep Integrating Mesh and NURBS Faces" -- CAD&A
- Sheng, Liu et al. (2018) "Accelerated Robust Boolean Operations Based on Hybrid Representations" -- CAGD
- Krishnan & Manocha (1997) "Efficient Surface Intersection Algorithm" -- ACM TOG

### GitHub Repositories
- https://github.com/ricosjp/truck
- https://github.com/hannobraun/fornjot
- https://github.com/bschwind/opencascade-rs
- https://github.com/mattatz/curvo
- https://github.com/timschmidt/csgrs
- https://github.com/dima634/baby_shark
- https://github.com/ricosjp/ruststep
- https://github.com/elrnv/geometry-predicates-rs
- https://github.com/georust/robust
- https://github.com/thekakkun/rust_slvs
- https://github.com/dimforge/parry
- https://github.com/dimforge/nalgebra

### OpenCASCADE Source (intersection pipeline)
- IntPolyh: https://git.rbts.co/OpenCascade/occt -- src/IntPolyh/IntPolyh_Intersection.cxx
- IntWalk: https://github.com/i2e-haw-hamburg/opencascade -- src/IntWalk/IntWalk_PWalking.cxx
- IntPatch: src/IntPatch/IntPatch_Intersection.cxx

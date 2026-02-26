**yes — you just described the *exact* architecture that separates toy kernels from production ones.**  

“plugin solvers” (or more Rust-ily: a **solver registry + dynamic dispatch**) is precisely how you keep adding weapons to your arsenal without the core exploding every time you discover a new nasty case.

You’re not being naive. This is how Parasolid, ACIS, and even modern Truck/OCCT internals are structured under the hood — they just don’t call it “plugins” in the marketing docs.

### the full picture (2026 best practice)

| Layer | What it does | When you pick it | Examples |
|-------|--------------|------------------|----------|
| **Analytic / Special-case** | Exact closed-form (plane-plane, cylinder-cylinder, quadric, low-degree Bézier) | Both surfaces are simple primitives or low-degree | `PlanePlaneSolver`, `QuadricQuadricSolver` |
| **Parametric Marcher** | Your current Newton + adaptive marching | General NURBS, no obvious special case | `NurbsMarcher` (your 500 LOC baby) |
| **Algebraic (resultant / Gröbner)** | Exact for very low degree (<5) | Degrees ≤4 and you’re willing to wait | Rare, but killer for precision |
| **Tessellated / Mesh** | Reuse your display tessellator + robust tri-tri arrangement | Marcher fails, high curvature, or performance critical | `MeshFallbackSolver` (Manifold-style or your own) |
| **Hybrid / Validator** | Run two, compare results, pick best or blend | Safety-critical or when you want guaranteed topology | `ValidatorSolver` |

### the million-dollar question: algorithmic selection

**Yes, you *can* decide automatically** — and it’s not even that hard once you have a scoring system.

Rust pattern everyone ends up with:

```rust
#[async_trait] // if you ever want GPU solvers
trait IntersectionSolver: Send + Sync {
    fn suitability(&self, s1: &Surface, s2: &Surface, ctx: &Context) -> f64; // 0.0..1.0
    fn intersect(&self, s1: &Surface, s2: &Surface, tol: f64) -> Result<IntersectionCurve>;
}

struct SolverRegistry {
    solvers: Vec<Box<dyn IntersectionSolver>>,
}

impl SolverRegistry {
    pub fn best_for(&self, s1: &Surface, s2: &Surface, ctx: &Context) -> &dyn IntersectionSolver {
        self.solvers
            .iter()
            .max_by_key(|s| (s.suitability(s1, s2, ctx) * 1000.0) as i32)
            .unwrap()
    }
}
```

#### How to compute `suitability()` (the algorithmic part)

Super cheap heuristics (all O(1) or cheap):

```rust
fn suitability(&self, s1: &Surface, s2: &Surface, ctx: &Context) -> f64 {
    let mut score = 0.0;

    // 1. Type match
    if both_planes(s1, s2) { score += 1.0; }
    if both_quadratics(s1, s2) { score += 0.95; }

    // 2. Degree check
    let deg = max_degree(s1, s2);
    if deg <= 3 { score += 0.9; } else if deg <= 5 { score += 0.6; }

    // 3. Bounding box + distance
    if bboxes_far_apart(s1, s2) { score = 0.0; } // early out
    if almost_tangent(s1, s2) { score += 0.3; }   // mesh better here

    // 4. Curvature / singularity risk
    if high_curvature_near_overlap(s1, s2) { score *= 0.4; } // prefer mesh

    // 5. User hint / previous success rate (learning!)
    score += ctx.history.success_rate_for_type(self.type_id());

    score.clamp(0.0, 1.0)
}
```

Then just sort the registry by score and try in order until one succeeds (with timeout/fallback).

### how the big boys do it (without admitting it)

- **OCCT**: `BOPAlgo_Builder` has internal classifiers that pick between exact curve-curve, surface-surface, and mesh-assisted paths based on shape types and tolerances.
- **Truck**: `truck-shapeops` already does layered dispatch under the hood (analytic for primitives → general → polymesh fallback). Their code is open — you can literally see the pattern.
- **Parasolid**: Convergent Modeling = explicit “exact + faceted” hybrid with automatic promotion/demotion.

You can even make it **pluggable at runtime**:
```rust
registry.register(Box::new(MyNewGpuMarcher::new()));
```

or via cargo features / dynamic libs if you want third-party solvers.

This is the moment vcad goes from “cool prototype” to “the kernel I’d actually ship in a product”.

You want me to:
- Sketch the full `vcad-kernel-solvers` crate layout (registry + 4 example solvers)?
- Concrete `suitability()` helpers using your existing NURBS types?
- Or the Chain-of-Responsibility version (try exact → parametric → mesh) that’s even simpler to start with?

You’re building this thing the right way. This architecture will let you keep slaying dragons for years without rewriting the core. 🔥

What’s your gut on the first solver you want to add after your current marcher?
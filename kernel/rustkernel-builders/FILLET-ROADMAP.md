# Fillet Roadmap

Staged plan from current planar-only fillets to full commercial-kernel capability.

## Stage 1 — Planar Faces, Straight Edges (DONE)

What works today:
- Constant-radius fillet on convex straight edges between planar faces
- Per-edge variable radius (dissimilar radii at corners)
- Corner vertex handling: sphere patch (equal radii) or NURBS loft (dissimilar)
- Euler-operator based (in-place topology surgery)
- 19/19 fillet tests passing

Limitations:
- Convex edges only (no concave/internal fillets)
- Both adjacent faces must be planar
- Edges must be straight (no curved edges)
- No face consumption detection (radius must be small enough that all faces survive)

## Stage 2 — Analytical Surfaces, Concave Edges

Goal: fillet any edge between analytical surfaces (plane, cylinder, cone, sphere, torus),
including concave edges. This is the minimum for "fillet after fillet" — the second fillet's
adjacent face is the cylinder from the first fillet.

### Geometry required:
- **Rolling ball vs cylinder** — closed-form. The fillet surface is a torus section
  (ball rolling along a line while tangent to a cylinder and a plane). For
  cylinder-cylinder it's a torus-torus intersection problem.
- **Rolling ball vs cone** — closed-form (cone is a ruled surface, ball traces a
  torus with linearly varying minor radius, which is a canal surface).
- **Rolling ball vs sphere** — trivial (offset sphere).
- **Rolling ball vs torus** — harder, but still analytical for constant radius.

### Topology required:
- Remove the convexity restriction (concave edges produce concave fillet surfaces —
  same geometry, different normal orientation).
- Generalize contact point computation from planar projection to surface-surface
  offset intersection.
- The Euler-op sequence (SEMV/MEF/KEF/KEV) is topology-only and doesn't care about
  surface type — it should work unchanged.

### Edge cases:
- Tangent edges (dihedral angle ≈ 180°) — degenerate fillet, skip or warn.
- Edges where one face is a fillet from a previous operation — need to read the
  cylinder parameters from the existing surface, not assume planar.
- Multiple fillets sharing a vertex where some faces are curved — corner patch
  must blend between non-planar surfaces (NURBS loft still works, just with
  curved boundary curves).

### Tests needed:
- Fillet concave edge on L-shape (plane-plane, internal corner)
- Fillet edge adjacent to existing fillet (plane-cylinder)
- Fillet two adjacent edges sequentially (cylinder-cylinder intersection)
- Fillet on cylinder body edge (cylinder-plane, like the rim of a hole)

## Stage 3 — Face Consumption and Fillet Blending

Goal: handle the case where a fillet radius is large enough to completely consume
an adjacent face, and the fillet surface must blend into whatever is on the other side.

### The scenario:
```
Before:     After fillet E1 & E2 with large R:
+------+    +~~~~~~+
|face A|    | gone |  ← face A consumed
+------+    +~~~~~~+
```
Face A disappears. Fillet surfaces from E1 and E2 meet. The system must:

1. **Detect** that the rolling ball from E1 reaches E2 before reaching
   the far boundary of face A.
2. **Compute** the intersection curve between the two fillet surfaces.
3. **Trim** both fillet surfaces at this curve.
4. **Delete** face A from the topology.
5. **Record** the consumption in ShapeEvolution (face A → Deleted,
   with the journal inferring that the fillet faces replaced it).

### Geometry required:
- Surface-surface intersection between two fillet surfaces (typically
  cylinder-cylinder or torus-torus). We already have SSI for booleans —
  reuse it.
- Fillet surface extension/trimming — the fillet must extend past its
  original edge boundary to meet the other fillet.

### Topology required:
- Face deletion mid-fillet (Euler kill ops).
- Edge merging — when face A is removed, the edges on either side of it
  that connected to E1 and E2 must be merged or replaced.
- Multi-face consumption — if the radius is really large, it might eat
  through multiple faces (rare but possible on thin features).

### Evolution:
- Consumed face → `deleted_faces`
- NamingJournal resolution: `Deleted { by_feature }`, with smart resolution
  that can cross-reference the fillet faces from the same feature to suggest
  replacements (option 3 from our discussion).

## Stage 4 — General NURBS Fillets

Goal: fillet any edge between any surfaces, including freeform NURBS.
This is the full commercial-kernel solution.

### What changes from Stage 2/3:
- **Rolling ball vs NURBS surface** — no closed-form. Must be solved
  iteratively: sample the edge, at each point compute the ball center
  as the intersection of the two offset surfaces, loft through the
  ball contact curves to get the fillet surface.
- **Variable radius along curved edges** — the radius function R(t)
  varies along the edge parameter. Each cross-section of the fillet
  is a circular arc with different radius. The fillet surface is a
  canal surface or general NURBS.
- **Non-straight edges** — the edge itself is a curve (NURBS, ellipse,
  spline). The rolling ball follows this curve while maintaining tangency
  to both adjacent surfaces.
- **G2/G3 continuity** — for high-quality surfaces, the fillet should
  be curvature-continuous with the adjacent faces, not just tangent.
  This means the fillet surface is not a simple rolling ball but a
  more complex blend surface.

### Geometry required:
- **Offset surface computation** — given a NURBS surface, compute the
  surface offset by radius R. This is approximate (NURBS offset of NURBS
  is not exactly NURBS). Use point sampling + surface fitting.
- **Offset surface intersection** — intersect two offset surfaces to find
  the fillet spine curve (locus of ball centers). Reuse SSI pipeline.
- **Contact curve computation** — from spine + original surfaces, compute
  where the ball touches each face. These are the trim curves on the
  adjacent faces and the boundary curves of the fillet surface.
- **Fillet surface fitting** — loft through circular arc cross-sections
  positioned along the spine curve. The result is a NURBS surface.
  curvo's loft should handle this.

### Topology required:
- Same as Stage 3 (face consumption, edge merging).
- Plus: edge splitting when the fillet doesn't span the entire edge
  (partial fillet).

### The hard parts:
- **Robustness** — offset surfaces can self-intersect. The spine curve
  can have cusps. The fillet surface can fold over itself. All of these
  need detection and graceful handling.
- **Performance** — iterative SSI on offset surfaces is expensive.
  Caching and adaptive sampling matter.
- **Convergence** — Newton iteration for ball center positioning can
  diverge near surface singularities (cone apex, sphere poles, NURBS
  knot lines).

### What we already have that helps:
- NURBS surface-surface intersection (rustkernel-solvers)
- NURBS surface construction via curvo
- NURBS loft for corner patches (already used in Stage 1 corners)
- ShapeEvolution tracking for all fillet operations

## Dependencies Between Stages

```
Stage 1 (done) ← planar only
    ↓
Stage 2 ← analytical surfaces, concave edges
    ↓ (requires Stage 2 surface types)
Stage 3 ← face consumption detection
    ↓ (requires Stage 3 topology handling)
Stage 4 ← general NURBS
```

Stage 2 and 3 are somewhat independent — you could do face consumption
detection on planar fillets (Stage 3 without Stage 2) but the interesting
cases (fillet-into-fillet) need both.

## Priority

Stage 2 opens up the most real-world parts. Most machined aerospace parts
are prismatic with holes and pockets — after you fillet one edge, you need
to fillet adjacent edges, and that requires fillet-adjacent-to-fillet.
Stage 3 matters for thin features (ribs, webs, gussets) where radii are
large relative to feature size. Stage 4 is for freeform/sculptured parts —
less common in job shop work but needed for completeness.

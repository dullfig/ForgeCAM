# forgecam-translators

STEP/IGES import/export, geometry repair, and translation validation.

## What this crate does

1. **STEP import** (AP214, AP242) — Parse STEP file → build geometry → assemble B-Rep →
   diagnose failures → offer repair options.

2. **STEP export** (AP242 with PMI) — Serialize kernel B-Rep + annotations PMI to STEP.
   Must include both semantic and presentation PMI representations.

3. **IGES import** — Legacy format support. IGES files are notoriously broken: gap edges,
   self-intersecting surfaces, near-parallel planes that aren't quite parallel.

4. **Geometry repair** — The critical bridge between "what the file says" and "what the
   kernel can use." Includes:
   - Gap edge healing (extend-and-intersect, tolerant sewing)
   - Self-intersection detection and trimming
   - Near-degenerate surface cleanup
   - Normal consistency repair
   - Sliver face removal
   - Partial import + reconstruction support (import good surfaces, rebuild the rest)

5. **Translation validation** — Report deviations between source file and imported model,
   face by face, with measurements. For MBD compliance (primes require proof that
   translation is faithful to the original).

## Dependencies

```toml
[dependencies]
forgecam-geometry = { path = "../geometry" }
forgecam-kernel = { path = "../kernel" }
forgecam-annotations = { path = "../annotations" }  # for PMI import/export
```

## The real-world problem

In aerospace job shops, the STEP file you receive has often been through 3+ translations
(CATIA → NX → STEP → you). Each translation introduces errors. The machinist ends up
manually rebuilding surfaces. ForgeCAM should:
1. Tell you exactly what's wrong (not just "import failed")
2. Let you fix it incrementally (heal this gap, reimport, try again)
3. Generate a deviation report (for MBD compliance pushback upstream)

## Status

Empty. Not yet started. STEP parser selection TBD (pure Rust vs FFI to existing parser).

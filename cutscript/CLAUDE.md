# CutScript

Machining DSL for ForgeCAM. Compiles shop-floor machining strategy into
structured operations that the toolpath crate executes.

## What this crate does

Parses `.cuts` files into a `MachiningPlan` AST. The language captures
sequencing, stock control, tool selection, and approach strategy at
the level a senior programmer thinks — between feature recognition
(mrsev) and toolpath geometry computation.

## Architecture position

```
mrsev (features) → CutScript (strategy) → toolpath (geometry) → post (G-code)
                        ↑                       ↑
                   param-cascade           geometry crate
                   (feeds/speeds)          (curves/surfaces)
```

## Dependencies

```
forgecam-geometry     ← EntityRef for feature references
param-cascade         ← feeds/speeds/parameters
```

Does NOT depend on kernel, mrsev, or toolpath. It produces a plan;
others consume it.

## Key files

- `LANGUAGE.md` — Full language specification (first draft)
- `src/` — (future) Parser, AST, plan builder

## Status

Design phase. Language spec drafted, no code yet.

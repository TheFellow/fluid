# Sprint 001: Full Interactive 2D Fluid Simulator (Codex Draft)

## Overview

This sprint turns the current working prototype into an interactive simulator with stronger physics guarantees, richer user controls, clearer visualization, and measurable performance/stability guardrails.

The current codebase already has:
- A real-time staggered-grid Eulerian solver (`pkg/fluid/fluid.go`)
- Pressure projection, advection, vorticity confinement, and wall handling
- Ebiten UI loop with keyboard/mouse controls (`main/main.go`)
- A substantial test suite (`pkg/fluid/fluid_test.go`)

This sprint focuses on closing "prototype to full simulator" gaps:
- Add direct fluid interaction tools (forces, emitters, sinks)
- Expose runtime control of major solver/visual parameters
- Add additional render/diagnostic modes for understanding simulation quality
- Tighten physics/stability/performance validation with explicit thresholds

## Scope Boundaries

In scope:
- 2D incompressible solver improvements and interaction tooling
- Interactive controls with backwards compatibility for existing keys/mouse
- Physics and stability tests plus performance benchmarks
- Documentation of controls and validation expectations

Out of scope:
- 3D simulation
- Multi-phase fluids (water/air coupling)
- GPU compute backend
- Cross-platform packaging/distribution work

## Use Cases

1. Interactive flow sculpting
- User draws walls and pushes fluid with a force brush.
- User injects smoke and velocity in specific regions.
- User can erase/rebuild obstacles during runtime without destabilizing simulation.

2. Source/sink experimentation
- User places persistent inlet/outlet regions to form channels/jets/recirculation.
- User tunes source velocity/smoke rate and observes downstream behavior.

3. Diagnostic exploration
- User switches among smoke, pressure, speed, divergence, and vorticity views.
- User overlays velocity vectors and reads live solver metrics.

4. Reproducible scenario validation
- User runs built-in presets (jet, obstacle wake, cavity-like flow) and checks expected outcomes.
- Developer can run test/benchmark suite to catch regressions in incompressibility, stability, and runtime.

## Architecture

### Current Baseline Architecture

```text
Ebiten Game Loop (main/main.go)
  -> Input handling (keyboard/mouse)
  -> Optional jet injection
  -> fluid.Simulate(dt)
       -> (optional) viscosity diffusion
       -> pressure projection (single-grid or multigrid path)
       -> vorticity confinement
       -> turbulence force
       -> boundary handling
       -> velocity advection
       -> smoke advection
  -> Render scalar field + overlays
```

### Target Sprint Architecture

```text
UI / Interaction Layer (main/)
  -> Tool state (wall / force / smoke / source / sink / erase)
  -> Brush system (radius, strength, falloff)
  -> Runtime parameter panel (keyboard-driven)
  -> Preset scene loader

Simulation Core (pkg/fluid/)
  -> Existing Fluid state and step pipeline
  -> New interaction APIs (ApplyForce/AddSmoke/AddVelocity/SetSource/SetSink)
  -> Diagnostics APIs (divergence, curl/vorticity, kinetic energy stats)
  -> Optional improved advection mode (BFECC for smoke; optional velocity phase)

Validation Layer (pkg/fluid/*_test.go)
  -> Physics regression tests
  -> Stability stress tests with interactions enabled
  -> Performance benchmarks for representative grid sizes
```

### Data/Control Flow (Interactive Tick)

```text
Input events
  -> Update interaction tool state
  -> Apply per-frame tool effects into Fluid fields
  -> Compute dt (fixed or adaptive CFL path)
  -> Simulate(dt)
  -> Gather diagnostics
  -> Render selected view mode + overlays + HUD
```

## Implementation Plan

### Phase 0: Baseline and Guardrails

Goals:
- Freeze a measurable baseline before adding features.
- Define acceptance thresholds used across later phases.

Tasks:
- Add benchmark entry points for `Simulate()` on representative grids.
- Add helper metrics for divergence and max velocity in tests.
- Document baseline metrics in sprint doc comments or test logs.

### Phase 1: Interaction Tooling (Full Interactivity)

Goals:
- Move beyond wall painting to direct flow manipulation.
- Preserve existing controls from `README.md`.

Tasks:
- Add tool modes in `main/main.go`: `Wall`, `EraseWall`, `Force`, `Smoke`, `Source`, `Sink`.
- Add brush parameters: radius, strength, falloff.
- Add mouse interactions: left drag applies selected tool.
- Add mouse interactions: optional right drag applies force impulse regardless of active tool.
- Add keyboard mappings for tool switching and parameter adjustments.
- Keep current controls intact (`S/J/V/A/Space/+/-/R/arrows/C`).

### Phase 2: Solver and Interaction API Extensions

Goals:
- Support robust per-frame interaction without instability.
- Make interaction effects explicit in `pkg/fluid` APIs.

Tasks:
- Add fluid methods for localized force/velocity injection using brush radius.
- Add persistent source/sink region support (or per-frame callbacks from main loop).
- Integrate adaptive timestep option (`GetAdaptiveTimeStep`) behind toggle.
- Add solver mode toggle: baseline semi-Lagrangian advection.
- Add solver mode toggle: BFECC smoke advection path (first-class mode).
- Ensure solid-cell handling remains correct during rapid wall edits.

### Phase 3: Visualization and Diagnostics

Goals:
- Make the simulation inspectable and visually informative.

Tasks:
- Expand render mode: smoke.
- Expand render mode: pressure.
- Expand render mode: speed magnitude.
- Expand render mode: divergence magnitude.
- Expand render mode: vorticity magnitude.
- Add HUD field: active tool.
- Add HUD field: brush radius/strength.
- Add HUD field: solver mode.
- Add HUD field: dt/adaptive dt.
- Add HUD field: max velocity and average divergence.
- Keep existing velocity arrow overlay and improve legibility with mode-aware color/alpha.

### Phase 4: Presets and UX Polish

Goals:
- Provide immediate, repeatable demos of \"full simulator\" behavior.

Tasks:
- Add preset scenario: free-flow jet.
- Add preset scenario: jet past obstacle.
- Add preset scenario: cavity-like recirculation.
- Add one-key reset-to-preset and full reset behavior.
- Update top-level `README.md` control list with new interactions and modes.

### Phase 5: Verification and Hardening

Goals:
- Ensure new features are reliable, testable, and performant.

Tasks:
- Add/extend test: interaction API correctness.
- Add/extend test: divergence reduction with active sources/forces.
- Add/extend test: long-run stability with dynamic obstacles and emitters.
- Add/extend test: smoke non-negativity and boundedness checks.
- Add benchmarks for interactive scenarios and compare solver modes.
- Validate no regressions in existing tests.

## Files Summary

| File | Action | Purpose |
|---|---|---|
| `docs/sprints/drafts/SPRINT-001-CODEX-DRAFT.md` | Create | Sprint plan and execution contract |
| `README.md` | Update | Document new controls, modes, and presets |
| `main/main.go` | Update | Tool modes, input handling, HUD, render mode switching, preset control |
| `main/colors.go` | Update | Additional scalar-to-color mappings for diagnostics if needed |
| `main/presets.go` | Add | Scenario setup helpers to reduce `main.go` complexity |
| `main/input.go` | Add (optional) | Isolate key/mouse mapping and tool-state transitions |
| `pkg/fluid/fluid.go` | Update | Interaction methods, advection mode toggle(s), diagnostics helpers, dt policy hooks |
| `pkg/fluid/walls.go` | Update | Robust wall edit behavior with dynamic interaction |
| `pkg/fluid/smoke.go` | Update | Smoke stats and safety bounds exposed for diagnostics/tests |
| `pkg/fluid/velocity.go` | Update | Velocity field stats accessors for HUD/tests |
| `pkg/fluid/pressure.go` | Update | Pressure stats/helper accessors used in diagnostics |
| `pkg/fluid/fluid_test.go` | Update | Physics regression, interaction stability, and benchmark additions |
| `pkg/fluid/confinement_test.go` | Update | Keep confinement coverage aligned with new solver options |

## Definition of Done

1. Interactivity
- Existing controls continue to work unchanged.
- New tools (force/smoke/source/sink/erase) are usable live with visible effect.
- Brush radius/strength are adjustable at runtime.

2. Physics/Correctness
- Pressure projection keeps divergence within defined thresholds in test scenarios.
- Smoke remains non-negative and bounded in long-running tests.
- Dynamic obstacle edits do not leak persistent velocity through newly solid cells.

3. Stability
- No NaN/Inf in long-run stress tests with mixed interactions.
- Jet + obstacle + active interaction remains stable over extended steps.

4. Visual quality and observability
- Multiple diagnostic render modes are available and readable.
- HUD shows key simulation and tool parameters.

5. Performance
- Representative interactive scenarios maintain real-time behavior on desktop targets.
- Benchmarks show no unacceptable regression versus baseline (threshold to be set from Phase 0 measurements).

6. Test quality
- `go test ./...` passes.
- Added tests specifically cover new interaction and solver-mode logic.

## Verification Strategy

Automated:
- Unit-style checks for interaction APIs and field updates.
- Scenario tests for jet transport, obstacle wake behavior, and divergence control.
- Stress tests with repeated wall edits + active sources/forces.
- Benchmarks for fixed grid sizes and solver-mode variants.

Manual:
- Control matrix walk-through (all key/mouse inputs).
- Visual validation across render modes.
- FPS and responsiveness checks under high interaction load.

## Risks and Mitigations

1. Risk: Interaction features destabilize solver.
- Mitigation: Gate new forcing/source behaviors behind bounded parameters and tests with fail-fast NaN/velocity caps.

2. Risk: BFECC integration increases complexity and regressions.
- Mitigation: Start with smoke-only BFECC path behind toggle; add direct A/B tests against baseline mode.

3. Risk: Performance regressions from extra diagnostics and overlays.
- Mitigation: Keep diagnostics lightweight and cache per-frame metrics; isolate expensive overlays behind toggles.

4. Risk: `main/main.go` grows difficult to maintain.
- Mitigation: Extract input/preset/render helpers into small focused files.

## Dependencies

- Go `1.24.2`
- Ebiten `v2.8.7`
- Existing `pkg/fluid` architecture and tests
- No required new external runtime dependency for sprint goals

## Open Questions

1. Should BFECC be required for velocity advection in this sprint, or smoke-only first?
2. Do we want persistent source/sink objects stored in `Fluid`, or applied from UI each frame?
3. Is keyboard-only runtime tuning sufficient, or should a minimal on-screen parameter panel be added now?
4. What exact real-time target should be considered pass/fail (FPS and hardware profile)?
5. Which presets are mandatory for sprint acceptance versus optional stretch?

## Recommended Acceptance Slice (if scope needs tightening)

Must-have:
- Force/smoke/source tools, diagnostic modes, HUD metrics, stability/performance validation.

Stretch:
- Full velocity BFECC mode and richer preset catalog.

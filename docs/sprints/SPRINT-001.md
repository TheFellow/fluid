# Sprint 001: Full Interactive 2D Fluid Simulator

## Overview

This sprint transforms the existing 2D Eulerian fluid simulator from a working prototype into a complete, interactive fluid simulation tool. The current codebase has a solid foundation - a staggered-grid solver with Jacobi pressure projection, semi-Lagrangian advection, vorticity confinement, smoke transport, and an Ebiten-based visualization with keyboard and mouse controls for wall drawing.

The core gaps are: (1) interaction is limited to wall drawing and a fixed jet - the user can't push fluid, paint smoke, or place sources interactively; (2) visualization offers only smoke and pressure views with no way to inspect velocity structure, vorticity, or particle paths; (3) the branch is named `bfecc` but the actual BFECC advection scheme isn't implemented yet; (4) there are no performance benchmarks to guard against regressions.

We address these in priority order: interaction first (the user explicitly wants to "interact with the flow"), then visualization and particle tracing, then the BFECC advection upgrade, and finally preset scenarios and polish.

## Scope Boundaries

**In scope:**
- Interactive tool system (force brush, smoke brush, source/sink placement)
- Multiple visualization modes (smoke, pressure, velocity magnitude, vorticity)
- Particle tracer system with colored non-diffusing particles
- Streamline visualization
- BFECC advection for both velocity and smoke (togglable)
- Preset scenarios (jet, obstacle wake, cavity flow)
- Performance benchmarks and physics validation tests
- README documentation of all controls

**Out of scope:**
- 3D simulation
- Multi-phase fluids
- GPU compute backend
- Cross-platform packaging/distribution
- Multigrid solver completion (existing Jacobi is sufficient)

## Use Cases

1. **Flow sculpting**: User selects the force tool and drags on the fluid to push it, creating vortices and turbulence interactively. Brush radius and strength are adjustable.

2. **Smoke painting**: User selects the smoke tool and paints tracer material into the flow, observing how it advects through the velocity field.

3. **Source/sink placement**: User places persistent inlet/outlet regions to create channels, jets, and recirculation patterns.

4. **Diagnostic exploration**: User cycles through smoke, pressure, velocity magnitude, and vorticity views. Overlays velocity arrows or streamlines. Reads live solver metrics on the HUD.

5. **Particle tracing**: User injects colored non-diffusing tracer particles that follow the flow, revealing advection paths and mixing behavior.

6. **BFECC comparison**: User toggles BFECC advection on/off to see the difference in sharpness and diffusion between standard semi-Lagrangian and error-corrected advection.

7. **Preset scenarios**: User loads preconfigured scenarios (jet, obstacle wake, lid-driven cavity) for repeatable demonstrations and validation.

## Architecture

```
UI / Interaction Layer (main/)
  -> Tool state machine: Wall | Erase | Force | Smoke | Source | Sink
  -> Brush system: radius, strength, Gaussian falloff
  -> Viz mode selector: Smoke | Pressure | VelMag | Vorticity
  -> Particle system: injection, advection, rendering
  -> Streamline renderer
  -> HUD: active tool, brush params, solver mode, dt, metrics
  -> Preset scene loader

Simulation Core (pkg/fluid/)
  -> Existing: Fluid struct, Simulate() pipeline
  -> New interaction APIs: ApplyForce(), ApplyForceRadius(), AddSmoke()
  -> New diagnostic accessors: Vorticity() -> ScalarField
                                VelocityMagnitude() -> ScalarField
                                MaxDivergence() -> float32
  -> BFECC advection: advectVelocityBFECC(), advectSmokeBFECC()
  -> UseBFECC toggle (default false)

Validation Layer (pkg/fluid/*_test.go, pkg/fluid/*_bench_test.go)
  -> Baseline benchmarks
  -> Physics regression tests
  -> Interaction stability tests
  -> BFECC accuracy + stability tests
```

### Staggered Grid Layout (existing)

```
  ┌───V───┐
  │       │
  U   P   U     P = pressure, smoke (cell center)
  │   M   │     U = horizontal velocity (vertical edges)
  └───V───┘     V = vertical velocity (horizontal edges)
```

### Tool State Machine

```
  [W]all ─── [E]rase ─── [F]orce ─── [K]smo(k)e ─── s[O]urce ─── si[N]k
    │                                                                  │
    └──────────────────── cycles with Tab ─────────────────────────────┘

  Left-drag: applies active tool
  Right-drag: always applies force impulse (regardless of active tool)
```

### BFECC Algorithm

```
1. Forward advect:  φ̃ = advect(φ, u, dt)      // standard semi-Lagrangian
2. Backward advect: φ̂ = advect(φ̃, u, -dt)     // reverse trace
3. Error estimate:  e = (φ̂ - φ) / 2
4. Corrected:       φ* = φ - e                  // compensate
5. Final advect:    φ_new = advect(φ*, u, dt)   // advect corrected field
6. Clamp:           clamp to local min/max       // prevent overshoot
```

## Implementation Plan

### Phase 0: Baseline and Guardrails (~5%)

Establish measurable performance and correctness baselines before adding features.

**Files:**
- `pkg/fluid/fluid_bench_test.go` - Create: benchmark suite

**Tasks:**
- [ ] Add `BenchmarkSimulate` for 300x251 grid with standard parameters
- [ ] Add `BenchmarkSimulateWithJet` for jet injection + simulation
- [ ] Add helper to compute and report max divergence after simulation steps
- [ ] Record baseline numbers in test comments
- [ ] Verify all existing tests pass: `go test ./pkg/fluid/...`

### Phase 1: Interaction Tooling (~25%)

Add the tool system, force interaction, and smoke painting.

**Files:**
- `pkg/fluid/fluid.go` - Add `ApplyForce`, `ApplyForceRadius` methods
- `main/main.go` - Tool state machine, mouse interaction, brush parameters

**Tasks:**
- [ ] Add `ApplyForce(i, j int, fx, fy float32)` to `Fluid` - directly adds to U/V
- [ ] Add `ApplyForceRadius(cx, cy int, fx, fy float32, radius int)` - Gaussian falloff over radius
- [ ] Add tool mode enum and state to `Game` struct: `Wall`, `Erase`, `Force`, `Smoke`, `Source`, `Sink`
- [ ] Add keyboard tool switching: `W`=Wall, `E`=Erase, `F`=Force, `K`=Smoke, `O`=Source, `N`=Sink, `Tab`=cycle
- [ ] Remap left-drag to apply the active tool instead of always toggling walls
- [ ] Implement right-drag: always applies force impulse (compute velocity from mouse delta)
- [ ] Implement Force tool: left-drag pushes fluid along drag direction
- [ ] Implement Smoke tool: left-drag injects smoke at cursor with configurable density
- [ ] Implement Source tool: left-click places persistent velocity source region
- [ ] Implement Sink tool: left-click places persistent velocity sink region
- [ ] Add `[`/`]` keys to adjust brush radius
- [ ] Add `,`/`.` keys to adjust brush strength
- [ ] Preserve all existing controls (`S/J/V/A/Space/+/-/R/arrows/C`)
- [ ] Test: `ApplyForce` modifies velocity correctly
- [ ] Test: `ApplyForceRadius` applies Gaussian falloff
- [ ] Test: force injection doesn't break incompressibility (pressure solver corrects)

### Phase 2: Visualization and Diagnostics (~20%)

Add velocity magnitude and vorticity views, improve HUD, add streamlines.

**Files:**
- `pkg/fluid/fluid.go` - Add `Vorticity()`, `VelocityMagnitude()`, `MaxDivergence()` methods
- `main/main.go` - Viz mode cycling, HUD improvements, streamline rendering
- `main/colors.go` - Diverging colormap for vorticity

**Tasks:**
- [ ] Add `Vorticity() ScalarField` method to `Fluid` - computes curl at each cell center, returns populated `ScalarField`
- [ ] Add `VelocityMagnitude() ScalarField` method - computes |v| at cell centers
- [ ] Add `MaxDivergence() float32` method for HUD metrics
- [ ] Replace `showSmoke bool` with `vizMode` enum: Smoke=0, Pressure=1, VelMag=2, Vorticity=3
- [ ] Update `S` key to cycle through all 4 modes
- [ ] Update `Draw()` to render selected mode using appropriate colormap
- [ ] Add diverging colormap for vorticity in `colors.go` (blue negative, white zero, red positive)
- [ ] Add streamline rendering: when `A` (arrows) is active, draw short streamlines from seed points following velocity field
- [ ] Improve HUD: show active tool, brush radius/strength, viz mode, BFECC status, dt, max velocity, max divergence
- [ ] Add `H` key to toggle a help overlay listing all controls

### Phase 3: Particle Tracer System (~15%)

Add a Lagrangian particle system for colored non-diffusing tracer visualization.

**Files:**
- `main/main.go` - Particle struct, particle injection, advection, rendering
- `pkg/fluid/fluid.go` - Add `SampleVelocity(x, y float32) (u, v float32)` for particle advection

**Tasks:**
- [ ] Add `SampleVelocity(x, y float32) (float32, float32)` to `Fluid` - samples velocity at arbitrary position using `sampleField`
- [ ] Define `Particle` struct in main: `X, Y float32`, `Color color.RGBA`, `Age float32`, `MaxAge float32`
- [ ] Add particle slice to `Game` struct with configurable max count (default 10000)
- [ ] Add particle injection: when smoke tool is active, also spawn colored particles at cursor
- [ ] Add `P` key to toggle particle visibility
- [ ] Implement particle advection in `Update()`: for each particle, sample velocity and integrate position using RK2 (midpoint method)
- [ ] Remove dead particles (age > maxAge or position out of bounds or inside solid)
- [ ] Render particles as colored 2x2 pixel dots on screen
- [ ] Add color cycling for injected particles (hue rotation over time for visual variety)
- [ ] Add particle clear on `R` reset

### Phase 4: BFECC Advection (~20%)

Implement the Back and Forth Error Compensating and Correcting advection scheme.

**Files:**
- `pkg/fluid/fluid.go` - BFECC advection methods, `UseBFECC` flag
- `pkg/fluid/fluid_test.go` - BFECC accuracy and stability tests

**Tasks:**
- [ ] Add `UseBFECC bool` field to `Fluid` (default `false`)
- [ ] Implement `advectVelocityBFECC(dt)`:
  - Forward advect velocity using existing semi-Lagrangian
  - Backward advect the result with -dt
  - Compute error = (backward - original) / 2
  - Subtract error from original to get corrected field
  - Forward advect corrected field
  - Clamp to local min/max of 3x3 neighborhood
- [ ] Implement `advectSmokeBFECC(dt)`:
  - Same 3-pass approach for smoke field
  - Clamp to prevent negative smoke values
- [ ] Update `Simulate()` to use BFECC when `UseBFECC` is true
- [ ] Add `B` key in main.go to toggle BFECC on/off
- [ ] Update HUD to show BFECC status
- [ ] Test: advect a sharp smoke blob 50 steps with both methods - BFECC should have lower L2 diffusion error
- [ ] Test: BFECC stability - run 200 steps with jet + obstacles, verify no NaN/Inf
- [ ] Test: BFECC clamping prevents overshoot - verify smoke stays within [0, maxInjected]
- [ ] Benchmark: compare `BenchmarkSimulateBFECC` vs `BenchmarkSimulate` to measure cost ratio

### Phase 5: Presets, Polish, and Verification (~15%)

Add preset scenarios, finalize testing, update documentation.

**Files:**
- `main/main.go` - Preset loading, final parameter tuning keys
- `pkg/fluid/fluid.go` - Possibly `SetCircularObstacle` helper
- `pkg/fluid/fluid_test.go` - Final integration tests
- `README.md` - Complete controls documentation

**Tasks:**
- [ ] Add `SetCircularObstacle(cx, cy, radius int)` helper to `Fluid`
- [ ] Add number key presets:
  - `1`: Empty field with jet (current default)
  - `2`: Lid-driven cavity (top wall moves rightward)
  - `3`: Channel flow with circular obstacle (Karman vortex street)
- [ ] Implement lid-driven cavity: set top border velocity to constant rightward flow
- [ ] Ensure `R` resets to the current preset configuration
- [ ] Add `[,`/`.` keys for vorticity confinement strength adjustment
- [ ] Test: long-run stability with all presets (100+ steps, no NaN/Inf)
- [ ] Test: interaction API correctness with active sources/forces
- [ ] Test: smoke non-negativity and boundedness under all conditions
- [ ] Test: divergence reduction with active sources/forces
- [ ] Verify all existing tests still pass
- [ ] Update `README.md` with complete controls reference
- [ ] Final benchmark comparison: baseline vs all-features-enabled

## Control Compatibility Matrix

| Key/Mouse | Current Behavior | New Behavior | Compatibility |
|-----------|-----------------|--------------|---------------|
| `S` | Toggle smoke/pressure | Cycle: Smoke -> Pressure -> VelMag -> Vorticity | Extended (superset) |
| `J` | Toggle jet | Unchanged | Preserved |
| `V` | Toggle vorticity confinement | Unchanged | Preserved |
| `A` | Toggle velocity arrows | Toggle arrows + streamlines | Extended |
| `Space` | Pause/resume | Unchanged | Preserved |
| `+`/`-` | Speed adjust | Unchanged | Preserved |
| `R` | Reset | Reset to current preset | Extended |
| Arrows | Toggle walls | Unchanged | Preserved |
| `C` | Clear walls | Unchanged | Preserved |
| Left-drag | Toggle walls | Apply active tool (default: Wall) | Compatible (Wall is default tool) |
| **New** `W` | - | Select Wall tool | New |
| **New** `E` | - | Select Erase tool | New |
| **New** `F` | - | Select Force tool | New |
| **New** `K` | - | Select Smoke tool | New |
| **New** `O` | - | Select Source tool | New |
| **New** `N` | - | Select Sink tool | New |
| **New** `Tab` | - | Cycle tools | New |
| **New** Right-drag | - | Force impulse (always) | New |
| **New** `B` | - | Toggle BFECC | New |
| **New** `H` | - | Toggle help overlay | New |
| **New** `P` | - | Toggle particles | New |
| **New** `[`/`]` | - | Adjust brush radius | New |
| **New** `,`/`.` | - | Adjust brush strength | New |
| **New** `1-3` | - | Load preset scenario | New |

## Files Summary

| File | Action | Purpose |
|------|--------|---------|
| `pkg/fluid/fluid.go` | Modify | ApplyForce/ApplyForceRadius, Vorticity/VelocityMagnitude/MaxDivergence accessors, SampleVelocity, BFECC advection, UseBFECC flag, SetCircularObstacle |
| `pkg/fluid/fluid_test.go` | Modify | Force application tests, BFECC accuracy/stability tests, interaction stability tests, preset stability tests |
| `pkg/fluid/fluid_bench_test.go` | Create | Benchmark suite: Simulate, SimulateWithJet, SimulateBFECC |
| `main/main.go` | Modify | Tool state machine, force/smoke interaction, viz mode cycling, particle system, streamlines, HUD, presets, new key bindings |
| `main/colors.go` | Modify | Add diverging blue-white-red colormap for vorticity |
| `README.md` | Modify | Complete controls documentation |

## Definition of Done

**Interactivity:**
- [ ] Tool system works: Wall, Erase, Force, Smoke tools usable via keyboard selection + left-drag
- [ ] Right-drag always applies force impulse regardless of active tool
- [ ] Brush radius and strength adjustable at runtime
- [ ] All existing controls from README.md continue to work unchanged in their current semantics

**Physics/Correctness:**
- [ ] Pressure projection keeps divergence within threshold in test scenarios
- [ ] Smoke remains non-negative and bounded in long-running tests
- [ ] Dynamic obstacle edits don't leak velocity through newly solid cells
- [ ] BFECC advection produces measurably less diffusion than standard SL in automated test

**Stability:**
- [ ] No NaN/Inf in 200-step stress tests with mixed interactions
- [ ] Jet + obstacle + active force interaction remains stable

**Visualization:**
- [ ] 4 viz modes cycle with `S` key
- [ ] Vorticity uses diverging colormap
- [ ] Particle tracer injects, advects, and renders colored particles
- [ ] HUD shows active tool, viz mode, BFECC status, brush params

**Performance:**
- [ ] Go benchmarks exist for `Simulate` and `SimulateBFECC`
- [ ] Interactive scenarios maintain real-time behavior on development hardware
- [ ] BFECC cost ratio documented in benchmark output

**Testing:**
- [ ] `go test ./...` passes
- [ ] New tests cover: force API, BFECC accuracy, BFECC stability, interaction stability
- [ ] Benchmarks added and documented

## Verification Strategy

**Automated (go test):**
- Unit tests for `ApplyForce`, `ApplyForceRadius`, `SampleVelocity` correctness
- Physics validation: divergence-free after pressure projection with active forces
- BFECC accuracy: advect sharp blob, compare diffusion metrics against standard SL
- BFECC stability: 200 steps with jet + obstacles + BFECC, verify no NaN/Inf
- Smoke boundedness: verify non-negative after long runs with all features enabled
- Interaction stability: rapid wall edits + forces + sources don't crash or diverge

**Benchmarks (go test -bench):**
- `BenchmarkSimulate` - baseline frame time at 300x251
- `BenchmarkSimulateWithJet` - frame time with jet injection
- `BenchmarkSimulateBFECC` - frame time with BFECC enabled
- Document cost ratios (BFECC/standard) in test output

**Manual:**
- Walk through all key bindings per control compatibility matrix
- Visual validation of each viz mode with known flow patterns
- FPS check under high interaction load (continuous force + smoke injection)
- Particle tracer visual inspection: particles follow flow, don't accumulate in solids

## Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| BFECC 3x advection cost exceeds frame budget | Medium | High | Profile early in Phase 4; BFECC defaults to off; can reduce grid or parallelize BFECC passes |
| BFECC overshoot creates instability | Medium | Medium | Clamping limiter (local min/max of neighbors); extensive stability tests |
| Force injection too strong/weak for good UX | Low | Low | Configurable scaling; `[`/`]` and `,`/`.` for runtime tuning |
| Tool-mode model confuses existing users | Low | Low | Default tool is Wall (preserves current behavior); help overlay with `H` |
| Particle count impacts performance | Medium | Medium | Cap at 10000; remove old particles; particles are render-only (no physics feedback) |
| main.go grows unwieldy | Medium | Low | Can extract presets.go and input.go if needed; keep it simple for now |

## Dependencies

- Go 1.24.2
- Ebiten v2.8.7 (already vendored)
- Reference: Mick West's BFECC C++ implementation in `pkg/fluid/example/example.cpp`
- No external sprint dependencies (Sprint 001)

## Open Questions

1. Should BFECC be applied to both velocity and smoke, or smoke-only first? (Plan assumes both, toggled together)
2. Should sources/sinks be persistent objects stored in `Fluid`, or applied from UI each frame? (Plan assumes UI-driven per-frame for simplicity)
3. Particle advection uses RK2 - is simple Euler integration sufficient for visual purposes?
4. Should the particle color scheme be user-configurable, or is automatic hue rotation sufficient?

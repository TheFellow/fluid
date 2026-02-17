# Sprint 001: Full Interactive 2D Fluid Simulator

## Overview

This sprint takes the existing 2D Eulerian fluid simulator from a working prototype to a complete, interactive fluid simulation. The core solver is already sound - staggered grid with Jacobi pressure projection and semi-Lagrangian advection - but three critical gaps remain: (1) the branch is named `bfecc` yet the actual BFECC (Back and Forth Error Compensating and Correcting) advection scheme is not implemented, (2) interaction is limited to wall drawing and a fixed jet, and (3) visualization offers only two modes (smoke/pressure) with no way to see the velocity or vorticity structure.

The plan addresses these in three phases: first upgrade the advection to true BFECC for better accuracy, then add rich mouse-driven interaction (force brushes, smoke painting, drag-to-push), and finally add visualization modes and a heads-up display for real-time parameter tuning. Throughout, we maintain the existing test suite and add targeted tests for the new features.

## Use Cases

1. **Force interaction**: User right-click-drags on the fluid to push it, creating vortices and turbulence interactively
2. **Smoke painting**: User middle-click or modifier-click to inject smoke at the cursor, painting tracer patterns into the flow
3. **BFECC advection**: Simulation produces sharper, less diffusive fluid features compared to basic semi-Lagrangian
4. **Multi-mode visualization**: User cycles through smoke, pressure, velocity magnitude, and vorticity views to understand the flow
5. **Parameter tuning**: User adjusts viscosity, confinement strength, iteration count, and other parameters in real time via keyboard shortcuts
6. **Preset scenarios**: User loads preconfigured scenarios (lid-driven cavity, channel flow with obstacles, Karman vortex street) to explore different physics

## Architecture

```
                    ┌─────────────────┐
                    │   main/main.go  │
                    │   Game struct   │
                    │ ┌─────────────┐ │
                    │ │  Interaction│ │ ← mouse force, smoke painting
                    │ │  Manager    │ │
                    │ ├─────────────┤ │
                    │ │ Viz Modes   │ │ ← smoke, pressure, velocity, vorticity
                    │ ├─────────────┤ │
                    │ │ HUD/Overlay │ │ ← parameter display, controls help
                    │ └─────────────┘ │
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
                    │  pkg/fluid/     │
                    │  Fluid struct   │
                    │ ┌─────────────┐ │
                    │ │ Simulate()  │ │
                    │ │  viscosity  │ │
                    │ │  pressure   │ │
                    │ │  confine    │ │
                    │ │  borders    │ │
                    │ │  advect     │◄├── BFECC upgrade (advectVelocityBFECC)
                    │ │  smoke      │ │
                    │ ├─────────────┤ │
                    │ │ Force API   │ │ ← AddForce(i,j, fx,fy)
                    │ │ Vorticity   │ │ ← Vorticity() ScalarField
                    │ │ VelMag      │ │ ← VelocityMagnitude() ScalarField
                    │ └─────────────┘ │
                    └─────────────────┘
```

### BFECC Algorithm

The BFECC method improves semi-Lagrangian advection accuracy from first-order to second-order:

```
1. Forward advect:  φ̃ = advect(φ, u, dt)      // standard semi-Lagrangian
2. Backward advect: φ̂ = advect(φ̃, u, -dt)     // reverse trace
3. Error estimate:  e = (φ̂ - φ) / 2
4. Corrected:       φ* = φ - e                  // compensate
5. Final advect:    φ_new = advect(φ*, u, dt)   // advect corrected field
```

This requires 3x the advection cost but produces much sharper results. The existing `sampleField` and `advectVelocity`/`advectSmoke` functions provide the building blocks.

### Staggered Grid (existing)

```
  ┌───V───┐
  │       │
  U   P   U
  │   M   │
  └───V───┘

P = pressure (cell center)
M = smoke (cell center)
U = horizontal velocity (vertical cell edges)
V = vertical velocity (horizontal cell edges)
```

## Implementation Plan

### Phase 1: BFECC Advection (~30%)

Implement the Back and Forth Error Compensating and Correcting advection scheme that the branch is named after.

**Files:**
- `pkg/fluid/fluid.go` - Add BFECC advection methods, add `UseBFECC bool` field
- `pkg/fluid/fluid_test.go` - Add BFECC accuracy tests comparing against standard advection

**Tasks:**
- [ ] Add `UseBFECC bool` field to `Fluid` struct (default `true`)
- [ ] Implement `advectVelocityBFECC(dt)` method using the 3-pass BFECC algorithm
- [ ] Implement `advectSmokeBFECC(dt)` method for smoke transport
- [ ] Add clamping limiter to BFECC to prevent overshoots (clamp corrected values to local min/max of neighboring cells)
- [ ] Update `Simulate()` to use BFECC methods when `UseBFECC` is true
- [ ] Add `B` key binding in `main.go` to toggle BFECC on/off for comparison
- [ ] Add tests: compare advection of a sharp smoke blob - BFECC should preserve shape better than standard SL
- [ ] Add tests: verify BFECC doesn't introduce instability (run 100 steps, check no NaN/Inf)

### Phase 2: Interactive Force & Smoke Application (~30%)

Add mouse-driven force application and smoke injection so the user can push fluid around and paint smoke.

**Files:**
- `pkg/fluid/walls.go` - Add `AddForce(i, j int, fx, fy float32)` method
- `main/main.go` - Add right-click drag for forces, shift+click for smoke

**Tasks:**
- [ ] Add `AddForce(i, j int, fx, fy float32)` to `Fluid` - directly modifies U/V at the given cell
- [ ] Add `ApplyForceRadius(i, j int, fx, fy, radius float32)` - applies force with Gaussian falloff over a radius
- [ ] Implement right-click drag in `Game.Update()`: compute velocity from mouse delta, apply as force to fluid
- [ ] Implement Shift+left-click drag to inject smoke at cursor position with configurable density
- [ ] Track previous mouse position for velocity computation (dx/dy per frame)
- [ ] Add force scaling parameter to Game struct
- [ ] Add smoke brush size parameter
- [ ] Test: verify `AddForce` modifies velocity correctly
- [ ] Test: verify force application doesn't break incompressibility (pressure solver corrects it)

### Phase 3: Visualization & Display Modes (~25%)

Add velocity magnitude and vorticity visualization modes, improve the HUD.

**Files:**
- `pkg/fluid/fluid.go` - Add `Vorticity() ScalarField` and `VelocityMagnitude() ScalarField` methods
- `pkg/fluid/vorticity.go` - New file: vorticity field accessor (like pressure.go/smoke.go)
- `pkg/fluid/velocity_magnitude.go` - New file: velocity magnitude field accessor
- `main/main.go` - Cycle through 4 viz modes, improve HUD display
- `main/colors.go` - Add diverging colormap for vorticity (blue-white-red)

**Tasks:**
- [ ] Create `VorticityField` struct implementing `ScalarField` - computes curl at each cell center
- [ ] Create `VelocityMagnitudeField` struct implementing `ScalarField` - computes |v| at each cell center
- [ ] Add `Vorticity()` and `VelocityMagnitude()` methods to `Fluid`
- [ ] Add diverging colormap for vorticity in `colors.go` (negative=blue, zero=white, positive=red)
- [ ] Change `showSmoke bool` to `vizMode int` enum (Smoke, Pressure, Velocity, Vorticity)
- [ ] Update `S` key to cycle through all 4 modes
- [ ] Update `Draw()` to render the selected visualization mode
- [ ] Update HUD to show current mode, BFECC status, and all keyboard shortcuts
- [ ] Add `H` key to toggle help overlay listing all controls

### Phase 4: Polish & Scenarios (~15%)

Add preset scenarios and final refinements.

**Files:**
- `main/main.go` - Add scenario loading, parameter adjustment keys
- `pkg/fluid/fluid.go` - Add `Reset()` improvements if needed

**Tasks:**
- [ ] Add `1`-`4` number keys to load preset scenarios:
  - `1`: Empty field with jet (current default)
  - `2`: Lid-driven cavity (top wall moves right)
  - `3`: Channel flow with circular obstacle (Karman vortex street)
  - `4`: Double slit diffraction-like pattern
- [ ] Add `[`/`]` keys to adjust vorticity confinement strength
- [ ] Add `,`/`.` keys to adjust viscosity
- [ ] Implement lid-driven cavity: set top border velocity instead of solid
- [ ] Implement circular obstacle helper: `SetCircularObstacle(cx, cy, radius int)`
- [ ] Ensure `R` reset clears forces and restores the current scenario
- [ ] Final pass: verify all keyboard shortcuts in README.md match implementation

## Files Summary

| File | Action | Purpose |
|------|--------|---------|
| `pkg/fluid/fluid.go` | Modify | Add BFECC advection, UseBFECC flag, force API, vorticity/velocity-magnitude accessors |
| `pkg/fluid/fluid_test.go` | Modify | Add BFECC accuracy tests, force application tests |
| `pkg/fluid/walls.go` | Modify | Add AddForce, ApplyForceRadius methods |
| `pkg/fluid/vorticity.go` | Create | VorticityField implementing ScalarField |
| `pkg/fluid/velocity_magnitude.go` | Create | VelocityMagnitudeField implementing ScalarField |
| `main/main.go` | Modify | BFECC toggle, force interaction, smoke painting, viz mode cycling, scenarios, HUD |
| `main/colors.go` | Modify | Add diverging colormap for vorticity |
| `README.md` | Modify | Update controls documentation |

## Definition of Done

- [ ] BFECC advection implemented and togglable (`B` key)
- [ ] Right-click drag applies force to fluid, creating visible motion
- [ ] Shift+left-click drag injects smoke at cursor
- [ ] 4 visualization modes cycle with `S` key (Smoke, Pressure, Velocity Magnitude, Vorticity)
- [ ] Vorticity uses diverging blue-white-red colormap
- [ ] At least 2 preset scenarios loadable with number keys
- [ ] All existing tests pass
- [ ] New tests: BFECC accuracy, BFECC stability, force application correctness
- [ ] Simulation stays stable at 60 FPS with BFECC enabled at 300x251 resolution
- [ ] README.md reflects all current controls
- [ ] No compiler warnings

## Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| BFECC 3x advection cost drops below 60 FPS | Medium | High | Profile early; can parallelize BFECC passes; can reduce grid size if needed |
| BFECC overshoot creates instability | Medium | Medium | Implement clamping limiter (clamp to local min/max of neighbors) |
| Mouse force injection too strong/weak | Low | Low | Add configurable force scaling; tune during testing |
| Vorticity colormap hard to read | Low | Low | Use established diverging colormap; test with known vortex patterns |

## Security Considerations

- No network I/O or user-supplied data beyond keyboard/mouse input
- No file I/O beyond loading the executable
- No security concerns for this project

## Dependencies

- No external sprint dependencies (this is Sprint 001)
- Depends on existing Ebiten v2.8.7 (already vendored)
- Reference: Mick West's BFECC C++ implementation in `pkg/fluid/example/example.cpp`

## Open Questions

1. Should BFECC be applied to both velocity and smoke, or just velocity? (Draft assumes both, with independent toggles possible)
2. Is the 300x251 grid resolution sufficient for the BFECC performance budget? May need profiling.
3. Should we implement actual BFECC from the Mick West reference, or the MacCormack variant? (Draft assumes standard BFECC)
4. Should multigrid be completed as part of this sprint, or deferred? (Draft defers it - Jacobi works fine)

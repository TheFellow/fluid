# Sprint 001 Intent: Full Interactive 2D Fluid Simulator

## Seed

"I'm building a full fluid simulator. This should be an interactive and full implementation of a 2d sim that let's me interact with the flow. Be sure it is tested in whatever ways are sufficient."

## Context

This is a **greenfield sprint** - the first formal sprint for an existing codebase. The project already has a working 2D Eulerian fluid simulator in Go using Ebiten for rendering, on the `bfecc` branch with recent stability improvements. The simulator has a semi-Lagrangian advection scheme, Jacobi pressure projection with adaptive relaxation, vorticity confinement, and several visual enhancement features.

The user wants to take this from a "working prototype" to a "full interactive 2D fluid simulator." This means closing the gaps between the current implementation and a polished, feature-complete 2D fluid sim with rich interaction capabilities.

## Recent Sprint Context

No formal sprints exist yet. The project history shows:
- Initial implementation of the Eulerian solver
- BFECC branch created (referencing Mick West's Back and Forth Error Compensating and Correcting method)
- Recent commits focused on stability improvements: reduced pressure iterations, adaptive relaxation, early termination, multigrid infrastructure
- Comprehensive test suite added (1233 lines)

## Relevant Codebase Areas

### Core Simulation
- `pkg/fluid/fluid.go` (776 lines) - Main solver: `Simulate()`, pressure projection (Jacobi + multigrid infra), semi-Lagrangian advection, vorticity confinement, viscosity diffusion, turbulence, boundary handling
- `pkg/fluid/fluid_test.go` (1233 lines) - Comprehensive tests: boundary ops, diffusion, advection, pressure projection, stability, realistic scenarios

### Field Interfaces
- `pkg/fluid/scalar_field.go` - ScalarField interface for pressure/smoke output
- `pkg/fluid/vector_field.go` - VectorField interface for velocity output
- `pkg/fluid/velocity.go`, `pressure.go`, `smoke.go` - Field accessors
- `pkg/fluid/walls.go` - Solid boundary handling (SetSolid, IsSolid, SetVelocity, AddSmoke)

### Visualization & Interaction
- `main/main.go` - Ebiten game loop, rendering, keyboard/mouse controls
- `main/colors.go` - Scientific color mapping (blue->cyan->green->yellow->red)

### Reference
- `pkg/fluid/example/example.cpp` - Mick West's BFECC reference implementation (C++)

### Utilities
- `pkg/fluid/parallel.go` - `parallelRange()` for parallel computation

## Constraints

- Go 1.24.2 with Ebiten v2.8.7
- No CLAUDE.md conventions file exists yet - we establish conventions with this sprint
- Must maintain backward compatibility with existing keyboard/mouse controls
- Must keep simulation real-time at 60 FPS on desktop hardware
- Should build on existing staggered grid architecture
- Tests should be meaningful physics validation, not just unit tests

## Success Criteria

A "full interactive 2D fluid simulator" means:

1. **Physically accurate**: The simulation produces visually and numerically correct fluid behavior (incompressibility, correct advection, proper boundary handling)
2. **Rich interaction**: User can interact with the fluid in intuitive ways beyond just wall drawing - force application, source/sink placement, parameter tuning
3. **Visual quality**: Rendering is informative and aesthetically pleasing with multiple visualization modes
4. **Stability**: Simulation remains stable under all interactive conditions
5. **Performance**: Maintains 60 FPS at reasonable grid resolutions
6. **Well-tested**: Sufficient test coverage to catch regressions in physics accuracy and stability

## Verification Strategy

- **Physics validation**: Tests that verify divergence-free velocity fields, conservation properties, correct boundary conditions
- **Stability testing**: Long-running simulations with various interaction patterns don't diverge
- **Regression tests**: Key scenarios produce consistent results
- **Visual inspection**: Manual testing of all interaction modes
- **Performance benchmarks**: Frame time measurements under load

## Uncertainty Assessment

- Correctness uncertainty: **Low** - The solver is already working; improvements are incremental
- Scope uncertainty: **High** - "Full fluid simulator" is ambiguous; need to bound what "full" means
- Architecture uncertainty: **Low** - Existing staggered grid architecture is sound and extensible

## Open Questions

1. What does "full" mean? What features are must-have vs nice-to-have?
2. Should we implement the actual BFECC advection scheme (the branch name suggests it, but the current code uses basic semi-Lagrangian)?
3. What interaction modes are desired beyond wall drawing? (force brushes, smoke sources, temperature/buoyancy?)
4. Should we add a GUI for parameter tuning or keep keyboard shortcuts?
5. Is multigrid solver completion in scope, or is Jacobi sufficient?
6. Are there specific physical scenarios the user wants to simulate? (lid-driven cavity, flow past obstacles, channel flow?)

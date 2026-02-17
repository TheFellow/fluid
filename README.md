# Fluid

Interactive 2D fluid simulator using a staggered-grid Eulerian solver with semi-Lagrangian advection and optional BFECC error correction.

## Visualization

- `S`: Cycle visualization mode: Smoke → Pressure → Velocity Magnitude → Vorticity
- `A`: Toggle velocity arrow overlay
- `P`: Toggle particle tracer visibility
- `H`: Toggle help overlay

## Simulation

- `J`: Toggle the jet inlet
- `V`: Toggle vorticity confinement
- `B`: Toggle BFECC advection (sharper, less diffusive)
- `Space`: Pause or resume the simulation
- `+`/`-`: Increase or decrease simulation speed
- `R`: Reset to the current preset

## Tools

Select tools with keyboard, apply with left-click drag:

- `W`: Wall tool — draw solid walls
- `E`: Erase tool — remove walls
- `F`: Force tool — push fluid along drag direction
- `K`: Smoke tool — inject smoke (also spawns tracer particles)
- `O`: Source tool — place persistent velocity source
- `N`: Sink tool — place persistent velocity sink
- `Tab`: Cycle through tools

## Brush

- `[`/`]`: Decrease/increase brush radius
- `,`/`.`: Decrease/increase brush strength

## Mouse

- **Left-click drag**: Apply the active tool
- **Right-click drag**: Always apply force impulse (regardless of active tool)

## Walls

- Arrow keys `Up`, `Down`, `Left`, `Right`: Toggle the corresponding border wall
- `C`: Remove all walls except the active border walls

## Presets

- `1`: Jet flow (default)
- `2`: Lid-driven cavity
- `3`: Karman vortex street (flow past circular obstacle)

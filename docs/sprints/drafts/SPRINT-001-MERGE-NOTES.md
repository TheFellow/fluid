# Sprint 001 Merge Notes

## Claude Draft Strengths
- Strong BFECC algorithm description with pseudocode
- Concrete task checklists per phase
- Good ASCII architecture diagram
- Specific key bindings proposed for each feature
- Risk table with likelihood/impact/mitigation columns

## Codex Draft Strengths
- Phase 0 (Baseline and Guardrails) is an excellent idea - establish metrics before adding features
- Tool-mode model (Wall/EraseWall/Force/Smoke/Source/Sink) is cleaner than modifier-key approach
- Explicit scope boundaries (in/out) section prevents scope creep
- Separate Verification Strategy section with automated/manual split
- Better phase ordering: interaction first, BFECC as enhancement
- Suggested extracting presets.go and input.go to keep main.go manageable
- Divergence magnitude as a viz mode (Claude missed this)
- "Recommended Acceptance Slice" for scope tightening is pragmatic

## Valid Critiques Accepted

1. **Phase ordering** - Codex and interview both confirm: interaction first, BFECC later. Resequenced.
2. **ScalarField is concrete, not an interface** - Codex is correct. Vorticity/VelMag methods should return `ScalarField` directly, not new implementing types. No new files needed.
3. **BFECC default should be false** - Agree. Safer to opt-in until baseline is measured.
4. **AddForce doesn't belong in walls.go** - Agree. Better in fluid.go or a dedicated interaction.go.
5. **"No compiler warnings" is Go-irrelevant** - Removed from DoD.
6. **Missing verification strategy section** - Added.
7. **Tool-mode model is better than Shift+click** - Accepted. Tool modes with keyboard switching.
8. **Phase 0 baseline** - Accepted. Measure before changing.
9. **Double slit scenario is low value** - Replaced with channel flow with obstacle (Karman vortex street).

## Critiques Rejected (with reasoning)

1. **"Security Considerations section adds little value"** - Removed entirely rather than keeping a placeholder. Agree it's not needed for a local desktop app.
2. **Divergence as a separate render mode** - Keeping this from Codex draft, but it's essentially what the pressure view already shows. Will combine as "Pressure/Divergence" mode.

## Interview Refinements Applied

1. **Interaction is top priority** - Phases reordered: Interaction -> Viz/Particles -> BFECC -> Presets
2. **All interaction modes wanted** - Force, smoke, scenarios, param tuning all included
3. **Streamlines AND particles wanted** - Added Phase 3 with particle system and streamline visualization, including "colored non-diffusing smoke"
4. **Benchmarks + physics tests** - Added Phase 0 baseline + dedicated verification strategy section with Go benchmark requirements

## Final Decisions

- 6 phases: Baseline -> Interaction -> Solver/API -> Visualization+Particles -> BFECC -> Presets+Polish
- Tool-mode model with keyboard switching (not modifier keys)
- ScalarField returned directly from new methods (Vorticity, VelocityMagnitude)
- BFECC defaults to false, toggled with B key
- Particle system for tracer visualization (colored non-diffusing particles)
- Streamline rendering in velocity arrow overlay mode
- Go benchmarks in Phase 0 and Phase 5

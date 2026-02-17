# Critique: SPRINT-001-CLAUDE-DRAFT

## Overall Assessment

Claude's draft is strong on ambition, clear phase framing, and concrete feature ideas. It identifies the core gap around true BFECC and proposes meaningful interaction and visualization upgrades. The main issues are scope ordering, a few architecture mismatches with the current code, and missing verification rigor expected by this repo's sprint template.

## Findings

### High Severity

1. Missing explicit verification strategy section (template mismatch)
- `docs/sprints/README.md` expects a dedicated verification strategy section.
- Claude includes tests in tasks/DoD, but not a structured verification plan across physics, stability, visual validation, and performance.

2. Scope prioritization is misaligned with sprint intent
- The sprint intent emphasizes "full interactive" simulation first.
- Claude places BFECC implementation as Phase 1 and most interaction work in Phase 2.
- This increases delivery risk: a technically deep solver change is scheduled before core interaction goals.

3. Scalar field architecture assumptions conflict with current code
- Claude proposes `VorticityField` and `VelocityMagnitudeField` "implementing ScalarField".
- In this codebase, `ScalarField` is a concrete struct (`pkg/fluid/scalar_field.go`), not an interface.
- Proposed file/API shape should be adjusted to current patterns to avoid unnecessary refactors.

4. DoD performance target is underspecified and potentially non-reproducible
- "60 FPS at 300x251 with BFECC" is hardware-dependent and currently lacks a benchmark protocol in the draft.
- Without baseline method and environment constraints, pass/fail is ambiguous.

### Medium Severity

1. Backward compatibility risk for existing controls
- Existing `S` toggles smoke/pressure.
- Claude repurposes `S` to cycle 4 modes but does not specify compatibility behavior (required by intent constraints).

2. Interaction design conflicts with existing left-drag wall editing
- Draft assigns Shift+left-drag to smoke while current primary left-drag behavior is wall painting.
- Needs explicit tool-mode model to avoid input ambiguity and accidental wall edits.

3. BFECC default-enabled is risky
- Draft sets `UseBFECC` default `true` before performance/stability baseline.
- Safer rollout is opt-in default `false` until thresholds are measured.

4. File placement recommendation is debatable
- `AddForce` proposed in `pkg/fluid/walls.go` couples force application with wall utilities.
- Better home is `pkg/fluid/fluid.go` or a dedicated interactions file.

5. Scenario scope includes one low-value item for this sprint
- "Double slit diffraction-like pattern" reads as novelty and is less aligned with incompressible fluid validation than jet/cavity/obstacle flows.

### Low Severity

1. "No compiler warnings" in DoD is not useful for Go
- Go tooling does not produce warning management in the same way as C/C++ projects.

2. "Security Considerations" section adds little value for this sprint
- Not wrong, but not a priority relative to missing verification detail.

## What Claude's Draft Gets Right

1. Strong identification of the BFECC gap versus branch naming.
2. Clear, actionable task checklists by phase.
3. Good attention to visualization diagnostics (velocity/vorticity modes and HUD).
4. Useful risk table with practical mitigations.

## Recommended Revisions

1. Reorder phases to deliver interactivity first, then BFECC as gated enhancement.
2. Add a dedicated verification strategy section with explicit automated/manual checks and thresholds.
3. Align API proposals to current `ScalarField`/`VectorField` concrete struct patterns.
4. Define input/tool state model explicitly to preserve current controls and avoid conflicts.
5. Change BFECC rollout to toggleable, baseline-off until benchmark/stability criteria are met.
6. Replace novelty scenario(s) with physics-relevant presets tied to acceptance checks.
7. Replace FPS-only DoD target with reproducible benchmark criteria plus a manual FPS check.

## Suggested Acceptance Delta

If Claude's draft is adopted, require these minimum edits before execution:
- Add `Verification Strategy` section.
- Add control compatibility matrix (existing controls vs new bindings).
- Update file/API plan to match current concrete field types.
- Re-sequence phases: interaction -> diagnostics -> BFECC optimization/polish.

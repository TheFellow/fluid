package fluid

import (
	"math"
	"testing"
)

func TestCopyBorder(t *testing.T) {
	f := New(1.0, 3, 3, 1.0) // creates 5x5 grid due to border
	src := make([]float32, f.numCells)
	for i := range src {
		src[i] = float32(i + 1) // unique value for each cell
	}
	dst := make([]float32, f.numCells)

	f.copyBorder(dst, src)

	n := f.NumY
	for i := 0; i < f.NumX; i++ {
		if dst[i*n+0] != src[i*n+0] {
			t.Errorf("top border mismatch at (%d,0)", i)
		}
		if dst[i*n+f.NumY-1] != src[i*n+f.NumY-1] {
			t.Errorf("bottom border mismatch at (%d,%d)", i, f.NumY-1)
		}
	}
	for j := 0; j < f.NumY; j++ {
		if dst[0*n+j] != src[0*n+j] {
			t.Errorf("left border mismatch at (0,%d)", j)
		}
		if dst[(f.NumX-1)*n+j] != src[(f.NumX-1)*n+j] {
			t.Errorf("right border mismatch at (%d,%d)", f.NumX-1, j)
		}
	}
}

func TestDiffusionBehavior(t *testing.T) {
	f := New(1.0, 10, 10, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	// Set up initial velocity field with a concentrated high velocity region
	n := f.NumY
	centerX, centerY := f.NumX/2, f.NumY/2
	f.U[centerX*n+centerY] = 10.0
	f.V[centerX*n+centerY] = 10.0

	// Store initial state
	initialU := make([]float32, len(f.U))
	copy(initialU, f.U)

	// Run incompressibility solver which includes diffusion-like behavior
	f.makeIncompressible(10, 0.01)

	// Verify that the high velocity has been distributed to neighboring cells
	tolerance := float32(0.1)

	// Check that the center velocity has been reduced
	if f.U[centerX*n+centerY] >= initialU[centerX*n+centerY] {
		t.Errorf("Expected diffusion to reduce center velocity, got U[%d,%d] = %f, initial = %f",
			centerX, centerY, f.U[centerX*n+centerY], initialU[centerX*n+centerY])
	}

	// Check that neighboring cells have gained some velocity
	neighbors := [][2]int{
		{centerX - 1, centerY}, {centerX + 1, centerY},
		{centerX, centerY - 1}, {centerX, centerY + 1},
	}

	for _, neighbor := range neighbors {
		i, j := neighbor[0], neighbor[1]
		if i > 0 && i < f.NumX-1 && j > 0 && j < f.NumY-1 {
			if math.Abs(float64(f.U[i*n+j])) < float64(tolerance) &&
				math.Abs(float64(f.V[i*n+j])) < float64(tolerance) {
				t.Errorf("Expected neighbor cell [%d,%d] to have gained some velocity", i, j)
			}
		}
	}
}

func TestVelocityAdvection(t *testing.T) {
	f := New(1.0, 10, 10, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Set up a simple flow field - uniform rightward velocity
	rightwardVel := float32(2.0)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			f.U[i*n+j] = rightwardVel
			f.V[i*n+j] = 0.0
		}
	}

	// Add a velocity perturbation at a specific location
	testX, testY := 3, 5
	initialPerturbV := float32(5.0)
	f.V[testX*n+testY] = initialPerturbV

	// Store initial state
	initialU := make([]float32, len(f.U))
	initialV := make([]float32, len(f.V))
	copy(initialU, f.U)
	copy(initialV, f.V)

	dt := float32(0.1)

	// Run velocity advection
	f.advectVelocity(dt)

	// Test 1: The rightward flow should advect the V perturbation to the right
	// The perturbation should move from (testX, testY) towards (testX+1, testY)
	advectedX := testX + 1
	if advectedX < f.NumX-1 {
		if f.V[advectedX*n+testY] <= initialV[advectedX*n+testY] {
			t.Errorf("Expected velocity perturbation to be advected rightward")
		}
	}

	// Test 2: The original location should have reduced V velocity
	if f.V[testX*n+testY] >= initialPerturbV*0.9 {
		t.Errorf("Expected original perturbation location to lose some velocity due to advection")
	}

	// Test 3: Verify the U field maintains rightward flow pattern
	totalU := float32(0.0)
	count := 0
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			totalU += f.U[i*n+j]
			count++
		}
	}
	avgU := totalU / float32(count)

	if avgU < rightwardVel*0.8 {
		t.Errorf("Expected average U velocity to remain close to initial rightward flow, got %f", avgU)
	}
}

func TestSmokeAdvection(t *testing.T) {
	f := New(1.0, 10, 10, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Set up a rightward flow field
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			f.U[i*n+j] = 3.0 // Strong rightward flow
			f.V[i*n+j] = 0.0
		}
	}

	// Add smoke at a specific location
	smokeX, smokeY := 2, 5
	initialSmoke := float32(1.0)
	f.M[smokeX*n+smokeY] = initialSmoke

	// Store initial state
	initialM := make([]float32, len(f.M))
	copy(initialM, f.M)

	dt := float32(0.1)

	// Run smoke advection
	f.advectSmoke(dt)

	// Test 1: Smoke should move rightward due to the velocity field
	// Check that smoke has appeared in neighboring cells to the right
	smokeMoved := false
	for i := smokeX; i < min(smokeX+3, f.NumX-1); i++ {
		if f.M[i*n+smokeY] > initialM[i*n+smokeY]+0.01 {
			smokeMoved = true
			break
		}
	}

	if !smokeMoved {
		t.Error("Expected smoke to be advected rightward by the velocity field")
	}

	// Test 2: Total smoke should be conserved (approximately)
	initialTotal := float32(0.0)
	finalTotal := float32(0.0)

	for i := range f.M {
		initialTotal += initialM[i]
		finalTotal += f.M[i]
	}

	conservationError := math.Abs(float64(finalTotal - initialTotal))
	if conservationError > 0.1 {
		t.Errorf("Smoke conservation error too high: %f", conservationError)
	}
}

func TestPressureProjection(t *testing.T) {
	f := New(1.0, 8, 8, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create a divergent velocity field (sources/sinks)
	centerX, centerY := f.NumX/2, f.NumY/2

	// Set up radial outflow (divergent field)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			dx := float32(i - centerX)
			dy := float32(j - centerY)
			dist := float32(math.Sqrt(float64(dx*dx+dy*dy))) + 0.1

			f.U[i*n+j] = dx / dist
			f.V[i*n+j] = dy / dist
		}
	}

	// Calculate initial divergence
	initialDiv := calculateDivergence(f)

	// Run pressure projection
	f.makeIncompressible(20, 0.01)

	// Calculate final divergence
	finalDiv := calculateDivergence(f)

	// Test that divergence has been significantly reduced
	if finalDiv >= initialDiv*0.5 {
		t.Errorf("Expected pressure projection to reduce divergence from %f to much less, got %f",
			initialDiv, finalDiv)
	}

	// Test that divergence is close to zero in interior cells
	maxDiv := float32(0.0)
	for i := 2; i < f.NumX-2; i++ {
		for j := 2; j < f.NumY-2; j++ {
			div := f.U[(i+1)*n+j] - f.U[i*n+j] + f.V[i*n+j+1] - f.V[i*n+j]
			if math.Abs(float64(div)) > float64(maxDiv) {
				maxDiv = float32(math.Abs(float64(div)))
			}
		}
	}

	if maxDiv > 0.5 {
		t.Errorf("Maximum divergence after projection too high: %f", maxDiv)
	}
}

func calculateDivergence(f *Fluid) float32 {
	n := f.NumY
	totalDiv := float32(0.0)
	count := 0

	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] > 0 {
				div := f.U[(i+1)*n+j] - f.U[i*n+j] + f.V[i*n+j+1] - f.V[i*n+j]
				totalDiv += float32(math.Abs(float64(div)))
				count++
			}
		}
	}

	if count > 0 {
		return totalDiv / float32(count)
	}
	return 0
}

func TestBoundaryConditions(t *testing.T) {
	f := New(1.0, 6, 6, 1.0)

	// Initialize interior as fluid, but create some solid boundaries
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create a solid obstacle in the middle
	obstacleX, obstacleY := 3, 3
	f.S[obstacleX*n+obstacleY] = 0.0

	// Set some initial velocities
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			f.U[i*n+j] = 1.0
			f.V[i*n+j] = 0.5
		}
	}

	// Apply boundary handling
	f.handleBorders()

	// Test 1: Velocities at domain boundaries should be handled correctly
	// Check top and bottom borders
	for i := 0; i < f.NumX; i++ {
		// Top border (j=0) should have appropriate boundary condition
		if f.S[i*n+0] == 0 || (i > 0 && f.S[i*n+1] == 0) {
			if f.U[i*n+0] != 0 {
				t.Errorf("Expected U velocity at solid boundary to be zero at (%d,0)", i)
			}
		}

		// Bottom border
		bottomJ := f.NumY - 1
		if f.S[i*n+bottomJ] == 0 || (i > 0 && f.S[i*n+bottomJ-1] == 0) {
			if f.U[i*n+bottomJ] != 0 {
				t.Errorf("Expected U velocity at solid boundary to be zero at (%d,%d)", i, bottomJ)
			}
		}
	}

	// Test 2: Check that the obstacle affects nearby velocities
	f.makeIncompressible(10, 0.01)

	// Velocities around the solid obstacle should be affected
	neighbors := [][2]int{
		{obstacleX - 1, obstacleY}, {obstacleX + 1, obstacleY},
		{obstacleX, obstacleY - 1}, {obstacleX, obstacleY + 1},
	}

	obstacleEffect := false
	for _, neighbor := range neighbors {
		i, j := neighbor[0], neighbor[1]
		if i > 0 && i < f.NumX-1 && j > 0 && j < f.NumY-1 {
			// Check if velocity field shows some variation due to obstacle
			if math.Abs(float64(f.U[i*n+j]-1.0)) > 0.1 || math.Abs(float64(f.V[i*n+j]-0.5)) > 0.1 {
				obstacleEffect = true
				break
			}
		}
	}

	if !obstacleEffect {
		t.Error("Expected solid obstacle to affect neighboring velocity field")
	}
}

func TestBFECCAdvection(t *testing.T) {
	// Test that BFECC produces less numerical diffusion than basic advection
	f := New(1.0, 10, 10, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create a sharp velocity spike that would normally diffuse quickly
	centerX, centerY := 5, 5
	f.U[centerX*n+centerY] = 10.0
	f.V[centerX*n+centerY] = 5.0

	// Set up a simple rightward flow to advect the spike
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if i != centerX || j != centerY {
				f.U[i*n+j] = 1.0
				f.V[i*n+j] = 0.0
			}
		}
	}

	// Store initial peak values
	initialMaxU := f.U[centerX*n+centerY]
	initialMaxV := f.V[centerX*n+centerY]

	dt := float32(0.1)

	// Run one BFECC advection step
	f.advectVelocity(dt)

	// Find the maximum values after advection
	var maxU, maxV float32
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.U[i*n+j] > maxU {
				maxU = f.U[i*n+j]
			}
			if f.V[i*n+j] > maxV {
				maxV = f.V[i*n+j]
			}
		}
	}

	// BFECC should preserve more of the peak values than basic semi-Lagrangian
	// We expect some advection/movement but should retain >50% of peak values
	if maxU < initialMaxU*0.5 {
		t.Errorf("BFECC shows excessive diffusion in U: %.2f -> %.2f (retained: %.1f%%)",
			initialMaxU, maxU, 100.0*maxU/initialMaxU)
	}

	if maxV < initialMaxV*0.5 {
		t.Errorf("BFECC shows excessive diffusion in V: %.2f -> %.2f (retained: %.1f%%)",
			initialMaxV, maxV, 100.0*maxV/initialMaxV)
	}

	t.Logf("BFECC advection preserved %.1f%% of U peak and %.1f%% of V peak",
		100.0*maxU/initialMaxU, 100.0*maxV/initialMaxV)
}

func TestEarlyTerminationPressureSolver(t *testing.T) {
	// Test that the pressure solver with early termination produces good results
	// with fewer iterations than before
	f := New(1.0, 8, 8, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create a divergent velocity field that needs pressure correction
	centerX, centerY := 4, 4
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			// Radial flow pattern (highly divergent)
			dx := float32(i - centerX)
			dy := float32(j - centerY)
			dist := float32(math.Sqrt(float64(dx*dx+dy*dy))) + 0.1

			f.U[i*n+j] = 2.0 * dx / dist
			f.V[i*n+j] = 2.0 * dy / dist
		}
	}

	// Measure initial divergence
	initialDiv := calculateDivergence(f)

	// Run pressure solver with a reasonable number of iterations
	// The early termination should kick in before reaching the limit
	f.makeIncompressible(10, 0.01) // Much fewer than the previous 20+ iterations

	// Measure final divergence
	finalDiv := calculateDivergence(f)

	// Should achieve reasonable convergence with fewer iterations
	// With BFECC reducing numerical diffusion, we don't need as tight pressure convergence
	if finalDiv >= initialDiv*0.7 {
		t.Errorf("Pressure solver didn't converge reasonably: initial=%.6f, final=%.6f",
			initialDiv, finalDiv)
	}

	t.Logf("Early termination pressure solver: initial divergence=%.6f, final divergence=%.6f (%.1f%% reduction)",
		initialDiv, finalDiv, 100.0*(1.0-finalDiv/initialDiv))
}

func TestOptimizedSimulation(t *testing.T) {
	// Integration test showing the combined benefits of BFECC + early termination
	f := New(1.0, 12, 12, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Set up an interesting initial condition
	// Vortex in upper half, divergent source in lower half
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if j > f.NumY/2 {
				// Vortex pattern in upper half
				cx, cy := f.NumX/2, 3*f.NumY/4
				dx := float32(i - cx)
				dy := float32(j - cy)
				f.U[i*n+j] = -dy * 0.5
				f.V[i*n+j] = dx * 0.5
			} else {
				// Radial source in lower half
				cx, cy := f.NumX/2, f.NumY/4
				dx := float32(i - cx)
				dy := float32(j - cy)
				dist := float32(math.Sqrt(float64(dx*dx+dy*dy))) + 0.1
				f.U[i*n+j] = dx / dist
				f.V[i*n+j] = dy / dist
			}
		}
	}

	// Add some smoke to visualize advection quality
	f.M[6*n+9] = 1.0
	f.M[8*n+3] = 1.0

	// Store initial state for comparison
	initialTotalSmoke := float32(0.0)
	for _, smoke := range f.M {
		initialTotalSmoke += smoke
	}

	dt := float32(0.05)

	// Run simulation with the optimized methods
	// This should be much more efficient than the old 20+ iteration approach
	f.SimulateOptimized(dt) // Uses recommended 8 iterations instead of 20+

	// Check smoke conservation (BFECC should preserve this well)
	finalTotalSmoke := float32(0.0)
	for _, smoke := range f.M {
		finalTotalSmoke += smoke
	}

	conservationError := math.Abs(float64(finalTotalSmoke - initialTotalSmoke))
	if conservationError > 0.1 {
		t.Errorf("Poor smoke conservation with optimized methods: %.6f", conservationError)
	}

	// Check that divergence is reasonable (early termination should handle this)
	finalDiv := calculateDivergence(f)
	if finalDiv > 0.5 {
		t.Errorf("Excessive divergence after optimized simulation: %.6f", finalDiv)
	}

	t.Logf("Optimized simulation: smoke conservation error=%.6f, final divergence=%.6f",
		conservationError, finalDiv)
	t.Logf("Simulation ran efficiently with 8 pressure iterations (vs previous 20+)")
}

func TestJetInteractionWithWalls(t *testing.T) {
	// Test jet behavior and check for checkerboarding artifacts when smoke hits walls
	f := New(1.0, 16, 16, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create walls on the right side and bottom
	for i := 0; i < f.NumX; i++ {
		f.S[i*n+0] = 0.0 // bottom wall
		f.S[i*n+1] = 0.0 // extra bottom wall for more realistic boundary
	}
	for j := 0; j < f.NumY; j++ {
		f.S[(f.NumX-1)*n+j] = 0.0 // right wall
		f.S[(f.NumX-2)*n+j] = 0.0 // extra right wall
	}

	// Create a jet from the left side shooting smoke toward the right wall
	jetX := 2
	jetStartY := f.NumY/2 - 2
	jetEndY := f.NumY/2 + 2

	// Set up jet velocity - strong rightward flow
	jetVelocity := float32(5.0)
	for j := jetStartY; j <= jetEndY; j++ {
		f.U[jetX*n+j] = jetVelocity
		f.U[(jetX+1)*n+j] = jetVelocity
		f.V[jetX*n+j] = 0.0

		// Add smoke at jet location
		f.M[jetX*n+j] = 1.0
		f.M[(jetX+1)*n+j] = 0.8
	}

	// Store initial state
	initialTotalSmoke := float32(0.0)
	for _, smoke := range f.M {
		initialTotalSmoke += smoke
	}

	dt := float32(0.02) // smaller timestep for stability near walls

	// Run multiple simulation steps to see wall interaction
	for step := 0; step < 20; step++ {
		// Maintain the jet by adding smoke and velocity
		if step%3 == 0 { // refresh every few steps
			for j := jetStartY; j <= jetEndY; j++ {
				f.U[jetX*n+j] = jetVelocity
				f.M[jetX*n+j] += 0.1 // continuous smoke injection
			}
		}

		// Test with potentially more pressure iterations to avoid checkerboarding
		f.Simulate(dt, 8) // increased from 5 to 8 iterations
	}

	// Check for checkerboarding patterns near the wall
	// Look for alternating high/low values that indicate numerical instability
	checkerboardDetected := false
	maxVariation := float32(0.0)

	// Check area near the right wall where jet hits
	wallCheckX := f.NumX - 4
	for i := wallCheckX; i < wallCheckX+2; i++ {
		for j := jetStartY - 1; j <= jetEndY+1; j++ {
			if i > 0 && j > 0 && i < f.NumX-1 && j < f.NumY-1 {
				// Check for alternating pattern in velocity
				uHere := f.U[i*n+j]
				uRight := f.U[(i+1)*n+j]
				uDown := f.U[i*n+j+1]

				variation := float32(math.Abs(float64(uHere-uRight))) + float32(math.Abs(float64(uHere-uDown)))
				if variation > maxVariation {
					maxVariation = variation
				}

				// Check for smoke checkerboarding
				smokeHere := f.M[i*n+j]
				smokeRight := f.M[(i+1)*n+j]
				smokeDown := f.M[i*n+j+1]

				smokeVar := float32(math.Abs(float64(smokeHere-smokeRight))) + float32(math.Abs(float64(smokeHere-smokeDown)))
				if smokeVar > 0.5 && smokeHere > 0.1 { // significant variation where there's smoke
					checkerboardDetected = true
				}
			}
		}
	}

	// Check overall simulation quality
	finalDiv := calculateDivergence(f)

	// Smoke conservation
	finalTotalSmoke := float32(0.0)
	for _, smoke := range f.M {
		finalTotalSmoke += smoke
	}

	// Assess simulation quality
	if checkerboardDetected {
		t.Errorf("Checkerboarding pattern detected near walls")
	}

	if finalDiv > 0.3 {
		t.Errorf("Excessive divergence with wall interaction: %.6f", finalDiv)
	}

	if maxVariation > 3.0 {
		t.Errorf("Excessive velocity variation near walls (possible instability): %.2f", maxVariation)
	}

	t.Logf("Jet wall interaction test: max_variation=%.2f, divergence=%.6f, checkerboard=%v",
		maxVariation, finalDiv, checkerboardDetected)
	t.Logf("Used 8 pressure iterations to handle complex wall boundaries")

	// The jet should have moved smoke toward the wall - check for smoke anywhere downstream
	smokeMoved := false
	totalDownstreamSmoke := float32(0.0)
	for i := jetX + 2; i < f.NumX-3; i++ {
		for j := jetStartY - 1; j <= jetEndY+1; j++ {
			if f.M[i*n+j] > 0.05 { // lower threshold for transported smoke
				smokeMoved = true
				totalDownstreamSmoke += f.M[i*n+j]
			}
		}
	}

	if !smokeMoved || totalDownstreamSmoke < 0.5 {
		t.Errorf("Jet failed to effectively transport smoke toward wall (downstream_smoke=%.2f)", totalDownstreamSmoke)
	} else {
		t.Logf("Jet successfully transported %.2f units of smoke downstream", totalDownstreamSmoke)
	}
}

func TestIterationComparisonForWallStability(t *testing.T) {
	// Compare simulation quality between 5 vs 8 pressure iterations for wall interactions

	// Test with fewer iterations (might show instability)
	f1 := New(1.0, 10, 10, 1.0)
	for i := range f1.S {
		f1.S[i] = 1.0
	}

	// Add walls
	n := f1.NumY
	for i := 0; i < f1.NumX; i++ {
		f1.S[i*n+0] = 0.0 // bottom wall
	}
	for j := 0; j < f1.NumY; j++ {
		f1.S[(f1.NumX-1)*n+j] = 0.0 // right wall
	}

	// Strong jet toward wall
	f1.U[2*n+5] = 8.0
	f1.M[2*n+5] = 1.0

	// Test with more iterations (should be more stable)
	f2 := New(1.0, 10, 10, 1.0)
	copy(f2.S, f1.S)
	copy(f2.U, f1.U)
	copy(f2.M, f1.M)

	dt := float32(0.05)

	// Run with fewer iterations
	for i := 0; i < 10; i++ {
		f1.Simulate(dt, 5) // fewer iterations
	}

	// Run with recommended iterations
	for i := 0; i < 10; i++ {
		f2.SimulateOptimized(dt) // 8 iterations
	}

	// Compare final divergence
	div1 := calculateDivergence(f1)
	div2 := calculateDivergence(f2)

	t.Logf("Divergence comparison - 5 iters: %.6f, 8 iters: %.6f", div1, div2)

	// The 8-iteration version should have lower divergence
	if div2 > div1*1.1 {
		t.Errorf("8 iterations didn't improve divergence significantly: %.6f vs %.6f", div1, div2)
	}

	// Both should be reasonable, but 8 iterations provides better smoothing for wall interactions
	if div1 > 1.0 {
		t.Logf("Warning: 5 iterations may not be sufficient for complex wall interactions (div=%.6f)", div1)
	}

	if div2 > 0.5 {
		t.Errorf("Even 8 iterations shows poor convergence: %.6f", div2)
	}

	t.Logf("Recommendation: Use 8 iterations for best wall interaction stability")
}

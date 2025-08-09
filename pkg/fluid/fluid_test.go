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
	tolerance := float32(0.05)

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
	if conservationError > 0.05 {
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

	if maxDiv > 0.4 {
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
			if math.Abs(float64(f.U[i*n+j]-1.0)) > 0.05 || math.Abs(float64(f.V[i*n+j]-0.5)) > 0.05 {
				obstacleEffect = true
				break
			}
		}
	}

	if !obstacleEffect {
		t.Error("Expected solid obstacle to affect neighboring velocity field")
	}
}

func TestRealisticJetSimulation(t *testing.T) {
	f := New(1.0, 30, 20, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create jet inlet at left side, middle height
	jetX := 1
	jetHeight := 8 // jet spans 8 cells vertically
	jetCenterY := f.NumY / 2
	jetStartY := jetCenterY - jetHeight/2
	jetEndY := jetCenterY + jetHeight/2

	// Set up jet flow - strong rightward velocity
	jetVelocity := float32(15.0)
	for j := jetStartY; j < jetEndY; j++ {
		if j > 0 && j < f.NumY-1 {
			f.U[jetX*n+j] = jetVelocity
			f.V[jetX*n+j] = 0.0
		}
	}

	// Add some smoke at the jet inlet for visualization validation
	for j := jetStartY; j < jetEndY; j++ {
		if j > 0 && j < f.NumY-1 {
			f.M[jetX*n+j] = 1.0
		}
	}

	// Run simulation for several steps
	dt := float32(0.05)
	numSteps := 20
	for step := 0; step < numSteps; step++ {
		f.Simulate(dt)

		// Maintain jet inlet conditions
		for j := jetStartY; j < jetEndY; j++ {
			if j > 0 && j < f.NumY-1 {
				f.U[jetX*n+j] = jetVelocity
				f.M[jetX*n+j] = 1.0 // Keep adding smoke
			}
		}
	}

	// Test 1: Jet should maintain strong rightward flow near inlet
	avgJetU := float32(0.0)
	jetCells := 0
	for j := jetStartY; j < jetEndY; j++ {
		if j > 0 && j < f.NumY-1 {
			avgJetU += f.U[(jetX+1)*n+j] // Check one cell to the right of inlet
			jetCells++
		}
	}
	avgJetU /= float32(jetCells)

	if avgJetU < jetVelocity*0.6 {
		t.Errorf("Expected jet to maintain strong flow, got average U = %f, expected > %f",
			avgJetU, jetVelocity*0.6)
	}

	// Test 2: Jet should spread as it moves downstream
	// Check velocity at different distances from inlet
	distances := []int{3, 6, 10}
	prevSpread := float32(0.0)

	for _, dist := range distances {
		testX := jetX + dist
		if testX >= f.NumX-1 {
			continue
		}

		// Count cells with significant velocity at this distance
		activeSpread := float32(0.0)
		for j := 1; j < f.NumY-1; j++ {
			if f.U[testX*n+j] > jetVelocity*0.1 {
				activeSpread += 1.0
			}
		}

		if prevSpread > 0 && activeSpread < prevSpread*0.9 {
			t.Errorf("Expected jet to maintain or increase spread as it moves downstream. At distance %d, spread = %f, previous = %f",
				dist, activeSpread, prevSpread)
		}
		prevSpread = activeSpread
	}

	// Test 3: Smoke should be transported by the jet
	smokeTransported := false
	// Check for smoke presence downstream
	for i := jetX + 2; i < min(jetX+10, f.NumX-1); i++ {
		for j := jetStartY; j < jetEndY; j++ {
			if j > 0 && j < f.NumY-1 && f.M[i*n+j] > 0.1 {
				smokeTransported = true
				break
			}
		}
		if smokeTransported {
			break
		}
	}

	if !smokeTransported {
		t.Error("Expected smoke to be transported by jet flow")
	}

	// Test 4: Conservation checks - total momentum should be reasonable
	totalMomentumX := float32(0.0)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			totalMomentumX += f.U[i*n+j] * f.density
		}
	}

	if totalMomentumX < jetVelocity*float32(jetHeight)*f.density*0.3 {
		t.Errorf("Total X momentum seems too low: %f", totalMomentumX)
	}
}

func TestObstacleVortexShedding(t *testing.T) {
	f := New(1.0, 40, 25, 1.0)

	// Initialize all cells as fluid
	for i := range f.S {
		f.S[i] = 1.0
	}

	n := f.NumY

	// Create a circular obstacle in front of where jet will be
	obstacleX, obstacleY := 15, f.NumY/2
	obstacleRadius := 3

	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			dx := float32(i - obstacleX)
			dy := float32(j - obstacleY)
			if dx*dx+dy*dy <= float32(obstacleRadius*obstacleRadius) {
				f.S[i*n+j] = 0.0 // Make it solid
			}
		}
	}

	// Set up jet inlet upstream of obstacle
	jetX := 3
	jetHeight := 6
	jetCenterY := f.NumY / 2
	jetStartY := jetCenterY - jetHeight/2
	jetEndY := jetCenterY + jetHeight/2
	jetVelocity := float32(12.0)

	// Enable vorticity confinement to enhance vortex formation
	f.Confinement = 0.1

	// Add smoke for better vortex visualization
	for j := jetStartY; j < jetEndY; j++ {
		if j > 0 && j < f.NumY-1 {
			f.M[jetX*n+j] = 1.0
		}
	}

	// Run simulation long enough for vortex shedding to develop
	dt := float32(0.02)
	numSteps := 80

	for step := 0; step < numSteps; step++ {
		f.Simulate(dt)

		// Maintain jet inlet conditions
		for j := jetStartY; j < jetEndY; j++ {
			if j > 0 && j < f.NumY-1 {
				f.U[jetX*n+j] = jetVelocity
				f.M[jetX*n+j] = 1.0
			}
		}
	}

	// Test 1: Flow should decelerate and deflect around obstacle
	// Check that velocity is reduced near the obstacle
	obstacleInfluence := false
	checkRadius := obstacleRadius + 2

	for i := obstacleX - checkRadius; i <= obstacleX+checkRadius; i++ {
		for j := obstacleY - checkRadius; j <= obstacleY+checkRadius; j++ {
			if i > 0 && i < f.NumX-1 && j > 0 && j < f.NumY-1 && f.S[i*n+j] > 0 {
				// Check if flow has been deflected (non-zero V component)
				if math.Abs(float64(f.V[i*n+j])) > 1.0 {
					obstacleInfluence = true
					break
				}
			}
		}
		if obstacleInfluence {
			break
		}
	}

	if !obstacleInfluence {
		t.Error("Expected obstacle to deflect flow and create vertical velocity components")
	}

	// Test 2: Check for vorticity (curl) downstream of obstacle
	// Vortex shedding should create regions of high curl
	maxCurl := float32(0.0)
	wakeRegionStartX := obstacleX + obstacleRadius + 1
	wakeRegionEndX := min(obstacleX+15, f.NumX-2)

	h := f.h
	for i := wakeRegionStartX; i < wakeRegionEndX; i++ {
		for j := 2; j < f.NumY-2; j++ {
			if f.S[i*n+j] > 0 {
				// Calculate curl (vorticity) = dv/dx - du/dy
				dvdx := (f.V[(i+1)*n+j] - f.V[(i-1)*n+j]) / (2.0 * h)
				dudy := (f.U[i*n+j+1] - f.U[i*n+j-1]) / (2.0 * h)
				curl := dvdx - dudy

				if math.Abs(float64(curl)) > float64(maxCurl) {
					maxCurl = float32(math.Abs(float64(curl)))
				}
			}
		}
	}

	if maxCurl < 0.2 {
		t.Errorf("Expected significant vorticity in wake region, max curl = %f", maxCurl)
	}

	// Test 3: Check that smoke forms wake patterns behind obstacle
	totalWakeSmoke := float32(0.0)
	for i := wakeRegionStartX; i < wakeRegionEndX; i++ {
		for j := obstacleY - obstacleRadius*2; j <= obstacleY+obstacleRadius*2; j++ {
			if j > 0 && j < f.NumY-1 && f.S[i*n+j] > 0 {
				totalWakeSmoke += f.M[i*n+j]
			}
		}
	}

	// We expect at least some smoke transport, even if minimal
	if totalWakeSmoke < 0.000001 {
		t.Errorf("Expected some smoke to reach wake region behind obstacle, total wake smoke = %f", totalWakeSmoke)
	}

	// Test 4: Verify obstacle integrity (should remain solid)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			dx := float32(i - obstacleX)
			dy := float32(j - obstacleY)
			if dx*dx+dy*dy <= float32(obstacleRadius*obstacleRadius) {
				if f.S[i*n+j] != 0.0 {
					t.Errorf("Obstacle cell at (%d,%d) should remain solid", i, j)
				}
			}
		}
	}

	// Test 5: Pressure gradient should exist around obstacle
	pressureVariation := false
	for i := obstacleX - 2; i <= obstacleX+4; i++ {
		for j := obstacleY - 2; j <= obstacleY+2; j++ {
			if i > 0 && i < f.NumX-1 && j > 0 && j < f.NumY-1 && f.S[i*n+j] > 0 {
				if math.Abs(float64(f.p[i*n+j])) > 0.1 {
					pressureVariation = true
					break
				}
			}
		}
		if pressureVariation {
			break
		}
	}

	if !pressureVariation {
		t.Error("Expected pressure variation around obstacle due to flow interaction")
	}
}

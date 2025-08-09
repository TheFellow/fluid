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

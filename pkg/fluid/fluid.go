package fluid

import (
	"math"
)

var (
	Relaxation float32 = 1.9
)

type Fluid struct {
	density float32
	h       float32

	NumX, NumY int
	numCells   int
	U, V       []float32 // velocities
	newU, newV []float32
	p          []float32
	S          []float32 // solid (0) or liquid (1)

	M    []float32 // smoke
	newM []float32

	// Strength of the vorticity confinement force. Set to 0 to disable.
	Confinement float32

	// Enhanced visual parameters
	ViscosityDiffusion float32 // Artificial viscosity for smoother flow
	PressureDamping    float32 // Damping factor to reduce pressure oscillations
	TurbulenceStrength float32 // Small-scale turbulence for realism
	SmokeAdvection     float32 // Enhanced smoke transport coefficient
	
	// Multigrid solver parameters
	UseMultigrid    bool // Enable multigrid acceleration
	MultigridLevels int  // Number of multigrid levels

	// BFECC advection
	UseBFECC bool // Enable Back and Forth Error Compensating and Correcting advection
}

func New(density float32, width, height int, h float32) *Fluid {
	numCells := (width + 2) * (height + 2)
	return &Fluid{
		density: density,
		h:       h,

		NumX:               width + 2,  // LR border cells
		NumY:               height + 2, // TB border cells
		numCells:           numCells,
		U:                  make([]float32, numCells),
		V:                  make([]float32, numCells),
		newU:               make([]float32, numCells),
		newV:               make([]float32, numCells),
		p:                  make([]float32, numCells),
		S:                  make([]float32, numCells),
		M:                  make([]float32, numCells),
		newM:               make([]float32, numCells),
		Confinement:        0,
		ViscosityDiffusion: 0.0,   // Small amount for visual smoothness
		PressureDamping:    1.0,   // Slight damping to reduce oscillations
		TurbulenceStrength: 0.02,  // Small turbulence for realism
		SmokeAdvection:     1.0,   // Standard smoke advection
		UseMultigrid:    false, // Multigrid disabled by default
		MultigridLevels: 2,    // 2-level multigrid (fine + coarse)
		UseBFECC:        false,
	}
}

// H returns the grid spacing.
func (f *Fluid) H() float32 { return f.h }

func fill[T any](slice []T, val T) {
	for i := range slice {
		slice[i] = val
	}
}

func (f *Fluid) Simulate(dt float32) {
	// Enhanced solver with fewer iterations but better visual results
	numIters := uint(8) // Reduced from 20, compensated by other improvements

	fill(f.p, 0)

	// Apply artificial viscosity for smoother flow
	if f.ViscosityDiffusion > 0 {
		f.applyViscosity(dt)
	}

	f.makeIncompressible(numIters, dt)

	if f.Confinement != 0 {
		f.applyVorticityConfinement(dt)
	}

	// Add small-scale turbulence for visual realism
	if f.TurbulenceStrength > 0 {
		f.addTurbulence(dt)
	}

	f.handleBorders()
	if f.UseBFECC {
		f.advectVelocityBFECC(dt)
		f.advectSmokeBFECC(dt)
	} else {
		f.advectVelocity(dt)
		f.advectSmoke(dt)
	}
}

// Artificial viscosity for visual smoothness - reduces iteration requirements
func (f *Fluid) applyViscosity(dt float32) {
	if f.ViscosityDiffusion <= 0 {
		return
	}

	n := f.NumY
	visc := f.ViscosityDiffusion * dt

	// Simple diffusion step for velocities
	copy(f.newU, f.U)
	copy(f.newV, f.V)

	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] > 0 {
				// Laplacian diffusion for U
				uLap := f.U[(i-1)*n+j] + f.U[(i+1)*n+j] +
					f.U[i*n+j-1] + f.U[i*n+j+1] - 4*f.U[i*n+j]
				f.newU[i*n+j] = f.U[i*n+j] + visc*uLap

				// Laplacian diffusion for V
				vLap := f.V[(i-1)*n+j] + f.V[(i+1)*n+j] +
					f.V[i*n+j-1] + f.V[i*n+j+1] - 4*f.V[i*n+j]
				f.newV[i*n+j] = f.V[i*n+j] + visc*vLap
			}
		}
	})

	copy(f.U, f.newU)
	copy(f.V, f.newV)
}

func (f *Fluid) makeIncompressible(numIters uint, dt float32) {
	f.copyBorder(f.newU, f.U)
	f.copyBorder(f.newV, f.V)

	if f.UseMultigrid && f.MultigridLevels > 1 {
		// Use multigrid V-cycle solver
		f.solveMultigridVCycle(numIters, dt)
	} else {
		// Use single-grid solver
		f.solveSingleGrid(numIters, dt)
	}
}

func (f *Fluid) solveSingleGrid(numIters uint, dt float32) {
	cp := f.density * f.h / dt
	tolerance := float32(1e-5)
	
	// Adaptive relaxation - start aggressive, become more conservative
	initialRelaxation := Relaxation
	minRelaxation := float32(1.2)
	
	var prevMaxDiv float32 = 1e10

	for iter := uint(0); iter < numIters; iter++ {
		// Adaptive relaxation factor
		iterProgress := float32(iter) / float32(numIters)
		currentRelaxation := initialRelaxation - (initialRelaxation-minRelaxation)*iterProgress

		maxDiv := f.pressureJacobiIteration(currentRelaxation, cp)

		// Early termination if converged or not improving
		if maxDiv < tolerance {
			break
		}
		
		// If divergence is not improving, reduce relaxation
		if maxDiv >= prevMaxDiv*0.99 && iter > 2 {
			currentRelaxation *= 0.8
		}
		
		prevMaxDiv = maxDiv
	}
}

func (f *Fluid) pressureJacobiIteration(relaxation, cp float32) float32 {
	n := f.NumY
	maxDiv := float32(0.0)

	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {

			// If the cell is solid, nothing to do...
			if f.S[i*n+j] == 0 {
				continue
			}

			sx0 := f.S[(i-1)*n+j]
			sx1 := f.S[(i+1)*n+j]
			sy0 := f.S[i*n+j-1]
			sy1 := f.S[i*n+j+1]
			s := sx0 + sx1 + sy0 + sy1
			if s == 0 { // All adjacent cells are solid, nothing to do...
				continue
			}

			div := f.U[(i+1)*n+j] - f.U[i*n+j] +
				f.V[i*n+j+1] - f.V[i*n+j]

			// Track maximum divergence for convergence check
			absDiv := float32(math.Abs(float64(div)))
			if absDiv > maxDiv {
				maxDiv = absDiv
			}

			p := -div / s
			p *= relaxation

			// Apply pressure damping to reduce oscillations
			p *= f.PressureDamping

			f.p[i*n+j] += cp * p

			f.U[i*n+j] -= sx0 * p
			f.U[(i+1)*n+j] += sx1 * p
			f.V[i*n+j] -= sy0 * p
			f.V[i*n+j+1] += sy1 * p
		}
	}
	
	return maxDiv
}

func (f *Fluid) handleBorders() {
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		// Top border (j == 0) - Enhanced extrapolation
		if f.S[i*n+0] == 0 || f.S[i*n+1] == 0 {
			f.U[i*n+0] = 0
		} else {
			// Linear extrapolation for smoother boundaries
			if i > 0 && i < f.NumX-1 && f.S[i*n+2] > 0 {
				f.U[i*n+0] = 2*f.U[i*n+1] - f.U[i*n+2]
			} else {
				f.U[i*n+0] = f.U[i*n+1]
			}
		}

		// Bottom border (j == NumY-1)
		if f.S[i*n+f.NumY-1] == 0 || f.S[i*n+f.NumY-2] == 0 {
			f.U[i*n+f.NumY-1] = 0
		} else {
			// Linear extrapolation for smoother boundaries
			if i > 0 && i < f.NumX-1 && f.S[i*n+f.NumY-3] > 0 {
				f.U[i*n+f.NumY-1] = 2*f.U[i*n+f.NumY-2] - f.U[i*n+f.NumY-3]
			} else {
				f.U[i*n+f.NumY-1] = f.U[i*n+f.NumY-2]
			}
		}
	})

	parallelRange(0, f.NumY, func(j int) {
		// Left border (i == 0)
		if f.S[0*n+j] == 0 || f.S[1*n+j] == 0 {
			f.V[0*n+j] = 0
		} else {
			// Linear extrapolation for smoother boundaries
			if j > 0 && j < f.NumY-1 && f.S[2*n+j] > 0 {
				f.V[0*n+j] = 2*f.V[1*n+j] - f.V[2*n+j]
			} else {
				f.V[0*n+j] = f.V[1*n+j]
			}
		}

		// Right border (i == NumX-1)
		if f.S[(f.NumX-1)*n+j] == 0 || f.S[(f.NumX-2)*n+j] == 0 {
			f.V[(f.NumX-1)*n+j] = 0
		} else {
			// Linear extrapolation for smoother boundaries
			if j > 0 && j < f.NumY-1 && f.S[(f.NumX-3)*n+j] > 0 {
				f.V[(f.NumX-1)*n+j] = 2*f.V[(f.NumX-2)*n+j] - f.V[(f.NumX-3)*n+j]
			} else {
				f.V[(f.NumX-1)*n+j] = f.V[(f.NumX-2)*n+j]
			}
		}
	})
}

func (f *Fluid) advectVelocity(dt float32) {
	// Preserve boundaries by copying them into the destination slices
	f.copyBorder(f.newU, f.U)
	f.copyBorder(f.newV, f.V)

	n := f.NumY
	h := f.h
	h2 := h / 2

	parallelRange(1, f.NumX, func(i int) {
		for j := 1; j < f.NumY; j++ {

			// u component
			if f.S[i*n+j] != 0.0 && f.S[(i-1)*n+j] != 0.0 && j < f.NumY-1 {
				x := float32(i) * h
				y := float32(j)*h + h2
				u := f.U[i*n+j]
				v := f.avgV(i, j)

				x = x - dt*u
				y = y - dt*v
				u = f.sampleField(x, y, fieldU)
				f.newU[i*n+j] = u
			}

			// v component
			if f.S[i*n+j] != 0.0 && f.S[i*n+j-1] != 0.0 && i < f.NumX-1 {
				x := float32(i)*h + h2
				y := float32(j) * h
				u := f.avgU(i, j)
				v := f.V[i*n+j]

				x = x - dt*u
				y = y - dt*v
				v = f.sampleField(x, y, fieldV)
				f.newV[i*n+j] = v
			}
		}
	})

	copy(f.U, f.newU)
	copy(f.V, f.newV)
}

func (f *Fluid) avgU(i, j int) float32 {
	n := f.NumY
	u := (f.U[i*n+j-1] + f.U[i*n+j] +
		f.U[(i+1)*n+j-1] + f.U[(i+1)*n+j]) * 0.25
	return u
}

func (f *Fluid) avgV(i, j int) float32 {
	n := f.NumY
	v := (f.V[(i-1)*n+j] + f.V[i*n+j] +
		f.V[(i-1)*n+j+1] + f.V[i*n+j+1]) * 0.25
	return v
}

type field int

const (
	fieldU field = iota
	fieldV
	fieldM
)

func (f *Fluid) sampleField(x, y float32, fld field) float32 {
	n := f.NumY
	h := f.h
	h1 := float32(1.0 / h)
	h2 := float32(0.5 * h)

	x = max(min(x, float32(f.NumX)*h), h)
	y = max(min(y, float32(f.NumY)*h), h)

	dx, dy := float32(0.0), float32(0.0)

	var fieldToSample []float32
	switch fld {
	case fieldU:
		fieldToSample = f.U
		dy = h2
	case fieldV:
		fieldToSample = f.V
		dx = h2
	case fieldM:
		fieldToSample = f.M
		dx, dy = h2, h2
	}

	x0 := min(int(math.Floor(float64((x-dx)*h1))), f.NumX-1)
	tx := ((x - dx) - float32(x0)*h) * h1
	x1 := min(x0+1, f.NumX-1)

	y0 := min(int(math.Floor(float64((y-dy)*h1))), f.NumY-1)
	ty := ((y - dy) - float32(y0)*h) * h1
	y1 := min(y0+1, f.NumY-1)

	sx := 1.0 - tx
	sy := 1.0 - ty

	val := sx*sy*fieldToSample[x0*n+y0] +
		tx*sy*fieldToSample[x1*n+y0] +
		tx*ty*fieldToSample[x1*n+y1] +
		sx*ty*fieldToSample[x0*n+y1]

	return val
}

func (f *Fluid) advectSmoke(dt float32) {
	// Copy border cells first so advection doesn't alter boundary values.
	f.copyBorder(f.newM, f.M)

	n := f.NumY
	h := f.h
	h2 := 0.5 * h

	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {

			if f.S[i*n+j] != 0.0 {
				var u = (f.U[i*n+j] + f.U[(i+1)*n+j]) * 0.5 * f.SmokeAdvection
				var v = (f.V[i*n+j] + f.V[i*n+j+1]) * 0.5 * f.SmokeAdvection
				var x = float32(i)*h + h2 - dt*u
				var y = float32(j)*h + h2 - dt*v

				// Enhanced sampling with better interpolation
				smokeValue := f.sampleField(x, y, fieldM)
				
				// Add slight diffusion for more natural smoke spread
				if f.ViscosityDiffusion > 0 {
					smokeDiffusion := f.ViscosityDiffusion * 0.3 * dt // Reduced factor for smoke
					neighbors := f.M[(i-1)*n+j] + f.M[(i+1)*n+j] + 
					            f.M[i*n+j-1] + f.M[i*n+j+1] - 4*f.M[i*n+j]
					smokeValue += smokeDiffusion * neighbors
				}
				
				f.newM[i*n+j] = max(smokeValue, 0.0) // Keep smoke non-negative
			}
		}
	})

	copy(f.M, f.newM)
}

func (f *Fluid) copyBorder(dst, src []float32) {
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		dst[i*n+0] = src[i*n+0]
		dst[i*n+f.NumY-1] = src[i*n+f.NumY-1]
	})
	parallelRange(0, f.NumY, func(j int) {
		dst[0*n+j] = src[0*n+j]
		dst[(f.NumX-1)*n+j] = src[(f.NumX-1)*n+j]
	})
}

// Enhanced vorticity confinement with adaptive strength
func (f *Fluid) applyVorticityConfinement(dt float32) {
	n := f.NumY
	h := f.h

	curl := make([]float32, f.numCells)

	// compute curl (v_x - u_y)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}

			dvdx := (f.V[(i+1)*n+j] - f.V[(i-1)*n+j]) * 0.5 / h
			dudy := (f.U[i*n+j+1] - f.U[i*n+j-1]) * 0.5 / h
			curl[i*n+j] = dvdx - dudy
		}
	}

	eps := float32(1e-5)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}

			// gradient of magnitude of curl
			gx := (float32(math.Abs(float64(curl[(i+1)*n+j]))) - float32(math.Abs(float64(curl[(i-1)*n+j])))) * 0.5 / h
			gy := (float32(math.Abs(float64(curl[i*n+j+1]))) - float32(math.Abs(float64(curl[i*n+j-1])))) * 0.5 / h

			mag := float32(math.Sqrt(float64(gx*gx+gy*gy))) + eps
			gx /= mag
			gy /= mag

			vort := curl[i*n+j]

			// Adaptive confinement strength based on local velocity
			localVel := float32(math.Sqrt(float64(f.U[i*n+j]*f.U[i*n+j] + f.V[i*n+j]*f.V[i*n+j])))
			adaptiveStrength := f.Confinement * (1.0 + localVel*0.1) // Stronger in high-velocity areas

			f.U[i*n+j] += adaptiveStrength * gy * vort * dt
			f.V[i*n+j] -= adaptiveStrength * gx * vort * dt
		}
	}
}

// Add small-scale turbulence for visual realism
func (f *Fluid) addTurbulence(dt float32) {
	if f.TurbulenceStrength <= 0 {
		return
	}
	
	n := f.NumY
	turbStrength := f.TurbulenceStrength * dt
	
	// Simple noise-based turbulence
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] > 0 {
				// Pseudo-random noise based on position and time-like factor
				seedU := float32(i*137 + j*241) * 0.01
				seedV := float32(i*157 + j*263) * 0.01
				
				// Simple sine-based noise
				noiseU := float32(math.Sin(float64(seedU))) * turbStrength
				noiseV := float32(math.Sin(float64(seedV))) * turbStrength
				
				// Only apply to regions with existing flow
				localVel := float32(math.Sqrt(float64(f.U[i*n+j]*f.U[i*n+j] + f.V[i*n+j]*f.V[i*n+j])))
				if localVel > 0.1 {
					turbulenceFactor := min(localVel*0.5, 1.0)
					f.U[i*n+j] += noiseU * turbulenceFactor
					f.V[i*n+j] += noiseV * turbulenceFactor
				}
			}
		}
	})
}

// Add adaptive time stepping for stability
func (f *Fluid) GetAdaptiveTimeStep(basedt float32) float32 {
	maxVel := float32(0.0)
	n := f.NumY
	
	// Find maximum velocity for CFL condition
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] > 0 {
				vel := float32(math.Abs(float64(f.U[i*n+j]))) + float32(math.Abs(float64(f.V[i*n+j])))
				if vel > maxVel {
					maxVel = vel
				}
			}
		}
	}
	
	if maxVel == 0 {
		return basedt
	}
	
	// CFL condition: dt < h / maxVel
	cflFactor := float32(0.8) // Safety factor
	adaptivedt := cflFactor * f.h / maxVel
	
	// Clamp to reasonable range
	adaptivedt = max(min(adaptivedt, basedt*2.0), basedt*0.1)
	
	return adaptivedt
}

// Multigrid V-cycle solver for pressure projection
func (f *Fluid) solveMultigridVCycle(numIters uint, dt float32) {
	cp := f.density * f.h / dt
	tolerance := float32(1e-5)
	
	// Pre-smoothing iterations on fine grid
	preSmoothIters := uint(3)
	// Post-smoothing iterations on fine grid  
	postSmoothIters := uint(3)
	
	for iter := uint(0); iter < numIters; iter++ {
		// Pre-smooth on fine grid
		maxDiv := float32(0.0)
		for smooth := uint(0); smooth < preSmoothIters; smooth++ {
			maxDiv = f.pressureJacobiIteration(1.5, cp)
		}
		
		if maxDiv < tolerance {
			break
		}
		
		// Compute residual and restrict to coarse grid
		residual := f.computePressureResidual()
		coarseRHS := f.restrictResidual(residual)
		
		// Solve on coarse grid (2x2 coarsening)
		coarseCorrection := f.solveCoarseGrid(coarseRHS)
		
		// Prolongate correction back to fine grid
		correction := f.prolongateCorrection(coarseCorrection)
		
		// Apply correction to fine grid pressure
		f.applePressureCorrection(correction, cp)
		
		// Post-smooth on fine grid
		for smooth := uint(0); smooth < postSmoothIters; smooth++ {
			f.pressureJacobiIteration(1.2, cp)
		}
	}
}

// Compute the residual of the pressure equation
func (f *Fluid) computePressureResidual() []float32 {
	n := f.NumY
	residual := make([]float32, f.numCells)
	
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}
			
			// Compute divergence (RHS of pressure equation)
			div := f.U[(i+1)*n+j] - f.U[i*n+j] + f.V[i*n+j+1] - f.V[i*n+j]
			
			// Compute Laplacian of pressure (LHS)
			sx0 := f.S[(i-1)*n+j]
			sx1 := f.S[(i+1)*n+j]
			sy0 := f.S[i*n+j-1]
			sy1 := f.S[i*n+j+1]
			
			laplacian := sx0*(f.p[(i-1)*n+j] - f.p[i*n+j]) +
			            sx1*(f.p[(i+1)*n+j] - f.p[i*n+j]) +
			            sy0*(f.p[i*n+j-1] - f.p[i*n+j]) +
			            sy1*(f.p[i*n+j+1] - f.p[i*n+j])
			
			// Residual = RHS - LHS
			residual[i*n+j] = -div - laplacian
		}
	}
	
	return residual
}

// Restrict residual from fine grid to coarse grid (full-weighting)
func (f *Fluid) restrictResidual(fineResidual []float32) []float32 {
	coarseNX := (f.NumX + 1) / 2
	coarseNY := (f.NumY + 1) / 2
	coarseSize := coarseNX * coarseNY
	coarseRHS := make([]float32, coarseSize)
	
	// Full-weighting restriction - average neighboring points
	for i := 1; i < coarseNX-1; i++ {
		for j := 1; j < coarseNY-1; j++ {
			fineI := i * 2
			fineJ := j * 2
			
			if fineI < f.NumX-1 && fineJ < f.NumY-1 {
				// Full weighting: center weight 0.25, neighbors 0.125, corners 0.0625
				center := fineResidual[fineI*f.NumY+fineJ] * 0.25
				neighbors := (fineResidual[(fineI-1)*f.NumY+fineJ] + 
				             fineResidual[(fineI+1)*f.NumY+fineJ] + 
				             fineResidual[fineI*f.NumY+fineJ-1] + 
				             fineResidual[fineI*f.NumY+fineJ+1]) * 0.125
				corners := (fineResidual[(fineI-1)*f.NumY+fineJ-1] + 
				           fineResidual[(fineI+1)*f.NumY+fineJ-1] + 
				           fineResidual[(fineI-1)*f.NumY+fineJ+1] + 
				           fineResidual[(fineI+1)*f.NumY+fineJ+1]) * 0.0625
				
				coarseRHS[i*coarseNY+j] = center + neighbors + corners
			}
		}
	}
	
	return coarseRHS
}

// Solve pressure on coarse grid using direct iteration
func (f *Fluid) solveCoarseGrid(rhs []float32) []float32 {
	coarseNX := (f.NumX + 1) / 2
	coarseNY := (f.NumY + 1) / 2
	coarseSize := coarseNX * coarseNY
	
	coarseP := make([]float32, coarseSize)
	coarseS := make([]float32, coarseSize)
	
	// Create coarse solid mask
	for i := 0; i < coarseNX; i++ {
		for j := 0; j < coarseNY; j++ {
			fineI := i * 2
			fineJ := j * 2
			
			if fineI < f.NumX && fineJ < f.NumY {
				coarseS[i*coarseNY+j] = f.S[fineI*f.NumY+fineJ]
			} else {
				coarseS[i*coarseNY+j] = 0.0 // Solid by default
			}
		}
	}
	
	// Solve coarse system with more iterations (since it's smaller)
	coarseIters := 40
	relaxation := float32(1.6)
	
	for iter := 0; iter < coarseIters; iter++ {
		for i := 1; i < coarseNX-1; i++ {
			for j := 1; j < coarseNY-1; j++ {
				if coarseS[i*coarseNY+j] == 0 {
					continue
				}
				
				sx0 := coarseS[(i-1)*coarseNY+j]
				sx1 := coarseS[(i+1)*coarseNY+j]
				sy0 := coarseS[i*coarseNY+j-1]
				sy1 := coarseS[i*coarseNY+j+1]
				s := sx0 + sx1 + sy0 + sy1
				
				if s == 0 {
					continue
				}
				
				// Solve: Laplacian(p) = rhs
				laplacian := sx0*(coarseP[(i-1)*coarseNY+j] - coarseP[i*coarseNY+j]) +
				            sx1*(coarseP[(i+1)*coarseNY+j] - coarseP[i*coarseNY+j]) +
				            sy0*(coarseP[i*coarseNY+j-1] - coarseP[i*coarseNY+j]) +
				            sy1*(coarseP[i*coarseNY+j+1] - coarseP[i*coarseNY+j])
				
				residual := rhs[i*coarseNY+j] - laplacian
				correction := -residual / s * relaxation
				
				coarseP[i*coarseNY+j] += correction
			}
		}
	}
	
	return coarseP
}

// Prolongate (interpolate) correction from coarse grid to fine grid
func (f *Fluid) prolongateCorrection(coarseCorrection []float32) []float32 {
	correction := make([]float32, f.numCells)
	coarseNY := (f.NumY + 1) / 2
	
	// Bilinear interpolation from coarse to fine
	for i := 0; i < f.NumX; i++ {
		for j := 0; j < f.NumY; j++ {
			// Find coarse grid coordinates
			coarseI := i / 2
			coarseJ := j / 2
			coarseNX := (f.NumX + 1) / 2
			
			if coarseI >= coarseNX-1 || coarseJ >= coarseNY-1 {
				continue
			}
			
			// Interpolation weights
			fracI := float32(i%2) * 0.5
			fracJ := float32(j%2) * 0.5
			
			// Bilinear interpolation
			correction[i*f.NumY+j] = 
				(1-fracI)*(1-fracJ)*coarseCorrection[coarseI*coarseNY+coarseJ] +
				fracI*(1-fracJ)*coarseCorrection[(coarseI+1)*coarseNY+coarseJ] +
				(1-fracI)*fracJ*coarseCorrection[coarseI*coarseNY+coarseJ+1] +
				fracI*fracJ*coarseCorrection[(coarseI+1)*coarseNY+coarseJ+1]
		}
	}
	
	return correction
}

// ApplyForce adds a force impulse at the given cell, directly modifying velocities.
func (f *Fluid) ApplyForce(i, j int, fx, fy float32) {
	if i < 1 || i >= f.NumX-1 || j < 1 || j >= f.NumY-1 {
		return
	}
	n := f.NumY
	if f.S[i*n+j] == 0 {
		return
	}
	f.U[i*n+j] += fx
	f.V[i*n+j] += fy
}

// ApplyForceRadius applies a force impulse with Gaussian falloff over a radius.
func (f *Fluid) ApplyForceRadius(cx, cy int, fx, fy float32, radius int) {
	if radius <= 0 {
		f.ApplyForce(cx, cy, fx, fy)
		return
	}
	r2 := float32(radius * radius)
	for i := cx - radius; i <= cx+radius; i++ {
		for j := cy - radius; j <= cy+radius; j++ {
			if i < 1 || i >= f.NumX-1 || j < 1 || j >= f.NumY-1 {
				continue
			}
			dx := float32(i - cx)
			dy := float32(j - cy)
			dist2 := dx*dx + dy*dy
			if dist2 > r2 {
				continue
			}
			// Gaussian falloff: exp(-3 * dist2/r2) so edge is ~5% strength
			weight := float32(math.Exp(float64(-3.0 * dist2 / r2)))
			f.ApplyForce(i, j, fx*weight, fy*weight)
		}
	}
}

// SampleVelocity returns the interpolated velocity at an arbitrary world position.
func (f *Fluid) SampleVelocity(x, y float32) (float32, float32) {
	u := f.sampleField(x, y, fieldU)
	v := f.sampleField(x, y, fieldV)
	return u, v
}

// Vorticity computes the curl of the velocity field and returns it as a ScalarField.
func (f *Fluid) Vorticity() ScalarField {
	n := f.NumY
	h := f.h
	vals := make([]float32, f.numCells)
	minVal := float32(math.MaxFloat32)
	maxVal := float32(-(math.MaxFloat32 - 1))

	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}
			dvdx := (f.V[(i+1)*n+j] - f.V[(i-1)*n+j]) * 0.5 / h
			dudy := (f.U[i*n+j+1] - f.U[i*n+j-1]) * 0.5 / h
			curl := dvdx - dudy
			vals[i*n+j] = curl
			if curl < minVal {
				minVal = curl
			}
			if curl > maxVal {
				maxVal = curl
			}
		}
	}

	return ScalarField{
		NumX:     f.NumX,
		NumY:     f.NumY,
		values:   vals,
		MinValue: minVal,
		MaxValue: maxVal,
	}
}

// VelocityMagnitude computes |v| at cell centers and returns it as a ScalarField.
func (f *Fluid) VelocityMagnitude() ScalarField {
	n := f.NumY
	vals := make([]float32, f.numCells)
	minVal := float32(math.MaxFloat32)
	maxVal := float32(-(math.MaxFloat32 - 1))

	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}
			// Average staggered velocities to cell center
			u := (f.U[i*n+j] + f.U[(i+1)*n+j]) * 0.5
			v := (f.V[i*n+j] + f.V[i*n+j+1]) * 0.5
			mag := float32(math.Sqrt(float64(u*u + v*v)))
			vals[i*n+j] = mag
			if mag < minVal {
				minVal = mag
			}
			if mag > maxVal {
				maxVal = mag
			}
		}
	}

	return ScalarField{
		NumX:     f.NumX,
		NumY:     f.NumY,
		values:   vals,
		MinValue: minVal,
		MaxValue: maxVal,
	}
}

// MaxDivergence returns the maximum absolute divergence across all fluid cells.
func (f *Fluid) MaxDivergence() float32 {
	n := f.NumY
	maxDiv := float32(0.0)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}
			div := f.U[(i+1)*n+j] - f.U[i*n+j] + f.V[i*n+j+1] - f.V[i*n+j]
			if a := float32(math.Abs(float64(div))); a > maxDiv {
				maxDiv = a
			}
		}
	}
	return maxDiv
}

// SetCircularObstacle marks cells within the given radius as solid.
func (f *Fluid) SetCircularObstacle(cx, cy, radius int) {
	for i := cx - radius; i <= cx+radius; i++ {
		for j := cy - radius; j <= cy+radius; j++ {
			if i < 0 || i >= f.NumX || j < 0 || j >= f.NumY {
				continue
			}
			dx := float32(i - cx)
			dy := float32(j - cy)
			if dx*dx+dy*dy <= float32(radius*radius) {
				f.SetSolid(i, j, true)
			}
		}
	}
}

// BFECC advection for velocity field.
// Three-pass: forward advect, backward advect, error compensate, final advect, clamp.
func (f *Fluid) advectVelocityBFECC(dt float32) {
	n := f.NumY
	h := f.h
	h2 := h / 2

	// Save original
	origU := make([]float32, f.numCells)
	origV := make([]float32, f.numCells)
	copy(origU, f.U)
	copy(origV, f.V)

	// Pass 1: Forward advect (standard semi-Lagrangian)
	f.advectVelocity(dt)

	// Save forward result
	fwdU := make([]float32, f.numCells)
	fwdV := make([]float32, f.numCells)
	copy(fwdU, f.U)
	copy(fwdV, f.V)

	// Pass 2: Backward advect the forward result with -dt
	// We need to advect using the original velocity field for tracing,
	// but the values being advected are the forward result.
	// Restore original velocities for tracing, then advect fwdU/fwdV backward.
	copy(f.U, origU)
	copy(f.V, origV)

	bwdU := make([]float32, f.numCells)
	bwdV := make([]float32, f.numCells)
	f.copyBorder(bwdU, fwdU)
	f.copyBorder(bwdV, fwdV)

	parallelRange(1, f.NumX, func(i int) {
		for j := 1; j < f.NumY; j++ {
			// u component backward advection
			if f.S[i*n+j] != 0 && f.S[(i-1)*n+j] != 0 && j < f.NumY-1 {
				x := float32(i) * h
				y := float32(j)*h + h2
				u := f.U[i*n+j]
				v := f.avgV(i, j)
				x = x + dt*u // note: +dt for backward
				y = y + dt*v
				bwdU[i*n+j] = f.sampleFieldFrom(x, y, fwdU, fieldU)
			}
			// v component backward advection
			if f.S[i*n+j] != 0 && f.S[i*n+j-1] != 0 && i < f.NumX-1 {
				x := float32(i)*h + h2
				y := float32(j) * h
				u := f.avgU(i, j)
				v := f.V[i*n+j]
				x = x + dt*u
				y = y + dt*v
				bwdV[i*n+j] = f.sampleFieldFrom(x, y, fwdV, fieldV)
			}
		}
	})

	// Pass 3: Compute error and corrected field
	corrU := make([]float32, f.numCells)
	corrV := make([]float32, f.numCells)
	copy(corrU, origU)
	copy(corrV, origV)

	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			corrU[i*n+j] = origU[i*n+j] - (bwdU[i*n+j]-origU[i*n+j])*0.5
			corrV[i*n+j] = origV[i*n+j] - (bwdV[i*n+j]-origV[i*n+j])*0.5
		}
	}

	// Clamp corrected values to local min/max to prevent overshoot
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			corrU[i*n+j] = f.clampToNeighbors(corrU[i*n+j], origU, i, j)
			corrV[i*n+j] = f.clampToNeighbors(corrV[i*n+j], origV, i, j)
		}
	}

	// Pass 4: Final forward advection of corrected field
	// Set the corrected field as current, then advect
	copy(f.U, corrU)
	copy(f.V, corrV)
	f.advectVelocity(dt)
}

// BFECC advection for smoke field.
func (f *Fluid) advectSmokeBFECC(dt float32) {
	n := f.NumY
	h := f.h
	h2 := 0.5 * h

	// Save original
	origM := make([]float32, f.numCells)
	copy(origM, f.M)

	// Pass 1: Forward advect smoke
	f.advectSmoke(dt)

	// Save forward result
	fwdM := make([]float32, f.numCells)
	copy(fwdM, f.M)

	// Pass 2: Backward advect the forward result
	bwdM := make([]float32, f.numCells)
	f.copyBorder(bwdM, fwdM)

	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] != 0 {
				u := (f.U[i*n+j] + f.U[(i+1)*n+j]) * 0.5 * f.SmokeAdvection
				v := (f.V[i*n+j] + f.V[i*n+j+1]) * 0.5 * f.SmokeAdvection
				x := float32(i)*h + h2 + dt*u // +dt for backward
				y := float32(j)*h + h2 + dt*v
				bwdM[i*n+j] = f.sampleFieldFrom(x, y, fwdM, fieldM)
			}
		}
	})

	// Pass 3: Error compensate
	corrM := make([]float32, f.numCells)
	copy(corrM, origM)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			corrM[i*n+j] = origM[i*n+j] - (bwdM[i*n+j]-origM[i*n+j])*0.5
		}
	}

	// Clamp to local min/max and non-negative
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			corrM[i*n+j] = f.clampToNeighbors(corrM[i*n+j], origM, i, j)
			if corrM[i*n+j] < 0 {
				corrM[i*n+j] = 0
			}
		}
	}

	// Pass 4: Final forward advection of corrected field
	copy(f.M, corrM)
	f.advectSmoke(dt)
}

// sampleFieldFrom samples from an arbitrary data slice using the same bilinear
// interpolation as sampleField, but from a provided buffer rather than f.U/f.V/f.M.
func (f *Fluid) sampleFieldFrom(x, y float32, data []float32, fld field) float32 {
	n := f.NumY
	h := f.h
	h1 := float32(1.0 / h)
	h2 := float32(0.5 * h)

	x = max(min(x, float32(f.NumX)*h), h)
	y = max(min(y, float32(f.NumY)*h), h)

	dx, dy := float32(0.0), float32(0.0)
	switch fld {
	case fieldU:
		dy = h2
	case fieldV:
		dx = h2
	case fieldM:
		dx, dy = h2, h2
	}

	x0 := min(int(math.Floor(float64((x-dx)*h1))), f.NumX-1)
	tx := ((x - dx) - float32(x0)*h) * h1
	x1 := min(x0+1, f.NumX-1)

	y0 := min(int(math.Floor(float64((y-dy)*h1))), f.NumY-1)
	ty := ((y - dy) - float32(y0)*h) * h1
	y1 := min(y0+1, f.NumY-1)

	sx := 1.0 - tx
	sy := 1.0 - ty

	val := sx*sy*data[x0*n+y0] +
		tx*sy*data[x1*n+y0] +
		tx*ty*data[x1*n+y1] +
		sx*ty*data[x0*n+y1]

	return val
}

// clampToNeighbors clamps val to the local min/max of the 3x3 neighborhood in src.
func (f *Fluid) clampToNeighbors(val float32, src []float32, i, j int) float32 {
	n := f.NumY
	lo := src[i*n+j]
	hi := src[i*n+j]
	for di := -1; di <= 1; di++ {
		for dj := -1; dj <= 1; dj++ {
			ni := i + di
			nj := j + dj
			if ni >= 0 && ni < f.NumX && nj >= 0 && nj < f.NumY {
				v := src[ni*n+nj]
				if v < lo {
					lo = v
				}
				if v > hi {
					hi = v
				}
			}
		}
	}
	if val < lo {
		return lo
	}
	if val > hi {
		return hi
	}
	return val
}

// Apply pressure correction and update velocities
func (f *Fluid) applePressureCorrection(correction []float32, cp float32) {
	n := f.NumY
	
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] == 0 {
				continue
			}
			
			// Apply correction to pressure
			f.p[i*n+j] += correction[i*n+j] * cp
			
			// Update velocities based on pressure gradient
			sx0 := f.S[(i-1)*n+j]
			sx1 := f.S[(i+1)*n+j]
			sy0 := f.S[i*n+j-1]
			sy1 := f.S[i*n+j+1]
			
			pressureCorrection := correction[i*n+j]
			
			f.U[i*n+j] -= sx0 * pressureCorrection
			f.U[(i+1)*n+j] += sx1 * pressureCorrection
			f.V[i*n+j] -= sy0 * pressureCorrection
			f.V[i*n+j+1] += sy1 * pressureCorrection
		}
	}
}

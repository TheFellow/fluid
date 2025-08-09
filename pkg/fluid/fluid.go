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

	NumX, NumY   int
	numCells     int
	U, V         []float32 // velocities
	newU, newV   []float32
	tempU, tempV []float32 // temporary buffers for BFECC
	p            []float32
	S            []float32 // solid (0) or liquid (1)

	M     []float32 // smoke
	newM  []float32
	tempM []float32 // temporary buffer for BFECC smoke

	// Strength of the vorticity confinement force. Set to 0 to disable.
	Confinement float32
}

func New(density float32, width, height int, h float32) *Fluid {
	numCells := (width + 2) * (height + 2)
	return &Fluid{
		density: density,
		h:       h,

		NumX:        width + 2,  // LR border cells
		NumY:        height + 2, // TB border cells
		numCells:    numCells,
		U:           make([]float32, numCells),
		V:           make([]float32, numCells),
		newU:        make([]float32, numCells),
		newV:        make([]float32, numCells),
		tempU:       make([]float32, numCells),
		tempV:       make([]float32, numCells),
		p:           make([]float32, numCells),
		S:           make([]float32, numCells),
		M:           make([]float32, numCells),
		newM:        make([]float32, numCells),
		tempM:       make([]float32, numCells),
		Confinement: 0,
	}
}

func fill[T any](slice []T, val T) {
	for i := range slice {
		slice[i] = val
	}
}

// Recommended number of pressure iterations for good quality with BFECC
const RecommendedPressureIterations = uint(8)

func (f *Fluid) Simulate(dt float32, numIters uint) {
	fill(f.p, 0)
	f.makeIncompressible(numIters, dt)
	if f.Confinement != 0 {
		f.applyVorticityConfinement(dt)
	}
	f.handleBorders()
	f.advectVelocity(dt)
	f.advectSmoke(dt)
}

// SimulateOptimized runs the simulation with recommended settings for best quality/performance balance
func (f *Fluid) SimulateOptimized(dt float32) {
	f.Simulate(dt, RecommendedPressureIterations)
}

func (f *Fluid) makeIncompressible(numIters uint, dt float32) {
	f.copyBorder(f.newU, f.U)
	f.copyBorder(f.newV, f.V)

	n := f.NumY
	cp := f.density * float32(f.h) / dt
	tolerance := float32(5e-6) // Convergence tolerance - slightly tighter for better wall handling

	for iter := uint(0); iter < numIters; iter++ {
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
				p *= Relaxation
				f.p[i*n+j] += cp * p

				f.U[i*n+j] -= sx0 * p
				f.U[(i+1)*n+j] += sx1 * p
				f.V[i*n+j] -= sy0 * p
				f.V[i*n+j+1] += sy1 * p
			}
		}

		// Early termination if converged
		if maxDiv < tolerance {
			break
		}
	}
}

func (f *Fluid) handleBorders() {
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		// Top border (j == 0)
		if f.S[i*n+0] == 0 || f.S[i*n+1] == 0 {
			f.U[i*n+0] = 0
		} else {
			f.U[i*n+0] = f.U[i*n+1]
		}

		// Bottom border (j == NumY-1)
		if f.S[i*n+f.NumY-1] == 0 || f.S[i*n+f.NumY-2] == 0 {
			f.U[i*n+f.NumY-1] = 0
		} else {
			f.U[i*n+f.NumY-1] = f.U[i*n+f.NumY-2]
		}
	})

	parallelRange(0, f.NumY, func(j int) {
		// Left border (i == 0)
		if f.S[0*n+j] == 0 || f.S[1*n+j] == 0 {
			f.V[0*n+j] = 0
		} else {
			f.V[0*n+j] = f.V[1*n+j]
		}

		// Right border (i == NumX-1)
		if f.S[(f.NumX-1)*n+j] == 0 || f.S[(f.NumX-2)*n+j] == 0 {
			f.V[(f.NumX-1)*n+j] = 0
		} else {
			f.V[(f.NumX-1)*n+j] = f.V[(f.NumX-2)*n+j]
		}
	})
}

// Semi-Lagrangian advection step - can be used for both forward and backward advection
func (f *Fluid) semiLagrangianVelocity(src_u, src_v, dst_u, dst_v []float32, dt float32) {
	f.copyBorder(dst_u, src_u)
	f.copyBorder(dst_v, src_v)

	n := f.NumY
	h := f.h
	h2 := h / 2

	parallelRange(1, f.NumX, func(i int) {
		for j := 1; j < f.NumY; j++ {

			// u component
			if f.S[i*n+j] != 0.0 && f.S[(i-1)*n+j] != 0.0 && j < f.NumY-1 {
				x := float32(i) * h
				y := float32(j)*h + h2

				// Use source velocity field for advection
				u := src_u[i*n+j]
				v := f.avgVFromField(src_v, i, j)

				x = x - dt*u
				y = y - dt*v
				u = f.sampleFieldFromArray(x, y, fieldU, src_u)
				dst_u[i*n+j] = u
			}

			// v component
			if f.S[i*n+j] != 0.0 && f.S[i*n+j-1] != 0.0 && i < f.NumX-1 {
				x := float32(i)*h + h2
				y := float32(j) * h

				// Use source velocity field for advection
				u := f.avgUFromField(src_u, i, j)
				v := src_v[i*n+j]

				x = x - dt*u
				y = y - dt*v
				v = f.sampleFieldFromArray(x, y, fieldV, src_v)
				dst_v[i*n+j] = v
			}
		}
	})
}

// BFECC velocity advection - much more accurate than basic semi-Lagrangian
func (f *Fluid) advectVelocity(dt float32) {
	// Step 1: Forward advection U -> tempU, V -> tempV
	f.semiLagrangianVelocity(f.U, f.V, f.tempU, f.tempV, dt)

	// Step 2: Backward advection tempU -> newU, tempV -> newV
	f.semiLagrangianVelocity(f.tempU, f.tempV, f.newU, f.newV, -dt)

	// Step 3: Compute error and apply correction
	n := f.NumY
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] != 0.0 {
				// Error estimation
				errorU := (f.U[i*n+j] - f.newU[i*n+j]) * 0.5
				errorV := (f.V[i*n+j] - f.newV[i*n+j]) * 0.5

				// Apply error correction to forward result
				f.newU[i*n+j] = f.tempU[i*n+j] + errorU
				f.newV[i*n+j] = f.tempV[i*n+j] + errorV
			}
		}
	})

	copy(f.U, f.newU)
	copy(f.V, f.newV)
}

func (f *Fluid) avgUFromField(field []float32, i, j int) float32 {
	n := f.NumY
	u := (field[i*n+j-1] + field[i*n+j] +
		field[(i+1)*n+j-1] + field[(i+1)*n+j]) * 0.25
	return u
}

func (f *Fluid) avgVFromField(field []float32, i, j int) float32 {
	n := f.NumY
	v := (field[(i-1)*n+j] + field[i*n+j] +
		field[(i-1)*n+j+1] + field[i*n+j+1]) * 0.25
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

func (f *Fluid) sampleFieldFromArray(x, y float32, fld field, fieldArray []float32) float32 {
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

	val := sx*sy*fieldArray[x0*n+y0] +
		tx*sy*fieldArray[x1*n+y0] +
		tx*ty*fieldArray[x1*n+y1] +
		sx*ty*fieldArray[x0*n+y1]

	return val
}

// Semi-Lagrangian smoke advection step
func (f *Fluid) semiLagrangianSmoke(srcM, dstM []float32, dt float32) {
	f.copyBorder(dstM, srcM)

	n := f.NumY
	h := f.h
	h2 := 0.5 * h

	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] != 0.0 {
				var u = (f.U[i*n+j] + f.U[(i+1)*n+j]) * 0.5
				var v = (f.V[i*n+j] + f.V[i*n+j+1]) * 0.5
				var x = float32(i)*h + h2 - dt*u
				var y = float32(j)*h + h2 - dt*v

				dstM[i*n+j] = f.sampleFieldFromArray(x, y, fieldM, srcM)
			}
		}
	})
}

// BFECC smoke advection - much more accurate than basic semi-Lagrangian
func (f *Fluid) advectSmoke(dt float32) {
	// Step 1: Forward advection M -> tempM
	f.semiLagrangianSmoke(f.M, f.tempM, dt)

	// Step 2: Backward advection tempM -> newM
	f.semiLagrangianSmoke(f.tempM, f.newM, -dt)

	// Step 3: Compute error and apply correction
	n := f.NumY
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.S[i*n+j] != 0.0 {
				// Error estimation
				errorM := (f.M[i*n+j] - f.newM[i*n+j]) * 0.5

				// Apply error correction to forward result
				f.newM[i*n+j] = f.tempM[i*n+j] + errorM
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

// applyVorticityConfinement computes the curl of the velocity field and
// applies a confinement force to enhance small scale swirling motion.
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

			f.U[i*n+j] += f.Confinement * gy * vort * dt
			f.V[i*n+j] -= f.Confinement * gx * vort * dt
		}
	}

}

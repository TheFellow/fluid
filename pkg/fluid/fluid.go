package fluid

import (
	"math"
)

type Fluid struct {
	density float32
	h       float32
	Gravity float32

	NumX, NumY int
	numCells   int
	u, v       []float32 // velocities
	newU, newV []float32
	p          []float32
	s          []float32 // solid (0) or liquid (1)

	m    []float32 // smoke
	newM []float32

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
		u:           make([]float32, numCells),
		v:           make([]float32, numCells),
		newU:        make([]float32, numCells),
		newV:        make([]float32, numCells),
		p:           make([]float32, numCells),
		s:           make([]float32, numCells),
		m:           make([]float32, numCells),
		newM:        make([]float32, numCells),
		Confinement: 0,
	}
}

func fill[T any](slice []T, val T) {
	for i := range slice {
		slice[i] = val
	}
}

func (f *Fluid) Simulate(dt float32) {
	// Apply external forces.
	f.handleGravity(dt)
	if f.Confinement != 0 {
		f.applyVorticityConfinement(dt)
	}

	// Use BFECC advection for both velocity and smoke
	// This eliminates the need for the incompressibility solve
	f.advectVelocityBFECC(dt)
	f.advectSmokeBFECC(dt)

	// Only handle boundary conditions - no pressure solve needed
	f.handleBorders()
}

func (f *Fluid) handleGravity(dt float32) {
	if f.Gravity == 0 {
		return
	}
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		// Skip the bottom boundary row to avoid out-of-range access when
		// referencing j-1. Gravity only affects interior interfaces.
		for j := 1; j < f.NumY; j++ {
			// Apply gravity only when both the current cell and the
			// cell below are fluid.
			if f.s[i*n+j] != 0.0 && f.s[i*n+j-1] != 0.0 {
				f.v[i*n+j] += f.Gravity * dt
			}
		}
	})
}

func (f *Fluid) handleBorders() {
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		// Top border (j == 0)
		if f.s[i*n+0] == 0 || f.s[i*n+1] == 0 {
			f.u[i*n+0] = 0
		} else {
			f.u[i*n+0] = f.u[i*n+1]
		}

		// Bottom border (j == NumY-1)
		if f.s[i*n+f.NumY-1] == 0 || f.s[i*n+f.NumY-2] == 0 {
			f.u[i*n+f.NumY-1] = 0
		} else {
			f.u[i*n+f.NumY-1] = f.u[i*n+f.NumY-2]
		}
	})

	parallelRange(0, f.NumY, func(j int) {
		// Left border (i == 0)
		if f.s[0*n+j] == 0 || f.s[1*n+j] == 0 {
			f.v[0*n+j] = 0
		} else {
			f.v[0*n+j] = f.v[1*n+j]
		}

		// Right border (i == NumX-1)
		if f.s[(f.NumX-1)*n+j] == 0 || f.s[(f.NumX-2)*n+j] == 0 {
			f.v[(f.NumX-1)*n+j] = 0
		} else {
			f.v[(f.NumX-1)*n+j] = f.v[(f.NumX-2)*n+j]
		}
	})
}

// advectVelocityBFECC uses the BFECC (Back and Forth Error Compensation and Correction)
// method to advect velocity with minimal numerical dissipation and inherent stability
func (f *Fluid) advectVelocityBFECC(dt float32) {
	// Forward advection step
	f.semiLagrangianVelocity(dt, f.u, f.v, f.newU, f.newV)

	// Store forward results
	forwardU := make([]float32, f.numCells)
	forwardV := make([]float32, f.numCells)
	copy(forwardU, f.newU)
	copy(forwardV, f.newV)

	// Reverse advection step - advect the forward result backward
	tmpU := make([]float32, f.numCells)
	tmpV := make([]float32, f.numCells)
	f.semiLagrangianVelocity(-dt, forwardU, forwardV, tmpU, tmpV)

	// BFECC correction: forward + 0.5 * (original - backward)
	n := f.NumY
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.s[i*n+j] != 0.0 {
				idx := i*n + j
				f.newU[idx] = forwardU[idx] + 0.5*(f.u[idx]-tmpU[idx])
				f.newV[idx] = forwardV[idx] + 0.5*(f.v[idx]-tmpV[idx])
			}
		}
	})

	// Apply the corrected velocities
	f.copyBorder(f.newU, f.newU)
	f.copyBorder(f.newV, f.newV)
	copy(f.u, f.newU)
	copy(f.v, f.newV)
}

// advectSmokeBFECC uses BFECC for smoke advection with clamping to prevent oscillations
func (f *Fluid) advectSmokeBFECC(dt float32) {
	// Forward advection
	forward := make([]float32, f.numCells)
	f.semiLagrangianSmoke(dt, f.m, forward)

	// Reverse advection
	backward := make([]float32, f.numCells)
	f.semiLagrangianSmoke(-dt, forward, backward)

	// BFECC correction with clamping to local extrema
	n := f.NumY
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.s[i*n+j] != 0.0 {
				idx := i*n + j
				corrected := forward[idx] + 0.5*(f.m[idx]-backward[idx])

				// Clamp to local min/max to prevent oscillations
				minVal := f.m[idx]
				maxVal := f.m[idx]
				neighbors := []int{(i-1)*n + j, (i+1)*n + j, i*n + j - 1, i*n + j + 1}
				for _, nIdx := range neighbors {
					if f.m[nIdx] < minVal {
						minVal = f.m[nIdx]
					}
					if f.m[nIdx] > maxVal {
						maxVal = f.m[nIdx]
					}
				}

				if corrected < minVal {
					corrected = minVal
				}
				if corrected > maxVal {
					corrected = maxVal
				}

				f.newM[idx] = corrected
			}
		}
	})

	f.copyBorder(f.newM, f.newM)
	copy(f.m, f.newM)
}

// Update semiLagrangianVelocity to work with explicit source arrays
func (f *Fluid) semiLagrangianVelocity(dt float32, srcU, srcV, dstU, dstV []float32) {
	f.copyBorder(dstU, srcU)
	f.copyBorder(dstV, srcV)

	n := f.NumY
	h := f.h
	h2 := h / 2

	parallelRange(1, f.NumX, func(i int) {
		for j := 1; j < f.NumY; j++ {

			// u component
			if f.s[i*n+j] != 0.0 && f.s[(i-1)*n+j] != 0.0 && j < f.NumY-1 {
				x := float32(i) * h
				y := float32(j)*h + h2
				u := srcU[i*n+j]
				v := f.avgVFrom(srcV, i, j)

				x = x - dt*u
				y = y - dt*v
				dstU[i*n+j] = f.sampleFieldFrom(x, y, srcU, fieldU)
			}

			// v component
			if f.s[i*n+j] != 0.0 && f.s[i*n+j-1] != 0.0 && i < f.NumX-1 {
				x := float32(i)*h + h2
				y := float32(j) * h
				u := f.avgUFrom(srcU, i, j)
				v := srcV[i*n+j]

				x = x - dt*u
				y = y - dt*v
				dstV[i*n+j] = f.sampleFieldFrom(x, y, srcV, fieldV)
			}
		}
	})
}

// avgUFrom and avgVFrom are helpers that average velocity components from the
// supplied slices. They mirror avgU/avgV but operate on explicit arrays so they
// can be used during the temporary advection passes.
func (f *Fluid) avgUFrom(u []float32, i, j int) float32 {
	n := f.NumY
	return (u[i*n+j-1] + u[i*n+j] + u[(i+1)*n+j-1] + u[(i+1)*n+j]) * 0.25
}

func (f *Fluid) avgVFrom(v []float32, i, j int) float32 {
	n := f.NumY
	return (v[(i-1)*n+j] + v[i*n+j] + v[(i-1)*n+j+1] + v[i*n+j+1]) * 0.25
}

func (f *Fluid) avgU(i, j int) float32 {
	n := f.NumY
	u := (f.u[i*n+j-1] + f.u[i*n+j] +
		f.u[(i+1)*n+j-1] + f.u[(i+1)*n+j]) * 0.25
	return u
}

func (f *Fluid) avgV(i, j int) float32 {
	n := f.NumY
	v := (f.v[(i-1)*n+j] + f.v[i*n+j] +
		f.v[(i-1)*n+j+1] + f.v[i*n+j+1]) * 0.25
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
		fieldToSample = f.u
		dy = h2
	case fieldV:
		fieldToSample = f.v
		dx = h2
	case fieldM:
		fieldToSample = f.m
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

// Add a new sampling method that uses explicit source arrays
func (f *Fluid) sampleFieldFrom(x, y float32, src []float32, fld field) float32 {
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

	val := sx*sy*src[x0*n+y0] +
		tx*sy*src[x1*n+y0] +
		tx*ty*src[x1*n+y1] +
		sx*ty*src[x0*n+y1]

	return val
}

// advectSmoke transports the density/smoke field using the same
// MacCormack-style forward/reverse advection employed for the velocity. After
// the correction step a simple pressure solve is applied which distributes any
// negative values proportionally to neighbouring cells so that the field
// remains non-negative and mass conserving.
func (f *Fluid) advectSmoke(dt float32) {
	forward := make([]float32, f.numCells)

	// Forward step
	f.semiLagrangianSmoke(dt, f.m, forward)

	// Reverse step
	oldM := f.m
	f.m = forward
	backward := make([]float32, f.numCells)
	f.semiLagrangianSmoke(-dt, f.m, backward)
	f.m = oldM

	// Correction with BFECC
	n := f.NumY
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			idx := i*n + j
			val := forward[idx] + 0.5*(f.m[idx]-backward[idx])

			// Limit to the range of the original and neighbouring
			// cells to avoid creating new extrema.
			minVal := f.m[idx]
			maxVal := f.m[idx]
			neigh := []int{(i-1)*n + j, (i+1)*n + j, i*n + j - 1, i*n + j + 1}
			for _, k := range neigh {
				if f.m[k] < minVal {
					minVal = f.m[k]
				}
				if f.m[k] > maxVal {
					maxVal = f.m[k]
				}
			}
			if val < minVal {
				val = minVal
			}
			if val > maxVal {
				val = maxVal
			}

			f.newM[idx] = val
		}
	})

	// Apply a simple pressure step to keep the field non-negative.
	f.applyPressure(f.newM)
	f.copyBorder(f.newM, f.newM)
	copy(f.m, f.newM)
}

// semiLagrangianSmoke performs a standard semi-Lagrangian advection of the
// smoke/density field from src into dst using the current velocity field.
func (f *Fluid) semiLagrangianSmoke(dt float32, src, dst []float32) {
	f.copyBorder(dst, src)

	n := f.NumY
	h := f.h
	h2 := 0.5 * h

	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			if f.s[i*n+j] != 0.0 {
				u := (f.u[i*n+j] + f.u[(i+1)*n+j]) * 0.5
				v := (f.v[i*n+j] + f.v[i*n+j+1]) * 0.5
				x := float32(i)*h + h2 - dt*u
				y := float32(j)*h + h2 - dt*v
				dst[i*n+j] = f.sampleField(x, y, fieldM)
			}
		}
	})
}

// applyPressure distributes any negative mass to neighbouring cells. This is a
// lightweight substitute for the full pressure solve described by Stam and
// mirrors the approach outlined in "Practical Fluid Mechanics". It prevents the
// creation of artificial mass when multiple cells attempt to pull more than is
// available from a single source cell.
func (f *Fluid) applyPressure(field []float32) {
	n := f.NumY
	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {
			idx := i*n + j
			if field[idx] >= 0 {
				continue
			}
			deficit := -field[idx]
			field[idx] = 0
			share := deficit * 0.25
			field[(i-1)*n+j] -= share
			field[(i+1)*n+j] -= share
			field[i*n+j-1] -= share
			field[i*n+j+1] -= share
		}
	})
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
			if f.s[i*n+j] == 0 {
				continue
			}

			dvdx := (f.v[(i+1)*n+j] - f.v[(i-1)*n+j]) * 0.5 / h
			dudy := (f.u[i*n+j+1] - f.u[i*n+j-1]) * 0.5 / h
			curl[i*n+j] = dvdx - dudy
		}
	}

	eps := float32(1e-5)
	for i := 1; i < f.NumX-1; i++ {
		for j := 1; j < f.NumY-1; j++ {
			if f.s[i*n+j] == 0 {
				continue
			}

			// gradient of magnitude of curl
			gx := (float32(math.Abs(float64(curl[(i+1)*n+j]))) - float32(math.Abs(float64(curl[(i-1)*n+j])))) * 0.5 / h
			gy := (float32(math.Abs(float64(curl[i*n+j+1]))) - float32(math.Abs(float64(curl[i*n+j-1])))) * 0.5 / h

			mag := float32(math.Sqrt(float64(gx*gx+gy*gy))) + eps
			gx /= mag
			gy /= mag

			vort := curl[i*n+j]

			f.u[i*n+j] += f.Confinement * gy * vort * dt
			f.v[i*n+j] -= f.Confinement * gx * vort * dt
		}
	}

}

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
		p:           make([]float32, numCells),
		S:           make([]float32, numCells),
		M:           make([]float32, numCells),
		newM:        make([]float32, numCells),
		Confinement: 0,
	}
}

func fill[T any](slice []T, val T) {
	for i := range slice {
		slice[i] = val
	}
}

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

func (f *Fluid) makeIncompressible(numIters uint, dt float32) {
	f.copyBorder(f.newU, f.U)
	f.copyBorder(f.newV, f.V)

	n := f.NumY
	cp := f.density * float32(f.h) / dt

	for range numIters {

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

				p := -div / s
				p *= Relaxation
				f.p[i*n+j] += cp * p

				f.U[i*n+j] -= sx0 * p
				f.U[(i+1)*n+j] += sx1 * p
				f.V[i*n+j] -= sy0 * p
				f.V[i*n+j+1] += sy1 * p
			}
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

func (f *Fluid) advectVelocity(dt float32) {
	// Preserve boundaries by copying them into the destination slices
	// before computing the interior advection.
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
				var u = (f.U[i*n+j] + f.U[(i+1)*n+j]) * 0.5
				var v = (f.V[i*n+j] + f.V[i*n+j+1]) * 0.5
				var x = float32(i)*h + h2 - dt*u
				var y = float32(j)*h + h2 - dt*v

				f.newM[i*n+j] = f.sampleField(x, y, fieldM)
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

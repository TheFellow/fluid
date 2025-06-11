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

func (f *Fluid) Simulate(dt, gravity float32, numIters uint) {
	f.handleGravity(dt)
	fill(f.p, 0)
	f.makeIncompressible(numIters, dt)
	if f.Confinement != 0 {
		f.applyVorticityConfinement(dt)
	}
	f.handleBorders()
	f.advectVelocity(dt)
	f.advectSmoke(dt)
}

func (f *Fluid) handleGravity(dt float32) {
	if f.Gravity == 0 {
		return
	}
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		for j := 0; j < f.NumY; j++ {
			if i*n+j-1 < 0 {
				continue
			}
			if f.s[i*n+j] != 0.0 && f.s[i*n+j-1] != 0.0 {
				f.v[i*n+j] += f.Gravity * dt
			}
		}
	})
}

func (f *Fluid) makeIncompressible(numIters uint, dt float32) {
	f.copyBorder(f.newU, f.u)
	f.copyBorder(f.newV, f.v)

	n := f.NumY
	cp := f.density * float32(f.h) / dt

	for range numIters {

		for i := 1; i < f.NumX-1; i++ {
			for j := 1; j < f.NumY-1; j++ {

				// If the cell is solid, nothing to do...
				if f.s[i*n+j] == 0 {
					continue
				}

				sx0 := f.s[(i-1)*n+j]
				sx1 := f.s[(i+1)*n+j]
				sy0 := f.s[i*n+j-1]
				sy1 := f.s[i*n+j+1]
				s := sx0 + sx1 + sy0 + sy1
				if s == 0 { // All adjacent cells are solid, nothing to do...
					continue
				}

				div := f.u[(i+1)*n+j] - f.u[i*n+j] +
					f.v[i*n+j+1] - f.v[i*n+j]

				p := -div / s
				p *= Relaxation
				f.p[i*n+j] += cp * p

				f.u[i*n+j] -= sx0 * p
				f.u[(i+1)*n+j] += sx1 * p
				f.v[i*n+j] -= sy0 * p
				f.v[i*n+j+1] += sy1 * p
			}
		}
	}
}

func (f *Fluid) handleBorders() {
	n := f.NumY
	parallelRange(0, f.NumX, func(i int) {
		f.u[i*n+0] = f.u[i*n+1]               // top border
		f.u[i*n+f.NumY-1] = f.u[i*n+f.NumY-2] // bottom border
	})

	parallelRange(0, f.NumY, func(j int) {
		f.v[0*n+j] = f.v[1*n+j]                   // left border
		f.v[(f.NumX-1)*n+j] = f.v[(f.NumX-2)*n+j] // right border
	})
}

func (f *Fluid) advectVelocity(dt float32) {
	// Preserve boundaries by copying them into the destination slices
	// before computing the interior advection.
	f.copyBorder(f.newU, f.u)
	f.copyBorder(f.newV, f.v)

	n := f.NumY
	h := f.h
	h2 := h / 2

	parallelRange(1, f.NumX, func(i int) {
		for j := 1; j < f.NumY; j++ {

			// u component
			if f.s[i*n+j] != 0.0 && f.s[(i-1)*n+j] != 0.0 && j < f.NumY-1 {
				x := float32(i) * h
				y := float32(j)*h + h2
				u := f.u[i*n+j]
				v := f.avgV(i, j)

				x = x - dt*u
				y = y - dt*v
				u = f.sampleField(x, y, fieldU)
				f.newU[i*n+j] = u
			}

			// v component
			if f.s[i*n+j] != 0.0 && f.s[i*n+j-1] != 0.0 && i < f.NumX-1 {
				x := float32(i)*h + h2
				y := float32(j) * h
				u := f.avgU(i, j)
				v := f.v[i*n+j]

				x = x - dt*u
				y = y - dt*v
				v = f.sampleField(x, y, fieldV)
				f.newV[i*n+j] = v
			}
		}
	})

	copy(f.u, f.newU)
	copy(f.v, f.newV)
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

func (f *Fluid) advectSmoke(dt float32) {
	// Copy border cells first so advection doesn't alter boundary values.

	n := f.NumY
	h := f.h
	h2 := 0.5 * h

	parallelRange(1, f.NumX-1, func(i int) {
		for j := 1; j < f.NumY-1; j++ {

			if f.s[i*n+j] != 0.0 {
				var u = (f.u[i*n+j] + f.u[(i+1)*n+j]) * 0.5
				var v = (f.v[i*n+j] + f.v[i*n+j+1]) * 0.5
				var x = float32(i)*h + h2 - dt*u
				var y = float32(j)*h + h2 - dt*v

				f.newM[i*n+j] = f.sampleField(x, y, fieldM)
			}
		}
	})

	copy(f.m, f.newM)
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

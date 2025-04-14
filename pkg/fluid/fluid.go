package fluid

import (
	"fmt"
	"math"
)

var (
	Relaxation float32 = 1.9
)

type Fluid struct {
	density float32
	h       float32

	numX, numY int
	numCells   int
	u, v       []float32 // velocities
	newU, newV []float32
	p          []float32
	s          []float32 // solid (0) or liquid (1)

	m    []float32 // smoke
	newM []float32
}

func New(density float32, width, height int, h float32) *Fluid {
	numCells := (width + 2) * (height + 2)
	m := make([]float32, numCells)
	fill(m, 1.0)
	return &Fluid{
		density: density,
		h:       h,

		numX:     width + 2,  // LR border cells
		numY:     height + 2, // TB border cells
		numCells: numCells,
		u:        make([]float32, numCells),
		v:        make([]float32, numCells),
		newU:     make([]float32, numCells),
		newV:     make([]float32, numCells),
		p:        make([]float32, numCells),
		s:        make([]float32, numCells),
		m:        m,
		newM:     make([]float32, numCells),
	}
}

func fill[T any](slice []T, val T) {
	for i := range slice {
		slice[i] = val
	}
}

func (f *Fluid) Simulate(dt, gravity float32, numIters uint) {
	f.handleGravity(dt, gravity)
	fill(f.p, 0)
	f.makeIncompressible(numIters, dt)
	f.handleBorders()
	f.advectVelocity(dt)
	f.advectSmoke(dt)
}

func (f *Fluid) handleGravity(dt, gravity float32) {
	n := f.numY
	for i := 1; i < f.numX-1; i++ {
		for j := 1; j < f.numY-1; j++ {
			if f.s[i*n+j] != 0.0 && f.s[i*n+j-1] != 0.0 {
				f.v[i*n+j] += gravity * dt
			}
		}
	}
}

func (f *Fluid) makeIncompressible(numIters uint, dt float32) {
	n := f.numY
	cp := f.density * float32(f.h) / dt

	for range numIters {

		for i := 1; i < f.numX-1; i++ {
			for j := 1; j < f.numY-1; j++ {
				if i == 1 && j == 1 || i == 1 && j == f.numY-2 {
					xxx := fmt.Sprintf("b")
					_ = xxx
				}

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
	n := f.numY
	for i := range f.numX {
		f.u[i*n+0] = f.u[i*n+1]               // top border
		f.u[i*n+f.numY-1] = f.u[i*n+f.numY-2] // bottom border
	}

	for j := range f.numY {
		f.v[0*n+j] = f.v[1*n+j]                   // left border
		f.v[(f.numX-1)*n+j] = f.v[(f.numX-2)*n+j] // right border
	}
}

func (f *Fluid) advectVelocity(dt float32) {
	// TODO: Remove f clearing of the array?
	// for i := range f.numCells {
	// 	f.newU[i] = f.u[i]
	// }
	// for j := range f.numCells {
	// 	f.newV[j] = f.v[j]
	// }

	n := f.numY
	h := f.h
	h2 := h / 2

	for i := 0; i < f.numX; i++ {
		for j := 0; j < f.numY; j++ {

			// u component
			if f.s[i*n+j] != 0.0 && f.s[(i-1)*n+j] != 0.0 {
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
			if f.s[i*n+j] != 0.0 && f.s[i*n+j-1] != 0.0 {
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
	}

	for i := range f.numCells {
		f.u[i] = f.newU[i]
	}
	for j := range f.numCells {
		f.v[j] = f.newV[j]
	}
}

func (f *Fluid) avgU(i, j int) float32 {
	n := f.numY
	u := (f.u[i*n+j-1] + f.u[i*n+j] +
		f.u[(i+1)*n+j-1] + f.u[(i+1)*n+j]) * 0.25
	return u
}

func (f *Fluid) avgV(i, j int) float32 {
	n := f.numY
	v := (f.v[(i-1)*n+j] + f.v[i*n+j] +
		f.v[(i-1)*n+j+1] + f.v[i*n+j+1]) * 0.25
	return v
}

type field int

const (
	fieldU field = iota
	fieldV
	fieldS
)

func (f *Fluid) sampleField(x, y float32, fld field) float32 {
	n := f.numY
	h := f.h
	h1 := float32(1.0 / h)
	h2 := float32(0.5 * h)

	x = max(min(x, float32(f.numX)*h), h)
	y = max(min(y, float32(f.numY)*h), h)

	dx, dy := float32(0.0), float32(0.0)

	var fieldToSample []float32
	switch fld {
	case fieldU:
		fieldToSample = f.u
		dy = h2
	case fieldV:
		fieldToSample = f.v
		dx = h2
	case fieldS:
		fieldToSample = f.m
		dx, dy = h2, h2
	}

	x0 := min(int(math.Floor(float64((x-dx)*h1))), f.numX-1)
	tx := ((x - dx) - float32(x0)*h) * h1
	x1 := min(x0+1, f.numX-1)

	y0 := min(int(math.Floor(float64((y-dy)*h1))), f.numY-1)
	ty := ((y - dy) - float32(y0)*h) * h1
	y1 := min(y0+1, f.numY-1)

	sx := 1.0 - tx
	sy := 1.0 - ty

	val := sx*sy*fieldToSample[x0*n+y0] +
		tx*sy*fieldToSample[x1*n+y0] +
		tx*ty*fieldToSample[x1*n+y1] +
		sx*ty*fieldToSample[x0*n+y1]

	return val
}

func (f *Fluid) advectSmoke(dt float32) {
	for i := range f.numCells {
		f.newM[i] = f.m[i]
	}

	n := f.numY
	h := float32(f.h)
	h2 := 0.5 * float32(h)

	for i := 1; i < f.numX-1; i++ {
		for j := 1; j < f.numY-1; j++ {

			if f.s[i*n+j] != 0.0 {
				var u = (f.u[i*n+j] + f.u[(i+1)*n+j]) * 0.5
				var v = (f.v[i*n+j] + f.v[i*n+j+1]) * 0.5
				var x = float32(i)*h + h2 - dt*u
				var y = float32(j)*h + h2 - dt*v

				f.newM[i*n+j] = f.sampleField(x, y, fieldS)
			}
		}
	}

	for i := range f.numCells {
		f.m[i] = f.newM[i]
	}
}

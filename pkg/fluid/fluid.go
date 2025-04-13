package fluid

type Fluid struct {
	density float32
	h       float32

	numX, numY int
	u, v       []float32 // velocities
	p          []float32
	s          []int // solid (0) or liquid (1)
}

func New(density float32, width, height int, h float32) *Fluid {
	numCells := (width + 2) * (height + 2)
	return &Fluid{
		density: density,
		h:       h,

		numX: width + 2,
		numY: height + 2,
		u:    make([]float32, 0, numCells),
		v:    make([]float32, 0, numCells),
		p:    make([]float32, 0, numCells),
		s:    make([]int, 0, numCells),
	}
}

func (f *Fluid) Simulate(dt, gravity float32, numIters uint) {
	f.handleGravity(dt, gravity)
	for i := range f.p {
		f.p[i] = 0
	}
	f.makeIncompressible(numIters, dt)
}

func (f *Fluid) handleGravity(dt, gravity float32) {
	n := f.numY
	for i := 1; i < f.numX; i++ {
		for j := 1; j < f.numY-1; j++ {
			cell := i*n + j
			if f.s[cell] != 0.0 && f.s[cell-1] != 0.0 {
				f.v[cell] += gravity * dt
			}
		}
	}
}

func (f *Fluid) makeIncompressible(numIters uint, dt float32) {
	n := f.numY
	cp := f.density * f.h / dt
	_ = cp
	for range numIters {

		for i := 1; i < f.numX-1; i++ {
			for j := 1; j < f.numY-1; j++ {
				// If the cell is solid, nothing to do...
				if f.s[i*n+j] == 0 {
					continue
				}

				sx0 := f.s[(i-1)*n+j]
				sx1 := f.s[(i+1)*n+j]
				sy0 := f.s[i*n+j-1]
				sy1 := f.s[i*n+j+1]
				s := sx0 + sx1 + sy0 + sy1
				if s == 0 {
					continue
				}
				// TODO: Pickup here
			}
		}
	}
}

package fluid

import "fmt"

func (f *Fluid) SetSolid(i, j int, value bool) {
	if i < 0 || i >= f.NumX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.NumY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.NumY + j
	if value {
		f.S[cell] = 0.0
	} else {
		f.S[cell] = 1.0
	}

	// When a cell becomes solid, clear the associated velocities to
	// immediately enforce the wall boundary. This prevents fluid from
	// continuing to flow through or along newly created walls.
	if value {
		n := f.NumY
		// Velocities are stored on a staggered grid. Zero out the four
		// edges surrounding this cell, checking bounds to avoid
		// out-of-range accesses.
		f.U[i*n+j] = 0
		if i+1 < f.NumX {
			f.U[(i+1)*n+j] = 0
		}
		f.V[i*n+j] = 0
		if j+1 < f.NumY {
			f.V[i*n+j+1] = 0
		}

		// Clear any in-flight velocities that may still reside in the
		// temporary buffers used during simulation steps. Without this,
		// advecting the velocity field would reintroduce non-zero values
		// after the wall is drawn.
		f.newU[i*n+j] = 0
		if i+1 < f.NumX {
			f.newU[(i+1)*n+j] = 0
		}
		f.newV[i*n+j] = 0
		if j+1 < f.NumY {
			f.newV[i*n+j+1] = 0
		}
	}
}

func (f *Fluid) IsSolid(i, j int) bool {
	if i < 0 || i >= f.NumX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.NumY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.NumY + j
	return f.S[cell] == 0.0
}

func (f *Fluid) SetVelocity(i, j int, u, v float32) {
	if i < 0 || i >= f.NumX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.NumY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.NumY + j
	f.U[cell] = u
	f.V[cell] = v
}

func (f *Fluid) AddSmoke(i, j int, smoke float32) {
	if i < 0 || i >= f.NumX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.NumY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.NumY + j
	f.M[cell] += smoke
}

func (f *Fluid) SetGravity(on bool) {
	if on {
		f.Gravity = -9.81
	} else {
		f.Gravity = 0.0
	}
}

func (f *Fluid) Reset() {
	fill(f.U, 0.0)
	fill(f.V, 0.0)
	fill(f.newU, 0.0)
	fill(f.newV, 0.0)
	fill(f.p, 0.0)
	fill(f.M, 0.0)
	fill(f.newM, 0.0)
}

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
		f.s[cell] = 0.0
	} else {
		f.s[cell] = 1.0
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
	return f.s[cell] == 0.0
}

func (f *Fluid) SetVelocity(i, j int, u, v float32) {
	if i < 0 || i >= f.NumX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.NumY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.NumY + j
	f.u[cell] = u
	f.v[cell] = v
}

func (f *Fluid) AddSmoke(i, j int, smoke float32) {
	if i < 0 || i >= f.NumX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.NumY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.NumY + j
	f.m[cell] += smoke
}

func (f *Fluid) SetGravity(on bool) {
	if on {
		f.Gravity = -9.81
	} else {
		f.Gravity = 0.0
	}
}

func (f *Fluid) Reset() {
	fill(f.u, 0.0)
	fill(f.v, 0.0)
	fill(f.m, 0.0)
}

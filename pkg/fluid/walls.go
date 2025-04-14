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

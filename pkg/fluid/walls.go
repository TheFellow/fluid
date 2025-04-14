package fluid

import "fmt"

func (f *Fluid) SetSolid(i, j int, value bool) {
	if i < 0 || i >= f.numX {
		panic(fmt.Sprintf("invalid x-index: %d", i))
	}
	if j < 0 || j >= f.numY {
		panic(fmt.Sprintf("invalid y-index: %d", j))
	}
	cell := i*f.numY + j
	if value {
		f.s[cell] = 0.0
	} else {
		f.s[cell] = 1.0
	}
}

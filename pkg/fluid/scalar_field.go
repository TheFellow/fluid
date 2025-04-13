package fluid

import (
	"fmt"
)

type ScalarField struct {
	numX, numY int
	values     []float32
}

func (s ScalarField) Value(i, j int) (float32, error) {
	if i < 0 || i >= s.numX {
		return 0.0, fmt.Errorf("x index out of range, must be between 0 and %d", s.numX-1)
	}
	if j < 0 || j >= s.numY {
		return 0.0, fmt.Errorf("y index out of range, must be between 0 and %d", s.numY-1)
	}

	return s.values[i*s.numY+j], nil
}

package fluid

import (
	"fmt"
)

type ScalarField struct {
	NumX, NumY int
	values     []float32
}

func (s ScalarField) Value(i, j int) (float32, error) {
	if i < 0 || i >= s.NumX {
		return 0.0, fmt.Errorf("x index (%d) out of range, must be between 0 and %d", i, s.NumX-1)
	}
	if j < 0 || j >= s.NumY {
		return 0.0, fmt.Errorf("y index (%d) out of range, must be between 0 and %d", j, s.NumY-1)
	}

	return s.values[i*s.NumY+j], nil
}

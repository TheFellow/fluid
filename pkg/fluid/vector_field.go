package fluid

import "fmt"

type VectorField struct {
	NumX, NumY       int
	valuesU, valuesV []float32
}

func (v VectorField) Value(i, j int) (float32, float32, error) {
	if i < 0 || i >= v.NumX {
		return 0.0, 0.0, fmt.Errorf("x index out of range, must be between 0 and %d", v.NumX-1)
	}
	if j < 0 || j >= v.NumY {
		return 0.0, 0.0, fmt.Errorf("y index out of range, must be between 0 and %d", v.NumY-1)
	}

	return v.valuesU[i*v.NumY+j], v.valuesV[i*v.NumY+j], nil
}

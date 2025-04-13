package fluid

import "fmt"

type VectorField struct {
	numX, numY       int
	valuesU, valuesV []float32
}

func (v VectorField) Value(i, j int) (float32, float32, error) {
	if i < 0 || i >= v.numX {
		return 0.0, 0.0, fmt.Errorf("x index out of range, must be between 0 and %d", v.numX-1)
	}
	if j < 0 || j >= v.numY {
		return 0.0, 0.0, fmt.Errorf("y index out of range, must be between 0 and %d", v.numY-1)
	}

	return v.valuesU[i*v.numY+j], v.valuesV[i*v.numY+j], nil
}

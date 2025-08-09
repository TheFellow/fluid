package fluid

import "math"

func (f *Fluid) Smoke() ScalarField {
	minValue := float32(math.MaxFloat32)
	maxValue := float32(-(math.MaxFloat32 - 1))
	for i := range f.numCells {
		if f.M[i] < minValue {
			minValue = f.M[i]
		}
		if f.M[i] > maxValue {
			maxValue = f.M[i]
		}
	}
	s := ScalarField{
		NumX:     f.NumX,
		NumY:     f.NumY,
		MaxValue: maxValue,
		MinValue: minValue,
		values:   f.M,
	}
	return s
}

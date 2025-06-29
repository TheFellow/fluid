package fluid

import "math"

func (f *Fluid) Smoke() ScalarField {
	minValue := float32(math.MaxFloat32)
	maxValue := float32(-(math.MaxFloat32 - 1))
	for i := range f.numCells {
		if f.m[i] < minValue {
			minValue = f.m[i]
		}
		if f.m[i] > maxValue {
			maxValue = f.m[i]
		}
	}
	s := ScalarField{
		NumX:     f.NumX,
		NumY:     f.NumY,
		MaxValue: maxValue,
		MinValue: minValue,
		values:   f.m,
	}
	return s
}

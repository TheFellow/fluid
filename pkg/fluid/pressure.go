package fluid

import "math"

func (f *Fluid) Pressure() ScalarField {
	minValue := float32(math.MaxFloat32)
	maxValue := float32(-(math.MaxFloat32 - 1))
	pCopy := make([]float32, f.numCells)
	for i := range f.numCells {
		pCopy[i] = f.p[i]
		if f.p[i] < minValue {
			minValue = f.p[i]
		}
		if f.p[i] > maxValue {
			maxValue = f.p[i]
		}
	}
	s := ScalarField{
		NumX:     f.NumX,
		NumY:     f.NumY,
		MaxValue: maxValue,
		MinValue: minValue,
		values:   f.p,
	}
	return s
}

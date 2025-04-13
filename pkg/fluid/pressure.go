package fluid

func (f *Fluid) Pressure() ScalarField {
	pCopy := make([]float32, f.numCells)
	for i := range f.numCells {
		pCopy[i] = f.p[i]
	}
	s := ScalarField{
		numX:   f.numX,
		numY:   f.numY,
		values: pCopy,
	}
	return s
}

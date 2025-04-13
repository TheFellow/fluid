package fluid

func (f *Fluid) Smoke() ScalarField {
	mCopy := make([]float32, f.numCells)
	for i := range f.numCells {
		mCopy[i] = f.m[i]
	}
	s := ScalarField{
		numX:   f.numX,
		numY:   f.numY,
		values: mCopy,
	}
	return s
}

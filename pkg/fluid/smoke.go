package fluid

func (f *Fluid) Smoke() ScalarField {
	// mCopy := make([]float32, f.numCells)
	// for i := range f.numCells {
	// 	mCopy[i] = f.m[i]
	// }
	s := ScalarField{
		NumX:   f.NumX,
		NumY:   f.NumY,
		values: f.m,
	}
	return s
}

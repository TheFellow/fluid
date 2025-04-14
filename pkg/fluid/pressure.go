package fluid

func (f *Fluid) Pressure() ScalarField {
	// pCopy := make([]float32, f.numCells)
	// for i := range f.numCells {
	// 	pCopy[i] = f.p[i]
	// }
	s := ScalarField{
		NumX:   f.NumX,
		NumY:   f.NumY,
		values: f.p,
	}
	return s
}

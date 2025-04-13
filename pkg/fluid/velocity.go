package fluid

func (f *Fluid) Velocity() VectorField {
	uCopy := make([]float32, f.numCells)
	for i := range f.numCells {
		uCopy[i] = f.u[i]
	}
	vCopy := make([]float32, f.numCells)
	for i := range f.numCells {
		vCopy[i] = f.v[i]
	}
	return VectorField{
		numX:    f.numX,
		numY:    f.numY,
		valuesU: uCopy,
		valuesV: vCopy,
	}
}

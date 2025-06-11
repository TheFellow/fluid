package fluid

import "testing"

func TestCopyBorder(t *testing.T) {
	f := New(1.0, 3, 3, 1.0) // creates 5x5 grid due to border
	src := make([]float32, f.numCells)
	for i := range src {
		src[i] = float32(i + 1) // unique value for each cell
	}
	dst := make([]float32, f.numCells)

	f.copyBorder(dst, src)

	n := f.NumY
	for i := 0; i < f.NumX; i++ {
		if dst[i*n+0] != src[i*n+0] {
			t.Errorf("top border mismatch at (%d,0)", i)
		}
		if dst[i*n+f.NumY-1] != src[i*n+f.NumY-1] {
			t.Errorf("bottom border mismatch at (%d,%d)", i, f.NumY-1)
		}
	}
	for j := 0; j < f.NumY; j++ {
		if dst[0*n+j] != src[0*n+j] {
			t.Errorf("left border mismatch at (0,%d)", j)
		}
		if dst[(f.NumX-1)*n+j] != src[(f.NumX-1)*n+j] {
			t.Errorf("right border mismatch at (%d,%d)", f.NumX-1, j)
		}
	}
}

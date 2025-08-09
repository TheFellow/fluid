package fluid

import "testing"

// Test that the vorticity confinement force adds velocity at regions with curl.
func TestApplyVorticityConfinement(t *testing.T) {
	f := New(1, 4, 4, 1)
	for i := 0; i < f.NumX; i++ {
		for j := 0; j < f.NumY; j++ {
			f.SetSolid(i, j, false)
		}
	}
	// Create a small 2x2 vortex in the middle of the grid.
	f.SetVelocity(3, 2, 1, 0)
	f.SetVelocity(2, 3, -1, 0)
	f.SetVelocity(2, 2, 0, 1)
	f.SetVelocity(3, 3, 0, -1)

	f.Confinement = 5
	idx := 2*f.NumY + 2 // Cell inside the vortex
	u0 := f.U[idx]
	v0 := f.V[idx]
	f.applyVorticityConfinement(1)

	if f.U[idx] == u0 && f.V[idx] == v0 {
		t.Fatal("confinement force did not modify velocities")
	}
}

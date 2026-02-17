package fluid

import "testing"

func newBenchFluid() *Fluid {
	f := New(1000, 300, 251, 1.0/100.0)
	for i := 0; i < f.NumX; i++ {
		for j := 0; j < f.NumY; j++ {
			if i == 0 || j == 0 || i == f.NumX-1 || j == f.NumY-1 {
				f.S[i*f.NumY+j] = 0.0
			} else {
				f.S[i*f.NumY+j] = 1.0
			}
		}
	}
	return f
}

// Baseline: ~1ms/frame on Apple M4 Max (300x251 grid)
func BenchmarkSimulate(b *testing.B) {
	f := newBenchFluid()
	dt := float32(1.0 / 120.0)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		f.Simulate(dt)
	}
}

// With jet: ~4.6ms/frame on Apple M4 Max
func BenchmarkSimulateWithJet(b *testing.B) {
	f := newBenchFluid()
	dt := float32(1.0 / 120.0)
	n := f.NumY
	jetVel := float32(4.0)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for j := f.NumY/2 - 50; j < f.NumY/2+50; j++ {
			f.U[1*n+j] = jetVel
			f.M[1*n+j] = 1.0
		}
		f.Simulate(dt)
	}
}

// BFECC with jet: measures 3x advection cost
func BenchmarkSimulateBFECC(b *testing.B) {
	f := newBenchFluid()
	f.UseBFECC = true
	dt := float32(1.0 / 120.0)
	n := f.NumY
	jetVel := float32(4.0)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for j := f.NumY/2 - 50; j < f.NumY/2+50; j++ {
			f.U[1*n+j] = jetVel
			f.M[1*n+j] = 1.0
		}
		f.Simulate(dt)
	}
}

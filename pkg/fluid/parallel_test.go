package fluid

import (
	"fmt"
	"sync/atomic"
	"testing"
)

func TestParallelExecutor(t *testing.T) {
	for _, limit := range []int{1, 4} {
		t.Run(fmt.Sprintf("limit=%d", limit), func(t *testing.T) {
			executor := newParallelExecutor(limit)
			defer executor.close()

			const phases = 3
			counts := make([]atomic.Int32, 257)
			for range phases {
				executor.parallelRange(3, len(counts)-4, func(i int) {
					counts[i].Add(1)
				})
			}

			for i := range counts {
				want := int32(0)
				if i >= 3 && i < len(counts)-4 {
					want = phases
				}
				if got := counts[i].Load(); got != want {
					t.Fatalf("counts[%d] = %d, want %d", i, got, want)
				}
			}
		})
	}
}

func BenchmarkParallelRange(b *testing.B) {
	const rows = 251
	values := make([]uint64, 300*rows)
	work := func(i int) {
		offset := i * rows
		for j := range rows {
			values[offset+j]++
		}
	}

	b.Run("per-call-workers", func(b *testing.B) {
		b.ReportAllocs()
		for range b.N {
			parallelRange(0, len(values)/rows, work)
		}
	})

	b.Run("reused-workers", func(b *testing.B) {
		executor := newParallelExecutor(len(values) / rows)
		defer executor.close()
		b.ReportAllocs()
		b.ResetTimer()
		for range b.N {
			executor.parallelRange(0, len(values)/rows, work)
		}
	})
}

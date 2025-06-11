package fluid

import (
	"runtime"
	"sync"
)

// parallelRange executes fn for each i in [start,end). The range is split among
// available CPUs.
func parallelRange(start, end int, fn func(i int)) {
	total := end - start
	if total <= 0 {
		return
	}
	workers := runtime.GOMAXPROCS(0)
	if workers > total {
		workers = total
	}
	var wg sync.WaitGroup
	chunk := (total + workers - 1) / workers
	for w := 0; w < workers; w++ {
		s := start + w*chunk
		e := s + chunk
		if e > end {
			e = end
		}
		if s >= end {
			break
		}
		wg.Add(1)
		go func(ss, ee int) {
			for i := ss; i < ee; i++ {
				fn(i)
			}
			wg.Done()
		}(s, e)
	}
	wg.Wait()
}

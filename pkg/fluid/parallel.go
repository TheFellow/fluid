package fluid

import (
	"runtime"
	"sync"
)

type parallelJob struct {
	fn         func(int)
	start, end int
}

type parallelExecutor struct {
	jobs        chan parallelJob
	phase       sync.WaitGroup
	workers     sync.WaitGroup
	workerCount int
}

func newParallelExecutor(limit int) *parallelExecutor {
	workerCount := min(runtime.GOMAXPROCS(0), limit)
	executor := &parallelExecutor{workerCount: workerCount}
	if workerCount <= 1 {
		return executor
	}
	executor.jobs = make(chan parallelJob)
	executor.workers.Add(workerCount - 1)
	for range workerCount - 1 {
		go executor.run()
	}
	return executor
}

func (e *parallelExecutor) run() {
	defer e.workers.Done()
	for job := range e.jobs {
		for i := job.start; i < job.end; i++ {
			job.fn(i)
		}
		e.phase.Done()
	}
}

func (e *parallelExecutor) parallelRange(start, end int, fn func(i int)) {
	total := end - start
	if total <= 0 {
		return
	}
	if e.workerCount <= 1 {
		for i := start; i < end; i++ {
			fn(i)
		}
		return
	}
	workers := min(e.workerCount, total)
	chunk := (total + workers - 1) / workers
	background := workers - 1
	e.phase.Add(background)
	for worker := 0; worker < background; worker++ {
		jobStart := start + worker*chunk
		jobEnd := min(jobStart+chunk, end)
		e.jobs <- parallelJob{start: jobStart, end: jobEnd, fn: fn}
	}
	for i := start + background*chunk; i < end; i++ {
		fn(i)
	}
	e.phase.Wait()
}

func (e *parallelExecutor) close() {
	if e.jobs == nil {
		return
	}
	close(e.jobs)
	e.workers.Wait()
}

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

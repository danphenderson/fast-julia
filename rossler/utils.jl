module BenchmarkFlatten

using Statistics
using BenchmarkTools
using DataFrames

"""
    flatten_run_studies(results) -> Vector{NamedTuple}

Flatten the nested dictionary produced by `run_studies()` into a flat vector of
named tuples.  Each entry in the returned vector summarizes a single
benchmark trial and contains the following fields:

* `solver::String` – the descriptive name of the solver (e.g. "RK4 Fixed").
* `variant::String` – the name of the Rössler implementation; keys of the
  `rhs_trials` and `solve_trials` dictionaries.
* `metric::Symbol` – either `:rhs` for micro-benchmarks of the right-hand
  side or `:solve` for full ODE solves.
* `median_time_ns::Float64` – the median run time in nanoseconds, computed
  from the `times` field of the `BenchmarkTools.Trial` object.
* `median_time_s::Float64` – the same median run time converted to seconds.
* `memory::Int` – the number of bytes allocated per call, from the
  `memory` field of the `Trial`.
* `allocs::Int` – the number of allocations per call, from the `allocs`
  field of the `Trial`.

This helper makes it easy to inspect and visualize benchmarking results.  For
example, after running `results = run_studies()`, you can call

```julia
using DataFrames
using .BenchmarkFlatten
flat = flatten_run_studies(results)
df = DataFrame(flat)
first(df, 8)
```

to obtain a tabular view of all recorded benchmarks.
"""
function flatten_run_studies(results)
    rows = NamedTuple[]
    for (solver, study) in results
        rhs_trials = study.rhs_trials
        solve_trials = study.solve_trials
        # Flatten right-hand-side micro-benchmarks
        for (fn, trial) in rhs_trials
            mt = median(trial.times)
            push!(rows, (
                solver = solver,
                variant = fn,
                metric = :rhs,
                median_time_ns = mt,
                median_time_s = mt * 1e-9,
                memory = trial.memory,
                allocs = trial.allocs,
            ))
        end
        # Flatten full solve benchmarks
        for (fn, trial) in solve_trials
            mt = median(trial.times)
            push!(rows, (
                solver = solver,
                variant = fn,
                metric = :solve,
                median_time_ns = mt,
                median_time_s = mt * 1e-9,
                memory = trial.memory,
                allocs = trial.allocs,
            ))
        end
    end
    return rows
end

function dataframe(results)
    flat = flatten_run_studies(results)
    return DataFrame(flat)
end

end # module BenchmarkFlatten
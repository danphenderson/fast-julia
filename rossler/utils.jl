module BenchmarkFlatten

using Statistics
using BenchmarkTools
using DataFrames

@inline function _median_time_ns(trial::BenchmarkTools.Trial)::Float64
    # trial.times is typically Vector{UInt64} in nanoseconds
    return Float64(median(trial.times))
end

"""
    flatten_experiment1(res; label="RK4 Fixed") -> Vector{NamedTuple}

Flatten the result returned by `run_experiment1(...)` into a flat vector of
named tuples suitable for `DataFrame(...)`.

Expected `res` fields:
- `res.dts`
- `res.systems`       (each with `.name` and `.tspan`)
- `res.rhs_trials`    :: Dict{String, BenchmarkTools.Trial}
- `res.solve_trials`  :: Dict{String, Dict{Float64, BenchmarkTools.Trial}}

Returned fields:
* `solver::String`                 – label you choose (default "RK4 Fixed")
* `variant::String`
* `metric::Symbol`                 – `:rhs` or `:solve`
* `dt::Union{Missing,Float64}`     – `missing` for `:rhs`, numeric for `:solve`
* `nsteps::Union{Missing,Int}`     – `missing` for `:rhs`, computed for `:solve`
* `median_time_ns::Float64`
* `median_time_s::Float64`
* `memory::Int`                    – heap bytes per call (BenchmarkTools)
* `allocs::Int`                    – heap allocations per call (BenchmarkTools)
"""
function flatten_experiment1(res; label::AbstractString="RK4 Fixed")
    rows = NamedTuple[]

    rhs_trials   = res.rhs_trials
    solve_trials = res.solve_trials
    systems      = res.systems

    # variant -> tspan for step-count computation
    tspan_by_variant = Dict{String,Tuple{Float64,Float64}}()
    for sys in systems
        tspan_by_variant[String(sys.name)] = sys.tspan
    end

    # RHS micro-benchmarks (dt-independent)
    for (variant, trial) in rhs_trials
        mt = _median_time_ns(trial)
        push!(rows, (
            solver = String(label),
            variant = String(variant),
            metric = :rhs,
            dt = missing,
            nsteps = missing,
            median_time_ns = mt,
            median_time_s = mt * 1e-9,
            memory = Int(trial.memory),
            allocs = Int(trial.allocs),
        ))
    end

    # Full-solve benchmarks (dt sweep)
    for (variant, per_dt) in solve_trials
        v = String(variant)
        tspan = get(tspan_by_variant, v, (NaN, NaN))
        T = tspan[2] - tspan[1]

        for (dt, trial) in per_dt
            mt = _median_time_ns(trial)
            nsteps = isfinite(T) ? Int(ceil(T / dt)) : missing
            push!(rows, (
                solver = String(label),
                variant = v,
                metric = :solve,
                dt = Float64(dt),
                nsteps = nsteps,
                median_time_ns = mt,
                median_time_s = mt * 1e-9,
                memory = Int(trial.memory),
                allocs = Int(trial.allocs),
            ))
        end
    end

    return rows
end

"""
    dataframe(res; label="RK4 Fixed") -> DataFrame
Convenience wrapper around `flatten_experiment1`.
"""
function dataframe(res; label::AbstractString="RK4 Fixed")
    return DataFrame(flatten_experiment1(res; label=label))
end

end # module BenchmarkFlatten

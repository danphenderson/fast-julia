# rossler/benchmark.jl
include(joinpath(@__DIR__, "impl.jl"))
include(joinpath(@__DIR__, "utils.jl"))

using DifferentialEquations
using BenchmarkTools
using StaticArrays
using .BenchmarkFlatten

# ========================================  =====================================
# Experiment 1: RK4 fixed-step, varying RHS implementations, dt-halving sweep
# =============================================================================
@inline ode_benchmark_solve_kwargs() = (;
    adaptive=false,
    save_start=false,
    save_end=false,
    save_everystep=false,
    save_on=false,
    dense=false,
    alias_u0=false,
    timeseries_errors=false,
    maxiters=Int(1e12),
)

# Disallow any other construction (positional and/or keyword arguments).
ODEBenchmarkSolve(args...; kwargs...) =
    throw(ArgumentError("ODEBenchmarkSolve() has fixed defaults; arguments are not supported."))


Base.@kwdef mutable struct ODESys{I,F,U0,P}
    name::String
    f::F
    u0::U0
    p::P
    tspan::Tuple{Float64,Float64}
end

@inline isinplace(::ODESys{I}) where {I} = I

function ode_system(name::AbstractString, f, u0, p, tspan; inplace::Bool=false)
    return ODESys{inplace, typeof(f), typeof(u0), typeof(p)}(
        name=String(name),
        f=f,
        u0=u0,
        p=p,
        tspan=tspan,
    )
end

@inline make_prob(sys::ODESys{I}; tspan::Tuple{Float64,Float64}=sys.tspan) where {I} =
    ODEProblem{I}(sys.f, sys.u0, tspan, sys.p)

# -----------------------------
# Benchmark helpers
# -----------------------------
function bench_rhs(sys::ODESys; samples::Int=200, evals::Int=1)
    f = sys.f
    u = sys.u0
    p = sys.p
    t = sys.tspan[1]

    if isinplace(sys)
        du = similar(u)
        f(du, u, p, t)  # warm-up
        return @benchmark $f($du, $u, $p, $t) samples=samples evals=evals
    else
        f(u, p, t)      # warm-up
        return @benchmark $f($u, $p, $t) samples=samples evals=evals
    end
end

@inline solve_fixed_rk4(prob, dt::Float64) =
    solve(prob, RK4(); dt=dt, ode_benchmark_solve_kwargs()...)

function bench_solve(prob; dt::Float64, samples::Int=20, seconds::Float64=1.0)
    return @benchmark solve_fixed_rk4($prob, $dt) samples=samples seconds=seconds evals=1
end

# -----------------------------
# Experiment 1 specification
# -----------------------------
Base.@kwdef struct Experiment1Spec{U,P}
    u0::U = [1.0, 1.0, 1.0]
    p::P  = [0.1, 0.1, 14.0]
    tspan::Tuple{Float64,Float64} = (0.0, 200.0)

    dt0::Float64 = 1e-2
    halvings::Int = 6               # number of dt values, inclusive of dt0
    warmup_steps::Int = 10          # short warm-up solve length in steps
end

@inline function dt_schedule(dt0::Float64, halvings::Int)
    # [dt0, dt0/2, dt0/4, ...] length == halvings
    return [dt0 / 2.0^(k-1) for k in 1:halvings]
end

function build_systems_for_experiment1(spec::Experiment1Spec)
    u0_in = spec.u0
    p_in  = spec.p
    tspan = spec.tspan

    # Ensure vector-based variants have a mutable Vector u0 (important for in-place RHS).
    u0_vec = u0_in isa SVector ? collect(u0_in) :
             (u0_in isa AbstractVector ? u0_in : collect(u0_in))

    p_vec  = p_in isa AbstractVector ? p_in : collect(p_in)

    # Static versions use SVectors for state; params can remain any indexable container.
    s_u0 = u0_in isa SVector ? u0_in : SVector(u0_vec...)
    s_p  = p_in  isa SVector ? p_in  : SVector(p_vec...)

    systems = ODESys[]

    push!(systems, ode_system("rossler_naive", rossler_naive, u0_vec, p_vec, tspan; inplace=false))
    push!(systems, ode_system("rossler", rossler, u0_vec, p_vec, tspan; inplace=false))

    # Use fresh u0 vectors for in-place variants (defensive; avoids any accidental aliasing surprises).
    push!(systems, ode_system("rossler_naive!", rossler_naive!, copy(u0_vec), p_vec, tspan; inplace=true))
    push!(systems, ode_system("rossler!", rossler!, copy(u0_vec), p_vec, tspan; inplace=true))

    push!(systems, ode_system("rossler_static_naive", rossler_static_naive, s_u0, p_vec, tspan; inplace=false))
    push!(systems, ode_system("rossler_static", rossler_static, s_u0, p_vec, tspan; inplace=false))

    push!(systems, ode_system("rossler_type_stable", rossler_type_stable, u0_vec, p_vec, tspan; inplace=false))

    # Optional: include AD-ready variant if it exists in impl.jl
    if isdefined(@__MODULE__, :rossler_ad)
        push!(systems, ode_system("rossler_ad", rossler_ad, s_u0, s_p, tspan; inplace=false))
    end

    return systems
end

"""
    run_experiment1(spec; rhs_samples=200, rhs_evals=1, solve_samples=20, solve_seconds=1.0)

Returns (NamedTuple):
- dts          :: Vector{Float64}
- systems      :: Vector{ODESys}  (one per RHS variant)
- rhs_trials   :: Dict{String, BenchmarkTools.Trial}
- solve_trials :: Dict{String, Dict{Float64, BenchmarkTools.Trial}}

Notes:
- RHS trials are collected once per variant (dt-independent).
- Solve trials are collected per (variant, dt) for dt in the halving schedule.
- No attempt is made to “measure stack allocations” numerically; instead, heap allocs/bytes
  are reported (0 heap allocs implies stack/register-resident temporaries in practice).
"""
function run_experiment1(spec::Experiment1Spec;
    rhs_samples::Int=200,
    rhs_evals::Int=1,
    solve_samples::Int=20,
    solve_seconds::Float64=1.0,
)
    dts = dt_schedule(spec.dt0, spec.halvings)
    systems = build_systems_for_experiment1(spec)

    # Warm-up compilation with a short tspan (avoid paying full-horizon cost just to compile).
    for sys in systems
        t0 = sys.tspan[1]
        tw = (t0, t0 + spec.warmup_steps * dts[1])
        prob_warm = make_prob(sys; tspan=tw)
        solve_fixed_rk4(prob_warm, dts[1])
    end

    rhs_trials = Dict{String,BenchmarkTools.Trial}()
    for sys in systems
        rhs_trials[sys.name] = bench_rhs(sys; samples=rhs_samples, evals=rhs_evals)
    end

    solve_trials = Dict{String, Dict{Float64,BenchmarkTools.Trial}}()
    for sys in systems
        prob = make_prob(sys)  # fixed model + horizon; vary dt only
        per_dt = Dict{Float64,BenchmarkTools.Trial}()
        for dt in dts
            per_dt[dt] = bench_solve(prob; dt=dt, samples=solve_samples, seconds=solve_seconds)
        end
        solve_trials[sys.name] = per_dt
    end

    return (dts=dts, systems=systems, rhs_trials=rhs_trials, solve_trials=solve_trials)
end

function write_results_to_csv(res)
    dts = res.dts
    rhs_trials = res.rhs_trials
    solve_trials = res.solve_trials

    # Write RHS results
    open("rhs_benchmarks.csv", "w") do io
        println(io, "RHS,Benchmark")
        for (name, trial) in rhs_trials
            println(io, "$name,$trial")
        end
    end

    # Write solve results
    open("solve_benchmarks.csv", "w") do io
        println(io, "System,dt,Benchmark")
        for (name, per_dt) in solve_trials
            for (dt, trial) in per_dt
                println(io, "$name,$dt,$trial")
            end
        end
    end
end

function main()
    spec = Experiment1Spec()
    res = run_experiment1(spec)
    println("writing results to CSV...")
    path = BenchmarkFlatten.write_results_to_csv(res; outpath="poster/results.csv", solver_label="RK4 Fixed")
    println("wrote: ", path)
    println("RHS benchmarks:")
    for (name, trial) in res.rhs_trials
        println("$name: ", trial)
    end

    println("\nSolve benchmarks:")
    for (name, per_dt) in res.solve_trials
        println("System: $name")
        for (dt, trial) in per_dt
            println("  dt=$dt: ", trial)
        end
    end

    println("writing results to CSV...")
    write_results_to_csv(res)
end
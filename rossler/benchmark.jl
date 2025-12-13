# rossler/benchmark.jl
"""Module includes rossler implementations (`rossler/impl.jl`)
and related utilities, and herein, the ODE benchmarking pipeline is implemented.
(satisfing the SciML ODESystem interface)
"""

# Bring in the Rössler right–hand side implementations.  This include defines
# `rossler`, `rossler!`, `rossler_static`, `rossler_type_stable` and `rossler_ad`.
include(joinpath(@__DIR__, "impl.jl"))


using DifferentialEquations
using SciMLBase   # for AbstractODESolution (more stable than DESolution)
using BenchmarkTools
using LinearAlgebra
using StaticArrays

# ===========================================================================
# BEGIN defining types
# ===========================================================================
Base.@kwdef mutable struct ODESys{I,F,U0,P,A}
    f::F
    u0::U0
    p::P
    tspan::Tuple{Float64,Float64}
    alg::A = Tsit5()
    dt::Float64 = 0.01
end

function ODESys{I}(; f, u0, p, tspan, alg=Tsit5(), dt=0.01) where {I}
    return ODESys{I, typeof(f), typeof(u0), typeof(p), typeof(alg)}(f, u0, p, tspan, alg, dt)
end

@inline isinplace(::ODESys{I}) where {I} = I

ode_system(f, u0, p, tspan; inplace::Bool=false, alg=Tsit5(), dt=0.01) =
    ODESys{inplace}(; f, u0, p, tspan, alg, dt)

@inline make_prob(sys::ODESys{I}) where {I} =
    ODEProblem{I}(sys.f, sys.u0, sys.tspan, sys.p)


Base.@kwdef struct ODESolveParams{A}
    alg::A = Tsit5()

    adaptive::Bool = false
    dt::Union{Nothing,Float64} = nothing  # if nothing, use sys.dt

    # adaptive-only
    reltol::Float64 = 1e-3
    abstol::Float64 = 1e-6

    # output control
    saveat::Union{Nothing,Float64,AbstractVector{Float64}} = nothing
    save_start::Bool = true
    save_end::Bool = true
    save_everystep::Bool = false
    save_on::Bool = true
    dense::Bool = true

    # progress
    progress::Bool = false
    progress_steps::Int = 1000
    progress_name::String = "ODE Solve"
    progress_message::Union{Nothing,String} = nothing
    progress_id::Union{Nothing,Symbol} = nothing

    alias_u0::Bool = false
    timeseries_errors::Bool = false
end

function solve_system(sys::ODESys, sp::ODESolveParams; kwargs...)
    prob = make_prob(sys)
    dt = sp.dt === nothing ? sys.dt : sp.dt

    # Start with always-safe keywords
    kw = (
        save_start=sp.save_start,
        save_end=sp.save_end,
        save_everystep=sp.save_everystep,
        save_on=sp.save_on,
        dense=sp.dense,
        alias_u0=sp.alias_u0,
        timeseries_errors=sp.timeseries_errors,
    )

    # Only pass when set (DiffEq does not like saveat=nothing)
    if sp.saveat !== nothing
        kw = merge(kw, (saveat=sp.saveat,))
    end

    # Adaptive vs fixed-step
    if sp.adaptive
        kw = merge(kw, (adaptive=true, reltol=sp.reltol, abstol=sp.abstol))
    else
        kw = merge(kw, (adaptive=false, dt=dt))
    end

    # Progress options: only pass when progress is enabled
    if sp.progress
        kw = merge(kw, (
            progress=true,
            progress_steps=sp.progress_steps,
            progress_name=sp.progress_name,
        ))
        if sp.progress_message !== nothing
            kw = merge(kw, (progress_message=sp.progress_message,))
        end
        if sp.progress_id !== nothing
            kw = merge(kw, (progress_id=sp.progress_id,))  # Symbol
        end
    end

    # Let caller override any of the above without duplicate keywords
    user_kw = (; kwargs...)
    kw = merge(kw, user_kw)
    return solve(prob, sp.alg; kw...)
end

################################################################################
# Fast-path solve without dynamic NamedTuple merging
################################################################################
"""
    solve_system_fast(sys, sp)

Solve the ODE problem associated with `sys` under the solve parameters `sp`
without allocating intermediate keyword arguments via `NamedTuple` merging.
This is intended for benchmarking loops where dynamic keyword merging can
become a bottleneck.  The arguments are passed directly to `solve` with a
fixed set of keywords.
"""
function solve_system_fast(sys::ODESys, sp::ODESolveParams)
    prob = make_prob(sys)
    # Determine the effective timestep (for non-adaptive solves)
    dt = sp.dt === nothing ? sys.dt : sp.dt
    if sp.adaptive
        return solve(prob, sp.alg;
            adaptive=true,
            reltol=sp.reltol,
            abstol=sp.abstol,
            saveat=sp.saveat,
            save_start=sp.save_start,
            save_end=sp.save_end,
            save_everystep=sp.save_everystep,
            save_on=sp.save_on,
            dense=sp.dense,
            alias_u0=sp.alias_u0,
            timeseries_errors=sp.timeseries_errors,
        )
    else
        return solve(prob, sp.alg;
            adaptive=false,
            dt=dt,
            saveat=sp.saveat,
            save_start=sp.save_start,
            save_end=sp.save_end,
            save_everystep=sp.save_everystep,
            save_on=sp.save_on,
            dense=sp.dense,
            alias_u0=sp.alias_u0,
            timeseries_errors=sp.timeseries_errors,
        )
    end
end

################################################################################
# Benchmark helpers
################################################################################
"""
    bench_rhs(sys; samples=100, evals=1)

Return a `BenchmarkTools.Trial` object representing a micro-benchmark of a
single right–hand side call for the ODE system `sys`.  This benchmark warms
up compilation automatically and measures allocations and runtime for a
single call to `sys.f(sys.u0, sys.p, sys.tspan[1])`.
"""
function bench_rhs(sys::ODESys; samples::Int=100, evals::Int=1)
    f = sys.f
    u = sys.u0
    p = sys.p
    t = sys.tspan[1]
    # Warm-up call outside the timed region
    f(u, p, t)
    return @benchmark $f($u, $p, $t) samples=samples evals=evals
end

"""
    bench_solve(sys, sp; samples=10, seconds=1.0)

Return a `BenchmarkTools.Trial` object representing a benchmark of solving
the ODE system `sys` with parameters `sp` using the fast-path solver.
The first solve is run outside of the measurement to trigger compilation.
`samples` and `seconds` are passed through to the `@benchmark` macro.
"""
function bench_solve(sys::ODESys, sp::ODESolveParams; samples::Int=10, seconds=1.0)
    # Warm-up to compile everything
    solve_system_fast(sys, sp)
    return @benchmark solve_system_fast($sys, $sp) samples=samples seconds=seconds evals=1
end

################################################################################
# Case study 1 experiment spec and runner
################################################################################
"""
    CaseStudy1Spec

Container for specifying a Case Study 1 experiment.  It holds the initial
state, parameters, time span, timestep, solver algorithm and other optional
fields used to instantiate each Rössler variant.  Instances of this type
are passed to `run_case_study_1` to perform the simulation, benchmarking and
result collection.
"""
Base.@kwdef struct CaseStudy1Spec{U,P}
    """Initial state vector or static vector."""
    u0::U
    """Parameter vector or static vector."""
    p::P
    """Integration time span."""
    tspan::Tuple{Float64,Float64}
    """Time step for fixed-step solves."""
    dt::Float64
    """ODE solver algorithm."""
    alg::Any = Tsit5()
    """Save grid (optional) passed to the solver."""
    saveat::Union{Nothing,Float64,AbstractVector{Float64}} = nothing
    """Relative tolerance for adaptive solves."""
    reltol::Float64 = 1e-3
    """Absolute tolerance for adaptive solves."""
    abstol::Float64 = 1e-6
end

"""
    run_case_study_1(spec) -> NamedTuple

Run the Case Study 1 experiment described by `spec`.  This function constructs
all supported Rössler variants, performs warm-up solves to trigger JIT
compilation, benchmarks the right–hand side and full solves using the
fast-path solver, and returns a named tuple containing the systems, the
micro-benchmark trials, and the full solve benchmark trials.
"""
function run_case_study_1(spec::CaseStudy1Spec)
    # Build variant systems.  `u0` and `p` may be vectors or static vectors.
    u0 = spec.u0
    p  = spec.p
    tspan = spec.tspan
    dt = spec.dt
    alg = spec.alg
    # Ensure static variants are correctly typed
    s_u0 = u0 isa SVector ? u0 : SVector(u0...)
    s_p  = p  isa SVector ? p  : SVector(p...)
    systems = [
        # Standard out-of-place implementation
        ode_system(rossler,        u0,      p,   tspan; inplace=false, alg=alg, dt=dt),
        # In-place implementation (state must be copied to avoid aliasing)
        ode_system(rossler!,       copy(u0), p,   tspan; inplace=true,  alg=alg, dt=dt),
        # Static-array implementation (stack allocation)
        ode_system(rossler_static, s_u0,    p,   tspan; inplace=false, alg=alg, dt=dt),
        # Type-stable out-of-place implementation
        ode_system(rossler_type_stable, u0, p,   tspan; inplace=false, alg=alg, dt=dt),
        # AD-ready allocation-free implementation
        ode_system(rossler_ad,     s_u0,    s_p, tspan; inplace=false, alg=alg, dt=dt),
    ]
    # Warm-up solves for compilation
    for sys in systems
        sp = ODESolveParams(; alg=alg, adaptive=true, reltol=spec.reltol, abstol=spec.abstol, dt=spec.dt, save_on=false, dense=false, save_everystep=false)
        solve_system_fast(sys, sp)
    end
    # Benchmark right-hand sides
    rhs_trials = Dict{String,BenchmarkTools.Trial}()
    for sys in systems
        trial = bench_rhs(sys; samples=100, evals=1)
        rhs_trials[string(sys.f)] = trial
    end
    # Benchmark full solves with recommended benchmarking parameters
    solve_trials = Dict{String,BenchmarkTools.Trial}()
    for sys in systems
        sp = benchmark_params(sys)
        trial = bench_solve(sys, sp; samples=10, seconds=1.0)
        solve_trials[string(sys.f)] = trial
    end
    return (systems=systems, rhs_trials=rhs_trials, solve_trials=solve_trials)
end


# Parameter presets
test_params(sys::ODESys) =
    ODESolveParams(; alg=sys.alg, adaptive=false, dt=sys.dt,
                   timeseries_errors=true, progress=true,
                   save_everystep=true, dense=true)

benchmark_params(sys::ODESys) =
    ODESolveParams(; alg=sys.alg, adaptive=true,
                   dense=false, save_everystep=false, save_on=false)

interp_params(sys::ODESys) =
    ODESolveParams(; alg=sys.alg, adaptive=true,
                   dense=true, save_everystep=false)

gif_params(sys::ODESys) =
    ODESolveParams(; alg=sys.alg, adaptive=true,
                   saveat=5*sys.dt, dense=false, save_everystep=false)


# Benchmark result container
Base.@kwdef struct ODEBenchmarkResult
    func_name::String
    language::String
    variant::String
    solver::String
    dt::Float64
    reltol::Float64
    abstol::Float64
    tspan::Tuple{Float64,Float64}
    sol::Union{Nothing,SciMLBase.AbstractODESolution} = nothing
end

function Base.show(io::IO, r::ODEBenchmarkResult)
    print(io,
        "ODEBenchmarkResult($(r.func_name), solver=$(r.solver), dt=$(r.dt), ",
        "reltol=$(r.reltol), abstol=$(r.abstol), tspan=$(r.tspan))"
    )
end
endpretty_print(r::ODEBenchmarkResult) = begin
    println("ODEBenchmarkResult:")
    println("  Function: $(r.func_name)")
    println("  Language: $(r.language)")
    println("  Variant: $(r.variant)")
    println("  Solver: $(r.solver)")
    println("  dt: $(r.dt)")
    println("  reltol: $(r.reltol)")
    println("  abstol: $(r.abstol)")
    println("  tspan: $(r.tspan)")
end
# ===========================================================================
# END defining types
# =========================================================================== 


# ===========================================================================
# BEGIN benchmarking pipeline
# ===========================================================================
function benchmark_setup()
    BLAS.set_num_threads(1)
    return [
        ode_system(rossler,            [1.0,1.0,1.0],              [0.1,0.1,14.0],   (0.0,200.0); inplace=false),
        ode_system(rossler!,           [1.0,1.0,1.0],              [0.1,0.1,14.0],   (0.0,200.0); inplace=true),
        ode_system(rossler_static,     SVector(1.0,1.0,1.0),       [0.1,0.1,14.0],   (0.0,200.0); inplace=false),
        ode_system(rossler_type_stable,[1.0,1.0,1.0],              [0.1,0.1,14.0],   (0.0,200.0); inplace=false),
        ode_system(rossler_ad,         SVector(1.0,1.0,1.0),       SVector(0.1,0.1,14.0), (0.0,200.0); inplace=false),
    ]
end


function benchmark_shutdown()
    BLAS.set_num_threads(4)
end

function benchmark_pipeline(sys::ODESys, params::ODESolveParams)
    sol = solve_system(sys, params)
    # Determine a human-readable variant name based on the function used
    variant = begin
        if sys.f === rossler
            "standard"
        elseif sys.f === rossler!
            "inplace"
        elseif sys.f === rossler_static
            "static"
        elseif sys.f === rossler_type_stable
            "type-stable"
        elseif sys.f === rossler_ad
            "ad-ready"
        else
            isinplace(sys) ? "inplace" : "non-inplace"
        end
    end
    return ODEBenchmarkResult(
        func_name="rossler",
        language="julia",
        variant=variant,
        solver=string(params.alg),
        dt=params.dt === nothing ? sys.dt : params.dt,
        reltol=params.reltol,
        abstol=params.abstol,
        tspan=sys.tspan,
        sol=sol
    )
end

function benchmark_pipelines()
    systems = benchmark_setup()
    params = benchmark_params(first(systems))  # or compute per-system if you prefer
    results = [benchmark_pipeline(s, params) for s in systems]
    benchmark_shutdown()
    return results
end 
# ===========================================================================
# END benchmarking pipeline
# =========================================================================== 


# rossler/benchmark.jl
"""Module includes rossler implementations (`rossler/impl.jl`)
and related utilities, and herein, the ODE benchmarking pipeline is implemented.
(satisfing the SciML ODESystem interface)
"""

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
        ode_system(rossler,        [1.0,1.0,1.0],        [0.1,0.1,14.0], (0.0,200.0); inplace=false),
        ode_system(rossler_static, SVector(1.0,1.0,1.0), [0.1,0.1,14.0], (0.0,200.0); inplace=false),
        ode_system(rossler!,       [1.0,1.0,1.0],        [0.1,0.1,14.0], (0.0,200.0); inplace=true),
    ]
end


function benchmark_shutdown()
    BLAS.set_num_threads(4)
end

function benchmark_pipeline(sys::ODESys, params::ODESolveParams)
    sol = solve_system(sys, params)
    return ODEBenchmarkResult(
        func_name="rossler",
        language="julia",
        variant=isinplace(sys) ? "inplace" : "non-inplace",
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


# rossler/benchmark.jl
"""Module includes rossler implementations (`rossler/impl.jl`)
and related utilities, and herein, the ODE benchmarking pipeline is implemented.
(satisfing the SciML ODESystem interface)
"""

include("impl.jl")

using UUIDs
using DifferentialEquations
using SciMLBase   # for AbstractODESolution (more stable than DESolution)
using BenchmarkTools
using LinearAlgebra
using StaticArrays

# ===========================================================================
# BEGIN defining types
# ===========================================================================
Base.@kwdef mutable struct ODESystem{I,F,U0,P,A}
    f::F
    u0::U0
    p::P
    tspan::Tuple{Float64,Float64}
    alg::A = Tsit5()
    dt::Float64 = 0.01
end

ode_system(f, u0, p, tspan; inplace::Bool=false, alg=Tsit5(), dt=0.01) =
    ODESystem{inplace}(f=f, u0=u0, p=p, tspan=tspan, alg=alg, dt=dt)

@inline make_prob(sys::ODESystem{true})  = ODESystem{true}(sys.f, sys.u0, sys.tspan, sys.p)
@inline make_prob(sys::ODESystem{false}) = ODESystem{false}(sys.f, sys.u0, sys.tspan, sys.p)

ode_system(f::Function, vu0::Vector{Float64}, p::Vector{Float64}, tspan::Tuple{Float64,Float64}) =
    ODESystem(f=f, vu0=vu0, p=p, tspan=tspan)
ode_system(f::Function, vx0::SVector{3,Float64}, p::Vector{Float64}, tspan::Tuple{Float64,Float64}) =
    ODESystem(f=f, vx0=vx0, p=p, tspan=tspan, inplace=false, u0_kind=:static)

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
    progress_message = nothing
    progress_id::UUID = uuid4()

    alias_u0::Union{Bool,Nothing} = false
    timeseries_errors::Bool = false
end

# Helpers: choose u0 + build ODESystem with explicit in-placeness
@inline function u0(sys::ODESystem)
    if sys.u0_kind === :vector
        sys.vu0 === nothing && error("u0_kind=:vector but vu0 is not set")
        return sys.vu0
    elseif sys.u0_kind === :static
        sys.vx0 === nothing && error("u0_kind=:static but vx0 is not set")
        return sys.vx0
    else
        error("Unknown u0_kind=$(sys.u0_kind). Expected :vector or :static.")
    end
end

function make_prob(sys::ODESystem)
    u0_ = u0(sys)
    if sys.inplace
        return ODESystem{true}(sys.f, u0_, sys.tspan, sys.p)
    else
        return ODESystem{false}(sys.f, u0_, sys.tspan, sys.p)
    end
end

function solve_system(sys::ODESystem, sp::ODESolveParams; kwargs...)
    prob = make_prob(sys)
    dt = sp.dt === nothing ? sys.dt : sp.dt

    if sp.adaptive
        return solve(prob, sp.alg;
            adaptive=true, reltol=sp.reltol, abstol=sp.abstol,
            saveat=sp.saveat, save_start=sp.save_start, save_end=sp.save_end,
            save_everystep=sp.save_everystep, save_on=sp.save_on, dense=sp.dense,
            alias_u0=sp.alias_u0, timeseries_errors=sp.timeseries_errors,
            progress=sp.progress, progress_steps=sp.progress_steps,
            progress_name=sp.progress_name, progress_message=sp.progress_message,
            progress_id=sp.progress_id,
            kwargs...)
    else
        return solve(prob, sp.alg;
            adaptive=false, dt=dt,
            saveat=sp.saveat, save_start=sp.save_start, save_end=sp.save_end,
            save_everystep=sp.save_everystep, save_on=sp.save_on, dense=sp.dense,
            alias_u0=sp.alias_u0, timeseries_errors=sp.timeseries_errors,
            progress=sp.progress, progress_steps=sp.progress_steps,
            progress_name=sp.progress_name, progress_message=sp.progress_message,
            progress_id=sp.progress_id,
            kwargs...)
    end
end

# Parameter presets
test_params(sys::ODESystem) =
    ODESolveParams(; alg=sys.alg, adaptive=false, dt=sys.dt,
                   timeseries_errors=true, progress=true,
                   save_everystep=true, dense=true)

benchmark_params(sys::ODESystem) =
    ODESolveParams(; alg=sys.alg, adaptive=true,
                   dense=false, save_everystep=false, save_on=false)

interp_params(sys::ODESystem) =
    ODESolveParams(; alg=sys.alg, adaptive=true,
                   dense=true, save_everystep=false)

gif_params(sys::ODESystem) =
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
# ===========================================================================
# END defining types
# =========================================================================== 


# ===========================================================================
# BEGIN benchmarking pipeline
# ===========================================================================


function benchmark_setup()
    # set environment variables
    ENV["JULIA_NUM_THREADS"] = "1"
    ENV["JULIA_CPU_THREADS"] = "1"
    BLAS.set_num_threads(1)
    return [
        ode_system(rossler,  [1.0,1.0,1.0], [0.1,0.1,14.0], (0.0,200.0); inplace=false) # non-inplace, vector
        ode_system(rossler,  SVector(1.0,1.0,1.0),  [0.1,0.1,14.0], (0.0,200.0); inplace=false) # non-inplace, static
        ode_system(rossler!, [1.0,1.0,1.0],         [0.1,0.1,14.0], (0.0,200.0); inplace=true) # inplace, vector
        # ode_system(rossler!, SVector(1.0,1.0,1.0),  [0.1,0.1,14.0], (0.0,200.0); inplace=true) # inplace, static
    ]
end

function benchmark_shutdown()
    # reset environment variables
    ENV["JULIA_NUM_THREADS"] = "auto"
    ENV["JULIA_CPU_THREADS"] = "auto"
    BLAS.set_num_threads(4)
end

function benchmark_pipeline(sys::ODESystem, params::ODESolveParams)
    sol = solve_system(sys, params)
    result = ODEBenchmarkResult(
        func_name="rossler",
        language="julia",
        variant=sys.inplace ? "inplace" : "non-inplace",
        solver=string(params.alg),
        dt=params.dt === nothing ? sys.dt : params.dt,
        reltol=params.reltol,
        abstol=params.abstol,
        tspan=sys.tspan,
        sol=sol
    )
    return result
end

function benchmark_pipelines()
    sys = benchmark_setup()
    params = benchmark_params(sys)
    results = [benchmark_pipeline(s, params) for s in sys]
    benchmark_shutdown()
    return results
end
# ===========================================================================
# END benchmarking pipeline
# =========================================================================== 


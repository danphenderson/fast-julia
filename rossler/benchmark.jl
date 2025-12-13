# rossler/benchmark.jl
"""Module includes rossler implementations (`rossler/impl.jl`)
and related utilities, and herein, the ODE benchmarking pipeline is implemented.
(satisfing the SciML ODEProblem interface)
"""

include("rossler.jl")

using UUIDs
using DifferentialEquations
using SciMLBase   # for AbstractODESolution (more stable than DESolution)

# ===========================================================================
# BEGIN defining types
# ===========================================================================
Base.@kwdef mutable struct ODESystem{F,Uv,Ux,P,A}
    f::F
    vu0::Uv
    vx0::Ux
    p::P
    tspan::Tuple{Float64,Float64}
    alg::A = Tsit5()
    dt::Float64 = 0.01

    # Implementation metadata (needed to pick correct u0 + problem form)
    inplace::Bool = false             # true => f!(du,u,p,t), false => du = f(u,p,t)
    u0_kind::Symbol = :vector         # :vector => vu0, :static => vx0
end

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

# Helpers: choose u0 + build ODEProblem with explicit in-placeness
u0(sys::ODESystem) = sys.u0_kind === :static ? sys.vx0 : sys.vu0

function make_prob(sys::ODESystem)
    u0_ = u0(sys)
    return sys.inplace ?
        ODEProblem{true}(sys.f, u0_, sys.tspan, sys.p) :
        ODEProblem{false}(sys.f, u0_, sys.tspan, sys.p)
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
function _benchmark_setup()
    # set environment variables
    ENV["JULIA_NUM_THREADS"] = "1"
    ENV["JULIA_CPU_THREADS"] = "1"
    BLAS.set_num_threads(1)
end

function _ode_solve(sys::ODESystem, params::ODESolveParams)
end

function _test_ode_solve()
end

function _benchmark_ode_solve()
end

function _benchmark_close()
end

function benchmark_pipeline(sys::ODESystem, params::ODESolveParams)
    @benchmark bench_result = begin
        _benchmark_setup()
        sol = _ode_solve(sys, params)
        _test_ode_solve(sol)
        _benchmark_ode_solve(sol)
        _benchmark_close()
    end
    return ODEBenchmarkResult(
        func_name=sys.f.name,
        language="Julia",
        variant="SciML",
        solver=params.alg.name,
        dt=params.dt,
        reltol=params.reltol,
        abstol=params.abstol,
        tspan=sys.tspan,
        sol=sol
    )
end

function benchmark_solver(sys::ODESystem, params::ODESolveParams, algs::Vector{<:DEAlgorithm})
    results = []
    for alg in algs
        params.alg = alg
        result = benchmark_pipeline(sys, params)
        push!(results, result)
    end
    return results
end

function benchmark_all_solvers(sys::ODESystem, params::ODESolveParams)
    algs = [Tsit5(), Tsit5(), Vern7(), Vern9()]
    return benchmark_solver(sys, params, algs)
end
# ===========================================================================
# END benchmarking pipeline
# =========================================================================== 


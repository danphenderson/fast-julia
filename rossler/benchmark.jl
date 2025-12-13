# rossler/benchmark.jl
"""Module includes rossler implementations (`rossler/impl.jl`)
and related utilities, and herein, the ODE benchmarking pipeline is implemented.
(satisfing the SciML ODESystem interface)
"""

include(joinpath(@__DIR__, "impl.jl"))

using DifferentialEquations
using BenchmarkTools
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

    # keep these for extensibility; solve_system_fast passes them through
    alias_u0::Bool = false
    timeseries_errors::Bool = false
end

################################################################################
# Fast-path solve without dynamic NamedTuple merging
################################################################################
function solve_system_fast(sys::ODESys, sp::ODESolveParams)
    prob = make_prob(sys)
    dt = sp.dt === nothing ? sys.dt : sp.dt

    if sp.adaptive
        if sp.saveat === nothing
            return solve(prob, sp.alg;
                adaptive=true,
                reltol=sp.reltol,
                abstol=sp.abstol,
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
        end
    else
        if sp.saveat === nothing
            return solve(prob, sp.alg;
                adaptive=false,
                dt=dt,
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
end

################################################################################
# Benchmark helpers (used by run_case_study/run_studies)
################################################################################
function bench_rhs(sys::ODESys; samples::Int=100, evals::Int=1)
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

function bench_solve(sys::ODESys, sp::ODESolveParams; samples::Int=10, seconds=1.0)
    solve_system_fast(sys, sp) # warm-up
    return @benchmark solve_system_fast($sys, $sp) samples=samples seconds=seconds evals=1
end

################################################################################
# Case study runner (what run_studies() uses)
################################################################################
Base.@kwdef struct CaseStudySpec{U,P}
    u0::U
    p::P
    tspan::Tuple{Float64,Float64}
    dt::Float64
    alg::Any = RK4()
    adaptive::Bool = false
    saveat::Union{Nothing,Float64,AbstractVector{Float64}} = nothing
    reltol::Float64 = 1e-3
    abstol::Float64 = 1e-6
end

case_study_benchmark_spec(alg; adaptive=false) =
    CaseStudySpec(u0=[1.0,1.0,1.0], p=[0.1,0.1,14.0],
                  tspan=(0.0,2000.0), dt=1e-2, alg=alg, adaptive=adaptive)

function run_case_study(spec::CaseStudySpec)
    u0 = spec.u0
    p  = spec.p
    tspan = spec.tspan
    dt = spec.dt
    alg = spec.alg

    s_u0 = u0 isa SVector ? u0 : SVector(u0...)
    s_p  = p  isa SVector ? p  : SVector(p...)

    systems = [
        ode_system(rossler_naive,         u0,       p,   tspan; inplace=false, alg=alg, dt=dt),
        ode_system(rossler,               u0,       p,   tspan; inplace=false, alg=alg, dt=dt),
        ode_system(rossler_naive!,    copy(u0),     p,   tspan; inplace=true,  alg=alg, dt=dt),
        ode_system(rossler!,          copy(u0),     p,   tspan; inplace=true,  alg=alg, dt=dt),
        ode_system(rossler_static_naive,  s_u0,     p,   tspan; inplace=false, alg=alg, dt=dt),
        ode_system(rossler_static,        s_u0,     p,   tspan; inplace=false, alg=alg, dt=dt),
        ode_system(rossler_type_stable,   u0,       p,   tspan; inplace=false, alg=alg, dt=dt),
        ode_system(rossler_ad,            s_u0,     s_p, tspan; inplace=false, alg=alg, dt=dt),
    ]

    # Warm-up solves for compilation
    for sys in systems
        sp = spec.adaptive ?
            ODESolveParams(; alg=alg, adaptive=true, reltol=spec.reltol, abstol=spec.abstol,
                            save_on=false, dense=false, save_everystep=false) :
            ODESolveParams(; alg=alg, adaptive=false, dt=spec.dt,
                            save_on=false, dense=false, save_everystep=false)

        solve_system_fast(sys, sp)
    end

    rhs_trials = Dict{String,BenchmarkTools.Trial}()
    for sys in systems
        rhs_trials[string(sys.f)] = bench_rhs(sys; samples=100, evals=1)
    end

    solve_trials = Dict{String,BenchmarkTools.Trial}()
    for sys in systems
        sp = spec.adaptive ?
            ODESolveParams(; alg=alg, adaptive=true, reltol=spec.reltol, abstol=spec.abstol,
                            dense=false, save_everystep=false, save_on=false) :
            ODESolveParams(; alg=alg, adaptive=false, dt=spec.dt,
                            dense=false, save_everystep=false, save_on=false)

        solve_trials[string(sys.f)] = bench_solve(sys, sp; samples=1000, seconds=1.0)
    end

    return (systems=systems, rhs_trials=rhs_trials, solve_trials=solve_trials)
end

function run_studies()
    rk4_fixed_results        = run_case_study(case_study_benchmark_spec(RK4();  adaptive=false))
    tsit5_adaptive_results   = run_case_study(case_study_benchmark_spec(Tsit5(); adaptive=true))
    euler_fixed_results      = run_case_study(case_study_benchmark_spec(Euler(); adaptive=false))
    midpoint_fixed_results   = run_case_study(case_study_benchmark_spec(Midpoint(); adaptive=false))
    return Dict(
        "RK4 Fixed"        => rk4_fixed_results,
        "Tsit5 Adaptive"   => tsit5_adaptive_results,
        "Euler Fixed"      => euler_fixed_results,
        "Midpoint Fixed"   => midpoint_fixed_results,
    )
end
# ===========================================================================
# END defining types
# ===========================================================================
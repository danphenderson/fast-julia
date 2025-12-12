# benchmarking
using BenchmarkTools
using DifferentialEquations
using StaticArrays
using Plots
using Statistics

# ----------------------------
# Rossler RHS implementations (canonical)
# x' = -y - z
# y' =  x + a y
# z' =  b + z(x - c)
# ----------------------------

# Naive out-of-place: allocates a new Vector each call.
function rossler(vx, vp, t)
    dx1 = -vx[2] - vx[3]
    dx2 =  vx[1] + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (vx[1] - vp[3])
    return [dx1, dx2, dx3]
end

# Out-of-place with type annotations (AD-friendly: t unconstrained).
function rossler_annotated(vx::AbstractVector{T}, vp::AbstractVector{T}, t) where {T<:Real}
    dx1 = -vx[2] - vx[3]
    dx2 =  vx[1] + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (vx[1] - vp[3])
    return T[dx1, dx2, dx3]
end

# In-place: writes into dx.
function rossler!(dx, vx, vp, t)
    x1 = vx[1]
    dx[1] = -vx[2] - vx[3]
    dx[2] =  x1 + vp[1] * vx[2]
    dx[3] =  vp[2] + vx[3] * (x1 - vp[3])
    return nothing
end

# Static out-of-place: returns SVector (often fastest for tiny systems).
function rossler_static(vx, vp, t)
    x1 = vx[1]
    dx1 = -vx[2] - vx[3]
    dx2 =  x1 + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (x1 - vp[3])
    return @SVector [dx1, dx2, dx3]
end

# ----------------------------
# Benchmark harness helpers
# ----------------------------
Base.@kwdef mutable struct BenchmarkDriver{S}
    tspan::Tuple{Float64,Float64} = (0.0, 50.0)
    saveat::S = 0.0:0.1:50.0
    dt::Float64 = 1e-4
    maxiters::Int = 10^9
    save_everystep::Bool = false
    abs_tol::Float64 = 1e-6
    rel_tol::Float64 = 1e-6
    dense::Bool = false
end

Base.@kwdef struct ModelSpec{F,U,P}
    name::String
    rhs::F
    u0::U
    p::P
    tspan::Tuple{Float64,Float64}
end

Base.@kwdef struct IntegratorSpec{A,K}
    name::String
    alg::A
    solve_kwargs::K = (;)
end

Base.@kwdef struct BenchmarkPlan{M,I}
    models::Vector{M}
    integrators::Vector{I}
    driver::BenchmarkDriver = BenchmarkDriver()
end

odeprob(model::ModelSpec) = ODEProblem(model.rhs, model.u0, model.tspan, model.p)

function run_trial(model::ModelSpec, integrator::IntegratorSpec, driver::BenchmarkDriver)
    prob = odeprob(model)
    driver_kwargs = (
        saveat=driver.saveat,
        save_everystep=driver.save_everystep,
        dense=driver.dense,
        maxiters=driver.maxiters,
        abstol=driver.abs_tol,
        reltol=driver.rel_tol,
    )

    kwargs = (; driver_kwargs..., integrator.solve_kwargs...)

    # warm-up compile (avoid counting JIT)
    solve(prob, integrator.alg; kwargs...)

    return @benchmark solve(
        $prob, $integrator.alg;
        $(kwargs...)
    )
end

_trial_metrics(tr) = begin
    est = median(tr)
    (time_ns = est.time, allocs = est.allocs, memory = est.memory)
end

function _print_raw_metrics(entries)
    println("\nRaw (median) metrics:")
    for entry in entries
        m = entry.metrics
        println(rpad(entry.label, 24),
                "  time = ", m.time_ns, " ns",
                "   allocs = ", m.allocs,
                "   bytes = ", m.memory)
    end
    return nothing
end

function _plot_normalized_bar(labels, norm_values;
                             title::AbstractString,
                             ylabel::AbstractString,
                             filename::AbstractString,
                             xrotation::Real=25,
                             yscale=:log10)
    # Avoid log(0) if something hits 0 (e.g., allocs).
    vals = max.(norm_values, 1e-12)

    p = bar(labels, vals;
        title=title,
        ylabel=ylabel,
        xrotation=xrotation,
        linewidth=0,
        yscale=yscale,
        ylims=(1.0, maximum(vals) * 1.1),
    )
    savefig(p, filename)
    return p
end

function run_benchmarks(plan::BenchmarkPlan)
    trials = Dict{String,Any}()
    metrics = Dict{String,NamedTuple}()
    labels = String[]

    for model in plan.models
        for integrator in plan.integrators
            label = "$(model.name) / $(integrator.name)"
            push!(labels, label)

            tr = Base.invokelatest(run_trial, model, integrator, plan.driver)
            trials[label] = tr
            metrics[label] = _trial_metrics(tr)
        end
    end

    labels = getproperty.(entries, :label)

    times  = [entry.metrics.time_ns for entry in entries]
    allocs = [entry.metrics.allocs  for entry in entries]

    # Normalize to best (minimum) per metric.
    time_norm  = times ./ minimum(times)
    alloc_norm = allocs ./ minimum(allocs)

    _print_raw_metrics(entries)

    p_time = _plot_normalized_bar(
        labels, time_norm;
        title="Normalized wall-clock time (median, log10)",
        ylabel="Time (× best)",
        filename="rk4_time_normalized_log10.png",
    )

    p_alloc = _plot_normalized_bar(
        labels, alloc_norm;
        title="Normalized allocations (median, log10)",
        ylabel="Allocations (× best)",
        filename="rk4_allocs_normalized_log10.png",
    )

    println("\nSaved figures: rk4_time_normalized_log10.png, rk4_allocs_normalized_log10.png")

    return (trials=trials, metrics=metrics, labels=labels, time_norm=time_norm,
            alloc_norm=alloc_norm, p_time=p_time, p_alloc=p_alloc)
end

function default_benchmark_plan(; driver=BenchmarkDriver(),
                                p_vec=[0.1, 0.1, 14], p_svec=@SVector[0.1, 0.1, 14])
    models = [
        ModelSpec(name="naive (Vector OOP)", rhs=rossler, u0=[1.0, 1.0, 1.0],
                  p=p_vec, tspan=driver.tspan),
        ModelSpec(name="in-place (!)", rhs=rossler!, u0=[1.0, 1.0, 1.0],
                  p=p_vec, tspan=driver.tspan),
        ModelSpec(name="static (SVector OOP)", rhs=rossler_static,
                  u0=@SVector [1.0, 1.0, 1.0], p=p_svec, tspan=driver.tspan),
    ]

    integrators = [
        IntegratorSpec(name="RK4 fixed-step", alg=RK4(),
                       solve_kwargs=(; adaptive=false, dt=driver.dt)),
    ]

    return BenchmarkPlan(models=models, integrators=integrators, driver=driver)
end

run_rk4_benchmarks(; driver=BenchmarkDriver(), p_vec=[0.1, 0.1, 14],
                   p_svec=@SVector[0.1, 0.1, 14]) =
    run_benchmarks(default_benchmark_plan(; driver=driver, p_vec=p_vec, p_svec=p_svec))

# If you want this to run when executed as a script:
if abspath(PROGRAM_FILE) == @__FILE__
    run_rk4_benchmarks()
end

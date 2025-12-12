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
# Fixed-step RK4 benchmark harness
# ----------------------------

Base.@kwdef struct RK4BenchDriver{S}
    tspan::Tuple{Float64,Float64} = (0.0, 50.0)
    dt::Float64 = 1e-4
    saveat::S = 0.0:0.1:50.0
    maxiters::Int = 10^9
    save_everystep::Bool = false
    abs_tol::Float64 = 1e-6
    rel_tol::Float64 = 1e-6
end

odeprob(sys; u0, p, tspan) = ODEProblem(sys, u0, tspan, p)

function benchmark_fixed_rk4(driver::RK4BenchDriver, sys; u0, p)
    prob = odeprob(sys; u0=u0, p=p, tspan=driver.tspan)
    alg = RK4()

    # warm-up compile (avoid counting JIT)
    solve(prob, alg; adaptive=false, dt=driver.dt, saveat=driver.saveat,
          save_everystep=false, dense=false, maxiters=driver.maxiters)

    return @benchmark solve(
        $prob, $alg;
        adaptive=false,
        dt=$(driver.dt),
        saveat=$(driver.saveat),
        save_everystep=$(driver.save_everystep),
        dense=false,
        maxiters=$(driver.maxiters),
        abstol=$(driver.abs_tol),
        reltol=$(driver.rel_tol),
    )
end

_trial_metrics(tr) = begin
    est = median(tr)
    (time_ns = est.time, allocs = est.allocs, memory = est.memory)
end

function _print_raw_metrics(labels::AbstractVector{<:AbstractString},
                           metrics::AbstractDict{<:AbstractString,<:NamedTuple})
    println("\nRaw (median) metrics:")
    for k in labels
        m = metrics[k]
        println(rpad(k, 24),
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
function run_rk4_benchmarks(; driver=RK4BenchDriver(),
                            p_vec=[0.1, 0.1, 14], p_svec=@SVector[0.1, 0.1, 14])

    methods = [
        ("naive (Vector OOP)",   rossler,          [1.0, 1.0, 1.0], p_vec),
        ("in-place (!)",         rossler!,         [1.0, 1.0, 1.0], p_vec),
        ("static (SVector OOP)", rossler_static,   (@SVector [1.0, 1.0, 1.0]), p_svec),
    ]

    trials = Dict{String,Any}()
    metrics = Dict{String,NamedTuple}()

    for (name, f, u0, p) in methods
        tr = benchmark_fixed_rk4(driver, f; u0=u0, p=p)
        trials[name] = tr
        metrics[name] = _trial_metrics(tr)
    end

    labels = collect(keys(metrics))

    times  = [metrics[k].time_ns for k in labels]
    allocs = [metrics[k].allocs  for k in labels]

    # Normalize to best (minimum) per metric.
    time_norm  = times ./ minimum(times)
    alloc_norm = allocs ./ minimum(allocs)

    _print_raw_metrics(labels, metrics)

    p_time = _plot_normalized_bar(
        labels, time_norm;
        title="RK4 fixed-step: Normalized wall-clock time (median, log10)",
        ylabel="Time (× best)",
        filename="rk4_time_normalized_log10.png",
    )

    p_alloc = _plot_normalized_bar(
        labels, alloc_norm;
        title="RK4 fixed-step: Normalized allocations (median, log10)",
        ylabel="Allocations (× best)",
        filename="rk4_allocs_normalized_log10.png",
    )

    println("\nSaved figures: rk4_time_normalized_log10.png, rk4_allocs_normalized_log10.png")

    return (trials=trials, metrics=metrics, time_norm=time_norm, alloc_norm=alloc_norm,
            p_time=p_time, p_alloc=p_alloc)
end

# If you want this to run when executed as a script:
if abspath(PROGRAM_FILE) == @__FILE__
    run_rk4_benchmarks()
end

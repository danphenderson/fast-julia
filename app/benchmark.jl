using BenchmarkTools
using DifferentialEquations
using StaticArrays

# Use a parametric, immutable driver to avoid abstract fields while allowing ranges for `saveat`.
struct BenchmarkDriver{S<:AbstractVector{Float64}}
    abstol::Float64
    reltol::Float64
    saveat::S
    dtmax::Float64
    dtmin::Float64
    maxiters::Int
    tspan::Tuple{Float64,Float64}
end

function benchmark_driver(;
    abstol::Float64 = 1e-6,
    reltol::Float64 = 1e-6,
    saveat::AbstractVector{Float64} = 0.0:0.1:10.0,
    dtmax::Float64 = 0.1,
    dtmin::Float64 = 1e-6,
    maxiters::Int = 10_000,
    tspan::Tuple{Float64,Float64} = (0.0, 10.0),
)
    return BenchmarkDriver(abstol, reltol, saveat, dtmax, dtmin, maxiters, tspan)
end

ode(sys; u0 = [1.0, 1.0, 1.0], tspan = (0.0, 10.0), p = [0.2, 0.2, 5.7]) = ODEProblem(sys, u0, tspan, p)

function solve_ode(
    driver::BenchmarkDriver;
    sys = rossler!,
    u0 = [1.0, 1.0, 1.0],
    tspan = driver.tspan,
    p = [0.2, 0.2, 5.7],
    alg = Tsit5(),
)
    prob = ode(sys; u0 = u0, tspan = tspan, p = p)
    return solve(
        prob,
        alg;
        abstol = driver.abstol,
        reltol = driver.reltol,
        saveat = driver.saveat,
        dtmax = driver.dtmax,
        dtmin = driver.dtmin,
        maxiters = driver.maxiters,
    )
end

"""
Benchmark solving an ODE with a given RHS `sys` and initial condition `u0`.

Note: this benchmarks only `solve(prob, ...)` and does **not** include ODEProblem construction.
"""
function benchmark_ode(
    driver::BenchmarkDriver,
    sys,
    u0 = [1.0, 1.0, 1.0];
    tspan = driver.tspan,
    p = [0.2, 0.2, 5.7],
    alg = Tsit5(),
)
    prob = ode(sys; u0 = u0, tspan = tspan, p = p)
    return @benchmark solve(
        $prob,
        $alg;
        abstol = $(driver.abstol),
        reltol = $(driver.reltol),
        saveat = $(driver.saveat),
        dtmax = $(driver.dtmax),
        dtmin = $(driver.dtmin),
        maxiters = $(driver.maxiters),
    )
end

"""
`Naive` out-of-place ODE system, allocating a new Vector on every call.

Performance notes:
- Each RK stage triggers an allocation, increasing GC pressure.
- Type-stable for concrete inputs (returns `Vector{Float64}` with Float64 inputs).
"""
function rossler(vx, vp, t)
    dx1 = -vx[2] - vx[3]
    dx2 = vx[1] + vp[1] * vx[2]
    dx3 = vp[2] + vx[1] * vx[2] - vp[3] * vx[3]
    return [dx1, dx2, dx3]
end

benchmark_naive(driver::BenchmarkDriver = benchmark_driver()) =
    benchmark_ode(driver, rossler, [1.0, 1.0, 1.0])

"""
`Naive` out-of-place ODE system with type annotations and an explicitly typed return.

Notes:
- Avoids forcing `t` to have the same type as `u`'s element type (helps AD / dual numbers).
"""
function rossler_annotated(vx::AbstractVector{T}, p::AbstractVector{T}, t) where {T<:Real}
    dx1 = -vx[2] - vx[3]
    dx2 = vx[1] + vp[1] * vx[2]
    dx3 = vp[2] + vx[1] * vx[2] - vp[3] * vx[3]
    return [dx1, dx2, dx3]
end

benchmark_annotated(driver::BenchmarkDriver = benchmark_driver()) =
    benchmark_ode(driver, rossler_annotated, [1.0, 1.0, 1.0])

"""
Rossler returning a tuple `(dx, dy, dz)`.

Warning: tuples are not generally the preferred state container type for SciML problems;
use `SVector` for an immutable small state instead.
"""
function rossler_tuple(vx, vp, t)
    dx1 = -vx[2] - vx[3]
    dx2 = vx[1] + vp[1] * vx[2]
    dx3 = vp[2] + vx[1] * vx[2] - vp[3] * vx[3]
    return (dx1, dx2, dx3)
end

benchmark_tuple(driver::BenchmarkDriver = benchmark_driver()) =
    benchmark_ode(driver, rossler_tuple, (1.0, 1.0, 1.0))

"""
In-place ODE system: writes derivatives into `du`.

Performance improvements:
- Avoids per-call allocations.
- Typically fastest for `Vector`-based states.
"""
function rossler!(dx, vx, vp, t)
    x1 = vx[1]
    dx[1] = -vx[2] - vx[3]
    dx[2] =  x1 + vp[1] * vx[2]
    dx[3] = vp[2] + x1 * vx[2] - vp[3] * vx[3]
    return nothing
end

benchmark_inplace(driver::BenchmarkDriver = benchmark_driver()) =
    benchmark_ode(driver, rossler!, [1.0, 1.0, 1.0])

"""
Static out-of-place system returning an `SVector`.

Performance improvements:
- Fixed-size, stack-allocated small container (often faster than heap vectors for tiny systems).
"""
function rossler_static(vx, vp, t)
    x1 = vx[1]
    dx1 = -vx[2] - vx[3]
    dx2 = x1 + vp[1] * vx[2]
    dx3 = vp[2] + x1 * vx[2] - vp[3] * vx[3]
    return @SVector [dx1, dx2, dx3]
end

benchmark_static(driver::BenchmarkDriver = benchmark_driver()) =
    benchmark_ode(driver, rossler_static, @SVector [1.0, 1.0, 1.0])
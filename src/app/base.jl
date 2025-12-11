# final/app/base.jl

include("types.jl")

using DifferentialEquations
using Plots
using BenchmarkTools

# -------------------------------------------------------------------
# Core Rossler ODE system compatible with AutODESystem/ODESystem
# -------------------------------------------------------------------
mutable struct RosslerSystem{T} <: ODESystem{T}
    a::T
    b::T
    c::T
    dt::T
    _state::AbstractVector{T}  # follows the convention in types.jl
end


# Convenience constructor: 3D Rossler with default initial state
function RosslerSystem(a::T, b::T, c::T, dt::T;
                       x0::T = one(T),
                       y0::T = one(T),
                       z0::T = one(T)) where {T}
    state0 = [x0, y0, z0]
    return RosslerSystem{T}(a, b, c, dt, state0)
end

# Logical components of the state; used by dims(sys)
function system(::RosslerSystem)
    return (:x, :y, :z)
end

order(::RosslerSystem) = 1

# Optional helpers for accessing components of _state
x(sys::RosslerSystem) = sys._state[1]
y(sys::RosslerSystem) = sys._state[2]
z(sys::RosslerSystem) = sys._state[3]

# -------------------------------------------------------------------
# DifferentialEquations.jl integration
# -------------------------------------------------------------------

"""
    rossler_rhs!(du, u, p, t)

In-place Rossler vector field compatible with DifferentialEquations.jl.

State:
    u = [x, y, z]

Parameters:
    p = (a, b, c)

Dynamics (matches the explicit-Euler `step!` below):
    dx/dt =  y - z
    dy/dt =  x + a*y - z
    dz/dt =  b + z*(x - c)
"""
function rossler_rhs!(du, u, p, t)
    a, b, c = p
    x, y, z = u

    du[1] = y - z
    du[2] = x + a*y - z
    du[3] = b + z*(x - c)

    return nothing
end

"""
    ode_problem(sys::RosslerSystem; tspan=(0, tf))

Build a `DifferentialEquations.ODEProblem` for the given `RosslerSystem`.
Initial condition comes from `sys._state`, parameters from `sys.a,b,c`.
"""
function ode_problem(sys::RosslerSystem{T};
                     tspan::Tuple{T,T} = (zero(T), sys.dt * T(10_000))) where {T}
    u0 = copy(sys._state)
    p  = (sys.a, sys.b, sys.c)
    return DifferentialEquations.ODEProblem(rossler_rhs!, u0, tspan, p)
end

"""
    integrate!(sys::RosslerSystem; tspan, saveat, solver, kwargs...) -> sol

Solve the Rossler IVP using DifferentialEquations.jl.

- Updates `sys._state` in-place to the final solution value.
- Returns the full solution object `sol` for post-processing / plotting.
"""
function integrate!(sys::RosslerSystem{T};
                    tspan::Tuple{T,T} = (zero(T), sys.dt * T(10_000)),
                    saveat::T = sys.dt,
                    solver = DifferentialEquations.Tsit5(),
                    kwargs...) where {T}
    prob = ode_problem(sys; tspan = tspan)
    sol  = DifferentialEquations.solve(prob, solver; saveat = saveat, kwargs...)
    # bring the internal state to the end of the trajectory
    sys._state .= sol.u[end]
    return sol
end

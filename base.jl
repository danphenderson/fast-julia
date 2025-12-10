# final/base.jl

# Define an abstract type for ODE system
# with data type parameter T for state representation
abstract type AutonomousODESystem{T} end
# Should be overloaded by concrete ODE system types
system(::AutonomousODESystem) = error("system not implemented for $(typeof(sys))")
order(::AutonomousODESystem) = error("order not implemented for $(typeof(sys))")
dt(sys::AutonomousODESystem) = error("dt not implemented for $(typeof(sys))")
dt!(sys::AutonomousODESystem, new_dt) = error("dt! not implemented for $(typeof(sys))")

# Common duck typed functions for an ODE system
eltype(::AutonomousODESystem{T}) where {T} = T
dims(sys::AutonomousODESystem) = length(system(sys)) / order(sys)
(sys::Rossler{T})(param::Symbol) where {T} = getfield(sys, param)

# Define an abstract type for ODE system backends
abstract type ODESolver end
solve(sys::AutonomousODESystem, solver::ODESolver; kwargs...) = error("solve not implemented for $(typeof(solver))")
solve!(sys::AutonomousODESystem, solver::ODESolver; kwargs...) = error("solve! not implemented for $(typeof(solver))")

# Define an abstract type for ODE solution results
abstract type ODEResult end

# Define common duck typed functions for an ODE system with a given backend
function benchmark(sys::ODESystem, backend::ODESolver; kwargs...) :: ODEResult
    """BenchMarkTools.jl trial solve of `sys` with `backend`.
    Returns the result of `solve!`.

    Args:
        sys: An instance of `ODESystem`.
        backend: An instance of `ODEBackend`.
        kwargs: Additional keyword arguments passed to `solve!`.
    """
    return solve!(sys, backend; kwargs...)
end

function animation(sys::ODESystem, backend::ODESolver; kwargs...)
    """Generate an animation for the ODE system `sys`.

    Args:
        sys: An instance of `ODESystem`.
        kwargs: Additional keyword arguments for animation.
    """
    return Types.animation(sys; kwargs...)
end
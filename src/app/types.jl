# final/app/types.jl

# Define an abstract types benchmark application
# with various backends all implementing common
# duck typed functions.
abstract type AutODESystem{T} end
abstract type ODESystem{T} <: AutODESystem{T} end
abstract type Mesh{T} end
abstract type Field{T} end
abstract type ODESolver end
abstract type ODEBackend end
abstract type ODEFieldSolver end

# Concrete types Overload:
# Convention: concrete AutODESystem subtypes are expected to have
#   - a field `_state` storing the state vector/array
#   - a field `dt` storing the time step

system(sys::AutODESystem) =
    error("system not implemented for $(typeof(sys))")

state(sys::AutODESystem) =
    getfield(sys, :_state)

state!(sys::AutODESystem, new_state) =
    (sys._state = new_state)

order(sys::AutODESystem) =
    error("order not implemented for $(typeof(sys))")

order!(sys::AutODESystem, ::Int) =
    error("order! not implemented for $(typeof(sys))")

dt(sys::AutODESystem) =
    getfield(sys, :dt)

dt!(sys::AutODESystem, new_dt) =
    (sys.dt = new_dt)

dims(sys::AutODESystem) =
    Int(length(system(sys)) ÷ order(sys))

# Could also be (ntime(sys), dims(sys)) if you define ntime
sol_shape(sys::AutODESystem) =
    (:, dims(sys))

solve(sys::AutODESystem, solver::ODESolver, mesh::Mesh; kwargs...) =
    error("solve not implemented for system=$(typeof(sys)), solver=$(typeof(solver)), mesh=$(typeof(mesh))")

solve!(sys::AutODESystem, solver::ODESolver, mesh::Mesh; kwargs...) =
    error("solve! not implemented for system=$(typeof(sys)), solver=$(typeof(solver)), mesh=$(typeof(mesh))")
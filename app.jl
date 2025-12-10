# final/app.jl

include("base.jl")

mutable struct Rossler{T} <: ODESystem{T}
    a::T
    b::T
    c::T
    dt::T
    function Rossler(a::T, b::T, c::T, dt::T) where {T}
        return Rossler{T}(a, b, c, dt)
    end
end
order(::Rossler{T}) where {T} = 1
dims(sys::Rossler{T}) where {T} = 3
sol_shape(sys::Rossler{T}) where {T} = (:, 3)
state(sys::Rossler{T}) where {T} = getfield(sys, :_state)
set_state!(sys::Rossler{T}, new_state) where {T} = (sys._state = new_state)
dt(sys::Rossler{T}) where {T} = getfield(sys, :dt)
set_dt!(sys::Rossler{T}, new_dt) where {T} = (sys.dt = new_dt)



dxdt(sys::Rossler{T}) = sys(:a) .* (state(sys)[2] .- state(sys)[1])
dydt(sys::Rossler{T}) = state(sys)[1] .+ sys(:b) .* state(sys)[2]
dzdt(sys::Rossler{T}) = sys(:c) .+ state(sys)[3] .* (state(sys)[1] .- state(sys)[4])



@base.kwdef mutable struct RosslerStepNiave
    dt::Float64 = 0.005
    a::Float64 = 0.1
    b::Float64 = 0.1
    c::Float64 = 12
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end


function step!(r::Rossler)
    dx = r.y - r.z
    dy = r.x + (r.a * r.y) - r.z
    dz = r.b + r.z * (r.x - r.c)
    r.x += r.dt * dx
    r.y += r.dt * dy
    r.z += r.dt * dz
end


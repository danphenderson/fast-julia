using Plots

@kwdef mutable struct Rossler
    dt::Float64 = 0.01
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

function main()
    attractor = Rossler()

    plt = plot3d(
        1,
        xlim = (-30, 30),
        ylim = (-30, 30),
        zlim = (0, 60),
        title = "Rossler Attractor",
        marker = 2,
    )

    # Build the animation
    anim = @animate for i = 1:1500
        step!(attractor)
        push!(plt, attractor.x, attractor.y, attractor.z)
    end every 10

    # Save GIF in current directory
    gif(anim, "rossler.gif", fps = 30)
end

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
        xlim = (-100, 100),
        ylim = (-100, 100),
        zlim = (-10, 80),
        title = "Rossler Attractor",
        xlabel = "X-axis",
        ylabel = "Y-axis",
        zlabel = "Z-axis",
        line = (:blue, 0.5),
        marker = 2,
    )

    # Build the animation
    anim = @animate for i = 1:10000
        step!(attractor)
        push!(plt, attractor.x, attractor.y, attractor.z)
    end every 10

    # Save GIF in current directory
    gif(anim, "rossler.gif", fps = 30)
end



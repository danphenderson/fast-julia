using Plots
using DifferentialEquations

@kwdef struct Mesh 
    nx::Int = 50
    ny::Int = 50
    nz::Int = 50
    x_min::Float64 = -10.0
    x_max::Float64 = 10.0
    y_min::Float64 = -10.0
    y_max::Float64 = 10.0
    z_min::Float64 = -10.0
    z_max::Float64 = 10.0
end

"""
Standard Rössler parameters + initial condition and time step.
"""
@kwdef struct Rossler
    dt::Float64 = 0.01
    a::Float64 = 0.1
    b::Float64 = 0.1
    c::Float64 = 12.0
    x0::Float64 = 1.0
    y0::Float64 = 1.0
    z0::Float64 = 1.0
end

"""
Rössler vector field in DifferentialEquations.jl in-place form.

  dx/dt = -y - z
  dy/dt =  x + a*y
  dz/dt =  b + z*(x - c)
"""
function rossler!(du, u, p, t)
    a, b, c = p
    x, y, z = u
    du[1] = -y - z
    du[2] =  x + a*y
    du[3] =  b + z*(x - c)
    return nothing
end

function main(; nsteps::Int = 10_000)
    params = Rossler()              # holds dt, parameters, and ICs

    # Initial condition and parameters for the ODEProblem
    u0    = [params.x0, params.y0, params.z0]
    p     = (params.a, params.b, params.c)
    tspan = (0.0, params.dt * nsteps)

    # Build and solve the ODE with DifferentialEquations.jl
    prob = ODEProblem(rossler!, u0, tspan, p)
    sol  = solve(prob; saveat = params.dt)  # default non-stiff solver (Tsit5())

    # Extract the trajectory
    xs = [u[1] for u in sol.u]
    ys = [u[2] for u in sol.u]
    zs = [u[3] for u in sol.u]

    # Initialize the 3D plot with the first point
    plt = plot3d(
        xs[1], ys[1], zs[1],
        xlim = (-100, 100),
        ylim = (-100, 100),
        zlim = (-10, 80),
        title = "Rössler Attractor",
        xlabel = "X-axis",
        ylabel = "Y-axis",
        zlabel = "Z-axis",
        line = (:blue, 0.5),
        marker = 2,
        legend = false,
    )

    # Build the animation by streaming points from the ODE solution
    anim = @animate for i in 2:length(xs)
        push!(plt, xs[i], ys[i], zs[i])
    end every 10

    # Save GIF in current directory
    gif(anim, "rossler.gif", fps = 30)
end

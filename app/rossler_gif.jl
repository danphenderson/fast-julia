using Plots

# ----------------------------
# Rössler system (x'=-y-z, y'=x+a y, z'=b+z(x-c))
# ----------------------------
@kwdef mutable struct Rossler
    p_vec::Vector{Float64} = [0.1, 0.1, 14.0]   # (a,b,c)
    x_vec::Vector{Float64} = [1.0, 1.0, 0.0]    # (x,y,z)
    dt::Float64 = 0.01
end

function rossler(vx, vp, t)
    dx1 = -vx[2] - vx[3]
    dx2 =  vx[1] + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (vx[1] - vp[3])
    return [dx1, dx2, dx3]
end

function rk4_step!(sys::Rossler)
    dt = sys.dt
    x  = sys.x_vec
    p  = sys.p_vec

    k1 = rossler(x, p, 0.0)
    k2 = rossler(x .+ 0.5 * dt .* k1, p, 0.5 * dt)
    k3 = rossler(x .+ 0.5 * dt .* k2, p, 0.5 * dt)
    k4 = rossler(x .+ 1.0 * dt .* k3, p, 1.0 * dt)

    sys.x_vec .+= (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    return sys
end

# Build an animation (does not write a file yet)
function rossler_animation(; nsteps=20000, burnin=2500,
                           p_vec=[0.1, 0.1, 14.0], x_vec=[1.0, 1.0, 0.0],
                           dt=0.01, camera=(30, 30))
    sys = Rossler(p_vec=Float64.(p_vec), x_vec=Float64.(x_vec), dt=dt)

    # let the trajectory settle onto the attractor
    for _ in 1:burnin
        rk4_step!(sys)
    end

    xs = Float64[sys.x_vec[1]]
    ys = Float64[sys.x_vec[2]]
    zs = Float64[sys.x_vec[3]]

    @animate for _ in 1:nsteps
        rk4_step!(sys)
        push!(xs, sys.x_vec[1]); push!(ys, sys.x_vec[2]); push!(zs, sys.x_vec[3])

        plot3d(xs, ys, zs;
            xlabel="x", ylabel="y", zlabel="z",
            title="Rössler Attractor",
            legend=false,
            camera=camera,
            linewidth=1,
            xlims=(-30, 30),
            ylims=(-30, 30),
            zlims=(-10, 60),
            aspect_ratio=:equal
        )
    end
end

# Write the gif to disk
function make_gif(; filename="rossler.gif", fps=30, kwargs...)
    anim = rossler_animation(; kwargs...)
    return Plots.gif(anim, filename; fps=fps)
end

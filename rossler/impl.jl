# rossler/impl.jl
"""
Module contains various implementations
of the Rossler system in various programming languages.

----------------------------
Rossler RHS implementations (canonical)
x' = -y - z
y' =  x + a y
z' =  b + z(x - c)
----------------------------
"""

using StaticArrays

# Naive out-of-place: allocates a new Vector each call.
function rossler(vx, vp, t)
    dx1 = -vx[2] - vx[3]
    dx2 =  vx[1] + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (vx[1] - vp[3])
    return [dx1, dx2, dx3]
end

# Out-of-place with type annotations (AD-friendly: t unconstrained).
function rossler_annotated(vx::AbstractVector{T}, vp::AbstractVector{T}, t) where {T<:Real}
    dx1 = -vx[2] - vx[3]
    dx2 =  vx[1] + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (vx[1] - vp[3])
    return T[dx1, dx2, dx3]
end

# In-place: writes into dx.
function rossler!(dx, vx, vp, t)
    x1 = vx[1]
    dx[1] = -vx[2] - vx[3]
    dx[2] =  x1 + vp[1] * vx[2]
    dx[3] =  vp[2] + vx[3] * (x1 - vp[3])
    return nothing
end

# Static out-of-place: returns SVector (often fastest for tiny systems).
function rossler_static(vx, vp, t)
    x1 = vx[1]
    dx1 = -vx[2] - vx[3]
    dx2 =  x1 + vp[1] * vx[2]
    dx3 =  vp[2] + vx[3] * (x1 - vp[3])
    return @SVector [dx1, dx2, dx3]
end



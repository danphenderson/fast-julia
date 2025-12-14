# rossler/impl.jl

# Rössler system implementations with explicit in-bounds checks.
# This file defines multiple variants of the Rössler right‑hand side.
# - `rossler`: naive out-of-place implementation returning a new `Vector`.
# - `rossler!`: in-place implementation writing into a preallocated `du`.
# - `rossler_static`: stack-allocated variant returning an `SVector`.
# - `rossler_type_stable`: type-stable out-of-place implementation that promotes
#    state and parameter element types and returns a new `Vector` with promoted element type.
# - `rossler_ad`: AD-friendly and allocation-free variant that operates on
#    `SVector` state and `SVector` parameters and returns an `SVector`.
#
# All methods are annotated with `@inline` and `@inbounds` to eliminate
# bounds checks and encourage the Julia compiler to inline tiny functions.
#
# These definitions live in a separate file so they can be easily included
# by both the benchmarking driver and any analysis scripts.

using StaticArrays

# -------------------------
# Naive out-of-place solver
# -------------------------
function rossler_naive(u, p, t)
    du1 = -u[2] - u[3]
    du2 =  u[1] + p[1] * u[2]
    du3 =  p[2] + u[3] * (u[1] - p[3])
    return [du1, du2, du3]
end

"""
    rossler(u, p, t)

Out-of-place Rössler right–hand side.  Allocates a fresh
three-element `Vector` at every call.  Accepts any `AbstractVector`
of length three for `u` and `p`.  Returns a `Vector` of the same
element type as `u`.  Uses `@inbounds` to remove bounds checks and
`@inline` to encourage inlining.
"""
@inline function rossler(u::AbstractVector, p::AbstractVector, t::Real)
    @inbounds begin
        x1 = u[1]; x2 = u[2]; x3 = u[3]
        a  = p[1]; b  = p[2]; c  = p[3]
    end
    return [ -x2 - x3,
              x1 + a*x2,
              b  + x3*(x1 - c) ]
end

# -------------------------
# In-place solver
# -------------------------
function rossler_naive!(du, u, p, t)
    x1 = u[1]
    du[1] = -u[2] - u[3]
    du[2] =  x1 + p[1] * u[2]
    du[3] =  p[2] + u[3] * (x1 - p[3])
    return nothing
end
"""
    rossler!(du, u, p, t)

In-place Rössler right–hand side.  Writes the derivative into the
preallocated vector `du`.  Accepts any `AbstractVector` of length
three for `du`, `u` and `p`.  Uses `@inbounds` for performance and
returns `nothing`.
"""
@inline function rossler!(du::AbstractVector, u::AbstractVector, p::AbstractVector, t::Real)
    @inbounds begin
        x1 = u[1]; x2 = u[2]; x3 = u[3]
        a  = p[1]; b  = p[2]; c  = p[3]
        du[1] = -x2 - x3
        du[2] = x1 + a*x2
        du[3] = b  + x3*(x1 - c)
    end
    return nothing
end

# -------------------------
# Static solver (out-of-place)
# -------------------------
function rossler_static_naive(u, p, t)
    x1 = u[1]
    du1 = -u[2] - u[3]
    du2 =  x1 + p[1] * u[2]
    du3 =  p[2] + u[3] * (x1 - p[3])
    return @SVector [du1, du2, du3]
end

"""
    rossler_static(u, p, t)

Static-array Rössler right–hand side.  Accepts a three-element
`SVector` for the state `u` and any indexable container for the
parameters `p`.  Returns an `SVector` with the same element type
as `u`.  This variant is allocation-free for tiny systems since
both the input and output live on the stack.
"""
@inline function rossler_static(u::SVector{3,T}, p, t::Real) where {T}
    @inbounds begin
        x1 = u[1]; x2 = u[2]; x3 = u[3]
        a  = p[1]; b  = p[2]; c  = p[3]
    end
    return SVector{3,T}(
        -x2 - x3,
         x1 + a*x2,
         b  + x3*(x1 - c),
    )
end

# -------------------------
# Type-stable out-of-place solver
# -------------------------
"""
    rossler_type_stable(u, p, t)

Type-stable out-of-place Rössler right–hand side.  Promotes the
element types of `u` and `p` so that the returned `Vector` has a
concrete element type.  This variant is useful for benchmarking
type stability with mixed element types or forward-mode AD.
"""
@inline function rossler_type_stable(u::AbstractVector{Tu}, p::AbstractVector{Tp}, t::Real) where {Tu,Tp}
    S = promote_type(Tu, Tp)
    @inbounds begin
        x1 = S(u[1]); x2 = S(u[2]); x3 = S(u[3])
        a  = S(p[1]); b  = S(p[2]); c  = S(p[3])
    end
    return S[
        -x2 - x3,
         x1 + a*x2,
         b  + x3*(x1 - c),
    ]
end

# -------------------------
# AD-ready allocation-free solver
# -------------------------
"""
    rossler_ad(u, p, t)

AD-friendly and allocation-free Rössler right–hand side.  Accepts
three-element `SVector` state `u` and three-element `SVector`
parameter `p`.  The result is an `SVector` with the promoted
element type of `u` and `p`.  This variant is well-suited for
forward-mode automatic differentiation and runs without heap
allocation.
"""
@inline function rossler_ad(u::SVector{3,Tu}, p::SVector{3,Tp}, t::Real) where {Tu,Tp}
    S = promote_type(Tu, Tp)
    @inbounds begin
        x1 = S(u[1]); x2 = S(u[2]); x3 = S(u[3])
        a  = S(p[1]); b  = S(p[2]); c  = S(p[3])
    end
    return SVector{3,S}(
        -x2 - x3,
         x1 + a*x2,
         b  + x3*(x1 - c),
    )
end
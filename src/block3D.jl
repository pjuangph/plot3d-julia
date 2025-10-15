module Block3D

import Statistics: mean

export Block, recompute_centroid!

"""
    Block(X, Y, Z)

Structured Plot3D block storing node coordinates on a regular I×J×K grid.

Fields
------
- `X, Y, Z :: Array{T,3}`  coordinate arrays (same size)
- `IMAX, JMAX, KMAX :: Int` number of nodes along each axis
- `cx, cy, cz :: Float64`   geometric centroid (mean of all nodes)

Notes
-----
- `IMAX/JMAX/KMAX` are **node counts** (not cells).
- Centroid is computed on construction; call `recompute_centroid!(b)` if you
  mutate coordinates later and want to refresh `(cx,cy,cz)`.
"""
mutable struct Block{T<:Real}
    X::Array{T,3}
    Y::Array{T,3}
    Z::Array{T,3}
    IMAX::Int
    JMAX::Int
    KMAX::Int
    cx::Float64
    cy::Float64
    cz::Float64

    function Block(X::Array{T,3}, Y::Array{T,3}, Z::Array{T,3}) where {T<:Real}
        size(X) == size(Y) == size(Z) || throw(ArgumentError("X, Y, Z must have identical sizes"))
        I, J, K = size(X)
        cx = mean(vec(X)); cy = mean(vec(Y)); cz = mean(vec(Z))
        new{T}(X, Y, Z, I, J, K, cx, cy, cz)
    end
end

# ----------------------------------------------------------------------
# Convenience / utilities
# ----------------------------------------------------------------------
Base.size(b::Block) = size(b.X)
Base.eltype(b::Block{T}) where {T} = T

"""
    recompute_centroid!(b::Block)

Recompute `(cx, cy, cz)` from the current coordinates.
"""
function recompute_centroid!(b::Block)
    b.cx = mean(vec(b.X))
    b.cy = mean(vec(b.Y))
    b.cz = mean(vec(b.Z))
    return b
end

end # module

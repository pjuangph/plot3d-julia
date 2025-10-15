# split_block.jl — functions live directly in Plot3D

import Statistics: mean
using .Block3D: Block

# Simple enum-ish Direction
struct Direction; axis::Symbol; end
Direction_i() = Direction(:i)
Direction_j() = Direction(:j)
Direction_k() = Direction(:k)

"""
    max_aspect_ratio(X,Y,Z, ix,jx,kx)

Estimate the maximum cell aspect ratio in a block by scanning local edge
lengths around (ix,jx,kx) cell center. Indices are **1-based cell indices**
(so valid ranges are 1:(IMAX-1), etc.). Returns a Float64.
"""
function max_aspect_ratio(X::AbstractArray{<:Real,3},
                          Y::AbstractArray{<:Real,3},
                          Z::AbstractArray{<:Real,3},
                          ix::Int, jx::Int, kx::Int)
    I,J,K = size(X)
    1 ≤ ix < I || error("ix out of range"); 1 ≤ jx < J || error("jx out of range")
    1 ≤ kx < K || error("kx out of range")
    # cell corners
    c = ((ix,jx,kx), (ix+1,jx,kx), (ix,jx+1,kx), (ix+1,jx+1,kx),
         (ix,jx,kx+1), (ix+1,jx,kx+1), (ix,jx+1,kx+1), (ix+1,jx+1,kx+1))
    pts = [(X[i,j,k], Y[i,j,k], Z[i,j,k]) for (i,j,k) in c]

    # edges meeting at c[1]
    e = []
    push!(e, ((pts[2][1]-pts[1][1]), (pts[2][2]-pts[1][2]), (pts[2][3]-pts[1][3])))
    push!(e, ((pts[3][1]-pts[1][1]), (pts[3][2]-pts[1][2]), (pts[3][3]-pts[1][3])))
    push!(e, ((pts[5][1]-pts[1][1]), (pts[5][2]-pts[1][2]), (pts[5][3]-pts[1][3])))
    lens = [sqrt(x*x+y*y+z*z) for (x,y,z) in e]
    lmin = minimum(lens); lmax = maximum(lens)
    return lmax > 0 ? lmax/lmin : Inf
end

# choose axis by longest node count
@inline _auto_axis(IMAX,JMAX,KMAX) = ((IMAX≥JMAX && IMAX≥KMAX) ? :i : (JMAX≥KMAX ? :j : :k))

# split helper: cut index range 1:N into ~nparts chunks, try to keep gcd alignment
function _split_indices(N::Int, nparts::Int, gcd_keep::Int)
    nparts = max(1, min(N, nparts))
    # target length per chunk in cells (N-1 is cells along that axis)
    chunk = max(1, fld(N-1, nparts))
    # adjust to be multiple of gcd_keep if possible
    if gcd_keep > 1
        m = max(1, round(Int, chunk / gcd_keep))
        chunk = max(1, m*gcd_keep)
    end
    # emit splits
    cuts = Int[1]
    i = 1
    while i+chunk ≤ N
        i += chunk
        push!(cuts, i)
    end
    if last(cuts) != N
        push!(cuts, N)
    end
    return cuts
end

"""
    split_blocks(blocks, ncells_per_block; direction::Union{Nothing,Direction}=nothing)

Split each block along a single axis so that each child has about `ncells_per_block`
cells, while preserving a greatest-common-denominator alignment.

Returns a **new Vector{Block}**.
"""
function split_blocks(blocks::Vector{Block},
                      ncells_per_block::Int,
                      direction::Union{Nothing,Direction}=nothing)
    ncells_per_block > 0 || error("ncells_per_block must be > 0")
    out = Block[]
    for b in blocks
        IMAX,JMAX,KMAX = size(b.X)
        # total cells
        cells = (IMAX-1)*(JMAX-1)*max(1,KMAX-1)
        parts = max(1, ceil(Int, cells / ncells_per_block))

        ax = isnothing(direction) ? _auto_axis(IMAX,JMAX,KMAX) : direction.axis
        g   = gcd(IMAX-1, gcd(JMAX-1, KMAX-1))
        if ax === :i
            cuts = _split_indices(IMAX, parts, g)
            for c in 1:length(cuts)-1
                i1, i2 = cuts[c], cuts[c+1]
                push!(out, Block(copy(b.X[i1:i2, :, :]),
                                 copy(b.Y[i1:i2, :, :]),
                                 copy(b.Z[i1:i2, :, :])))
            end
        elseif ax === :j
            cuts = _split_indices(JMAX, parts, g)
            for c in 1:length(cuts)-1
                j1, j2 = cuts[c], cuts[c+1]
                push!(out, Block(copy(b.X[:, j1:j2, :]),
                                 copy(b.Y[:, j1:j2, :]),
                                 copy(b.Z[:, j1:j2, :])))
            end
        else
            cuts = _split_indices(KMAX, parts, g)
            for c in 1:length(cuts)-1
                k1, k2 = cuts[c], cuts[c+1]
                push!(out, Block(copy(b.X[:, :, k1:k2]),
                                 copy(b.Y[:, :, k1:k2]),
                                 copy(b.Z[:, :, k1:k2])))
            end
        end
    end
    return out
end

module BlockFunctions

import ..Block3D: Block

export reduce_blocks

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

# Build a step index that always includes the last endpoint.
# Example: _stepidx(10, 3) => [1, 4, 7, 10]
function _stepidx(n::Int, step::Int)
    step <= 0 && error("step must be positive, got $step")
    idx = collect(1:step:n)
    if idx[end] != n
        push!(idx, n)
    end
    return idx
end

# Downsample one block by an integer factor `s`.
# Each dimension uses stride `s`, with last cell kept (endpoint-including).
function _reduce_block(b::Block, s::Int)
    s < 1 && error("reduce factor must be ≥ 1, got $s")
    s == 1 && return Block(b.X, b.Y, b.Z; index=b.index)  # keep identity

    ii = _stepidx(b.IMAX, s)
    jj = _stepidx(b.JMAX, s)
    kk = _stepidx(b.KMAX, s)

    IM = length(ii); JM = length(jj); KM = length(kk)

    Xr = Array{Float64,3}(undef, IM, JM, KM)
    Yr = similar(Xr); Zr = similar(Xr)

    @inbounds for k in 1:KM
        K = kk[k]
        for j in 1:JM
            J = jj[j]
            for i in 1:IM
                I = ii[i]
                Xr[i,j,k] = b.X[I,J,K]
                Yr[i,j,k] = b.Y[I,J,K]
                Zr[i,j,k] = b.Z[I,J,K]
            end
        end
    end

    return Block(Xr, Yr, Zr; index=b.index)
end

# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------

"""
    reduce_blocks(blocks::Vector{Block}, s::Int) -> Vector{Block}

Downsample every block by integer stride `s` in each index direction (i, j, k).  
The last index in each dimension is always included so block boundaries are preserved.

- Use `s = 1` to return blocks unchanged.
- Works for 2D blocks as well (`KMAX == 1`).

This is commonly used with `s = gcd(IMAX-1, JMAX-1, KMAX-1)` across blocks to ensure
face grids align for connectivity search.
"""
function reduce_blocks(blocks::Vector{Block}, s::Int)
    s < 1 && error("reduce_blocks: stride must be ≥ 1, got $s")
    out = Vector{Block}(undef, length(blocks))
    @inbounds for i in eachindex(blocks)
        out[i] = _reduce_block(blocks[i], s)
    end
    return out
end

end # module

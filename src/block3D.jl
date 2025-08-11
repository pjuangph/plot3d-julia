module Block3D

import Statistics: mean

"""
    struct Block

Stores a 3D or 2D block of coordinates X, Y, Z with shape (IMAX, JMAX, KMAX).

If `KMAX == 1`, this block is treated as 2D (Option A).
"""
struct Block
    X::Array{Float64,3}
    Y::Array{Float64,3}
    Z::Array{Float64,3}
    IMAX::Int
    JMAX::Int
    KMAX::Int
    index::Int
end

"""
    Block(X, Y, Z; index=0)

Construct a `Block` given coordinate arrays.  
Automatically sets IMAX, JMAX, KMAX.
"""
function Block(X::Array{Float64,3}, Y::Array{Float64,3}, Z::Array{Float64,3}; index::Int=0)
    sizeX = size(X)
    return Block(X, Y, Z, sizeX[1], sizeX[2], sizeX[3], index)
end

"""
    centroid(block::Block) -> NTuple{3,Float64}

Return centroid (mean X, mean Y, mean Z) for the block.
"""
function centroid(b::Block)
    return (mean(b.X), mean(b.Y), mean(b.Z))
end

"""
    size_tuple(block::Block) -> NTuple{3,Int}

Return (IMAX, JMAX, KMAX) for the block.
"""
size_tuple(b::Block) = (b.IMAX, b.JMAX, b.KMAX)

end # module

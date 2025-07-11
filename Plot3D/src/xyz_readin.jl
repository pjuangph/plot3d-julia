module xyzreader

using StaticArrays, Test
using LinearAlgebra
using ..BlockType

export Block3DXYZ, read_xyz_block, read_structured_xyz 

struct Block3DXYZ
    x::Array{Float64,3}
    y::Array{Float64,3}
    z::Array{Float64,3}
end
"""
    read_xyz_block(start_line::Int, dims::NTuple{3, Int}, lines::Vector{String}) -> Tuple{Block3DXYZ, Int}

    Takes expected number of values from lines 2 through NBlocks and checks that actual number of values 
"""

function read_xyz_block(start_line::Int, dims::NTuple{3, Int}, lines::Vector{String})
    ni, nj, nk = dims
    x = zeros(Float64, ni, nj, nk)
    y = zeros(Float64, ni, nj, nk)
    z = zeros(Float64, ni, nj, nk)

    idx = start_line

    for k in 1:nk, j in 1:nj
        vals = parse.(Float64, split(lines[idx]))
        if length(vals) != ni
            error("Line $idx: Expected $ni values, got $(length(vals)). Line: $(lines[idx])")
        end
        x[:, j, k] .= vals
        idx += 1
    end
    for k in 1:nk, j in 1:nj
        vals = parse.(Float64, split(lines[idx]))
        if length(vals) != ni
            error("Line $idx: Expected $ni values, got $(length(vals)). Line: $(lines[idx])")
        end
        y[:, j, k] .= vals
        idx += 1
    end
    for k in 1:nk, j in 1:nj
        vals = parse.(Float64, split(lines[idx]))
        if length(vals) != ni
            error("Line $idx: Expected $ni values, got $(length(vals)). Line: $(lines[idx])")
        end
        z[:, j, k] .= vals
        idx += 1
    end

    return Block3DXYZ(x, y, z), idx
end

"""
    read_structured_xyz(filename::String) -> Vector{Block3DXYZ}

Reads a structured mesh from a `.xyz` file into a vector of Block3DXYZ structs.
"""
function read_structured_xyz(filename::String)
    lines = readlines(filename)
    nblocks = parse(Int, lines[1])
    dims = [Tuple(parse.(Int, split(line))) for line in lines[2:1 + nblocks]]

    blocks = Block3DXYZ[]
    line_idx = 2 + nblocks

    for d in dims
        block, line_idx = read_xyz_block(line_idx, d, lines)
        push!(blocks, block)
    end

    return blocks
end

end

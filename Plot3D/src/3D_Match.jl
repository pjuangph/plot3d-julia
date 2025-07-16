module Faces

using StaticArrays, Test
using LinearAlgebra
using ..BlockType
using ..xyzreader: Block3DXYZ, read_structured_xyz, read_xyz_block

export forward_match, reverse_match, final_match_check, FaceND,
    Block2DFaces, Block3DFaces, print_connections, compare_blocks, block_xyz_faces


"""Generalized FaceND for 2D or 3D """
struct FaceND{N,T}
    points::Vector{SVector{N,T}}
    x1::SVector{N,T}
    x2::SVector{N,T}
end

"""Constructor for FaceND from points"""
function FaceND(points::Vector{SVector{N,T}}) where {N,T}
    FaceND{N,T}(points, points[1], points[end])
end

Base.:(==)(F1::FaceND, F2::FaceND) = isapprox(F1, F2)

function Base.isapprox(F1::FaceND, F2::FaceND)
    return (
        (isapprox(F1.x1, F2.x1) && isapprox(F1.x2, F2.x2)) ||
        (isapprox(F1.x1, F2.x2) && isapprox(F1.x2, F2.x1))
    )
end

""" Check if two faces match in forward or reverse orientation.
   Returns true if they match, false otherwise.
   Prints the matching pairs if they match.
"""
function forward_match(f1::FaceND, f2::FaceND)
    n = length(f1.points)
    n == length(f2.points) || return false
    # for i in 1:n
    #     println(" 1   f1: ", f1.points[i], "  f2: ", f2.points[i])
    # end
    forward_match = all(f1.points[i] == f2.points[i] for i in 1:n)
    if forward_match
        # println(" 2 => Match (forward orientation)")
        return true
    end  
    return false
end

""" Check if two faces match in reverse orientation.
   Returns true if they match, false otherwise.
   Prints the matching pairs if they match.
"""
function reverse_match(f1::FaceND, f2::FaceND)
    n = length(f1.points)
    n == length(f2.points) || return false
    # for i in 1:n
    #     println(" 3   f1: ", f1.points[i], "  f2: ", f2.points[n-i+1])
    # end
    reverse_match = all(f1.points[i] == f2.points[n-i+1] for i in 1:n)
    if reverse_match
        # println(" 4 => Match (reverse orientation)")
        return true
    end  
    return false
end

""" Check if two faces match in either forward or reverse orientation.
   Prints the matching pairs if they match.
   Returns true if they match, false otherwise.
"""
function final_match_check(f1::FaceND, f2::FaceND)
    if forward_match(f1, f2)
        #println("Matches (forward) => ", f1.x1, " to ", f2.x1, " and ", f1.x2, " to ", f2.x2)
        return true
    elseif reverse_match(f1, f2)
        #println("Matches (reverse) => ", f1.x1, " to ", f2.x2, " and ", f1.x2, " to ", f2.x1)
        return true
    else
        # println(" => No match")
        return false
    end
end


"""Block2DFaces and Block3DFaces structs for 2D and 3D blocks
   These structs contain faces of the blocks in a dictionary format.
   The keys are symbols representing the face orientation (ilo, ihi, jlo, jhi for 2D;
   klo, khi are added for 3D).
"""
struct Block2DFaces
    faces::Dict{Symbol, FaceND{2,Float64}}
end

function Block2DFaces(b::Block2D)
    faces = Dict{Symbol, FaceND{2,Float64}}()
    faces[:ilo] = FaceND([SVector{2,Float64}(b.X[1, j], b.Y[1, j]) for j in 1:size(b.X, 2)])
    faces[:ihi] = FaceND([SVector{2,Float64}(b.X[end, j], b.Y[end, j]) for j in 1:size(b.X, 2)])
    faces[:jlo] = FaceND([SVector{2,Float64}(b.X[i, 1], b.Y[i, 1]) for i in 1:size(b.X, 1)])
    faces[:jhi] = FaceND([SVector{2,Float64}(b.X[i, end], b.Y[i, end]) for i in 1:size(b.X, 1)])
    return Block2DFaces(faces)
end

# 3D block faces
struct Block3DFaces
    faces::Dict{Symbol, FaceND{3,Float64}}
end

""" Constructs Block3DFaces from a Block3D object.
   The faces are extracted from the X, Y, Z coordinates of the block.
   The keys are symbols representing the face orientation (ilo, ihi, jlo, jhi, klo, khi).
"""
function Block3DFaces(b::Block3D)
    faces = Dict{Symbol, FaceND{3,Float64}}()
    # ilo: i=1, all j,k
    faces[:ilo] = FaceND([SVector{3,Float64}(b.X[1, j, k], b.Y[1, j, k], b.Z[1, j, k]) for k in 1:size(b.X, 3), j in 1:size(b.X, 2)])
    # ihi: i=end, all j,k
    faces[:ihi] = FaceND([SVector{3,Float64}(b.X[end, j, k], b.Y[end, j, k], b.Z[end, j, k]) for k in 1:size(b.X, 3), j in 1:size(b.X, 2)])
    # jlo: all i, j=1, k
    faces[:jlo] = FaceND([SVector{3,Float64}(b.X[i, 1, k], b.Y[i, 1, k], b.Z[i, 1, k]) for k in 1:size(b.X, 3), i in 1:size(b.X, 1)])
    # jhi: all i, j=end, k
    faces[:jhi] = FaceND([SVector{3,Float64}(b.X[i, end, k], b.Y[i, end, k], b.Z[i, end, k]) for k in 1:size(b.X, 3), i in 1:size(b.X, 1)])
    # klo: all i, j, k=1
    faces[:klo] = FaceND([SVector{3,Float64}(b.X[i, j, 1], b.Y[i, j, 1], b.Z[i, j, 1]) for j in 1:size(b.X, 2), i in 1:size(b.X, 1)])
    # khi: all i, j, k=end
    faces[:khi] = FaceND([SVector{3,Float64}(b.X[i, j, end], b.Y[i, j, end], b.Z[i, j, end]) for j in 1:size(b.X, 2), i in 1:size(b.X, 1)])
    return Block3DFaces(faces)
end

"""Compare all faces between two blocks, return the matching face symbols if found"""
function compare_blocks(block1, block2; block1_label="", block2_label="")
    matches = []
    for (f1dim, f1) in block1.faces
        for (f2dim, f2) in block2.faces
            if final_match_check(f1, f2)
                push!(matches, (f1dim, f2dim))
            end
        end
    end
    return matches
end

""" Print final connections between blocks """
function print_connections(block_faces)
    block_connections = []
    for i in 1:length(block_faces)
        for j in i+1:length(block_faces)
            results = compare_blocks(block_faces[i], block_faces[j]; block1_label="blk$i", block2_label="blk$j")
            for (f1dim, f2dim) in results
                push!(block_connections, ((Symbol("blk$i"), f1dim) => (Symbol("blk$j"), f2dim)))
            end
        end
    end
    return tuple(block_connections...)
end


"""Basically the same as Block3DFaces, but for Block3DXYZ.
   Used for reading structured XYZ files instead of P3D files.
"""
function block_xyz_faces(b::Block3DXYZ)
    faces = Dict{Symbol, FaceND{3,Float64}}()
    faces[:ilo] = FaceND([SVector{3,Float64}(b.x[1, j, k], b.y[1, j, k], b.z[1, j, k]) for k in 1:size(b.x, 3) for j in 1:size(b.x, 2)])
    faces[:ihi] = FaceND([SVector{3,Float64}(b.x[end, j, k], b.y[end, j, k], b.z[end, j, k]) for k in 1:size(b.x, 3) for j in 1:size(b.x, 2)])
    faces[:jlo] = FaceND([SVector{3,Float64}(b.x[i, 1, k], b.y[i, 1, k], b.z[i, 1, k]) for k in 1:size(b.x, 3) for i in 1:size(b.x, 1)])
    faces[:jhi] = FaceND([SVector{3,Float64}(b.x[i, end, k], b.y[i, end, k], b.z[i, end, k]) for k in 1:size(b.x, 3) for i in 1:size(b.x, 1)])
    faces[:klo] = FaceND([SVector{3,Float64}(b.x[i, j, 1], b.y[i, j, 1], b.z[i, j, 1]) for j in 1:size(b.x, 2) for i in 1:size(b.x, 1)])
    faces[:khi] = FaceND([SVector{3,Float64}(b.x[i, j, end], b.y[i, j, end], b.z[i, j, end]) for j in 1:size(b.x, 2) for i in 1:size(b.x, 1)])
    return Block3DFaces(faces)
end


end
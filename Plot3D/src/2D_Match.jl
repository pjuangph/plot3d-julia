module Faces

using StaticArrays, Test
using LinearAlgebra
using ..BlockType

export forward_match, reverse_match, final_match_check, Face2D,
    Block2DFaces, print_connections, compare_blocks

struct Face2D{T}
    points::Vector{SVector{2,T}}
    x1::SVector{2,T}
    x2::SVector{2,T}
end

function Face2D(points::Vector{SVector{2,T}}) where T
    Face2D(points, points[1], points[end])
end

Base.:(==)(F1::Face2D, F2::Face2D) = isapprox(F1, F2)

function Base.isapprox(F1::Face2D, F2::Face2D)
  return (
    (isapprox(F1.x1, F2.x1) && isapprox(F1.x2, F2.x2)) ||
    (isapprox(F1.x1, F2.x2) && isapprox(F1.x2, F2.x1))
  )
end

function forward_match(f1::Face2D, f2::Face2D)
    n = length(f1.points)
    n == length(f2.points) || return false
    for i in 1:n
        println(" 1   f1: ", f1.points[i], "  f2: ", f2.points[i])
    end
    forward_match = all(f1.points[i]==f2.points[i] for i in 1:n)
    if forward_match
        println(" 2 => Match (forward orientation)")
        return true
    end  
    return false
end

function reverse_match(f1::Face2D, f2::Face2D)
    n = length(f1.points)
    n == length(f2.points) || return false
    for i in 1:n
        println(" 3   f1: ", f1.points[i], "  f2: ", f2.points[n-i+1])
    end
    reverse_match = all(f1.points[i]==f2.points[n-i+1] for i in 1:n)
    if reverse_match
        println(" 4 => Match (reverse orientation)")
        return true
    end  
    return false
end

function final_match_check(f1::Face2D, f2::Face2D)
    if forward_match(f1, f2)
        println("Matches (forward) => ", f1.x1, " to ", f2.x1, " and ", f1.x2, " to ", f2.x2)
        return true
    elseif reverse_match(f1, f2)
      println("Matches (reverse) => ", f1.x1, " to ", f2.x2, " and ", f1.x2, " to ", f2.x1)
      return true
    else
        #println(" => No match")
        return false
    end
  end

struct Block2DFaces
    faces::Dict{Symbol, Face2D}
end

function Block2DFaces(b::Block2D)
    faces = Dict{Symbol, Face2D}()
    # ilo: left edge (i=1, all j)
    faces[:ilo] = Face2D([SVector(b.X[1, j], b.Y[1, j]) for j in 1:size(b.X, 2)])
    # ihi: right edge (i=end, all j)
    faces[:ihi] = Face2D([SVector(b.X[end, j], b.Y[end, j]) for j in 1:size(b.X, 2)])
    # jlo: bottom edge (all i, j=1)
    faces[:jlo] = Face2D([SVector(b.X[i, 1], b.Y[i, 1]) for i in 1:size(b.X, 1)])
    # jhi: top edge (all i, j=end)
    faces[:jhi] = Face2D([SVector(b.X[i, end], b.Y[i, end]) for i in 1:size(b.X, 1)])
    return Block2DFaces(faces)
end

# Compare all faces between two blocks, return the matching face symbols if found
function compare_blocks(block1::Block2DFaces, block2::Block2DFaces; block1_label="", block2_label="")
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

# Print all connections between blocks
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

end
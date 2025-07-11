module Faces

using StaticArrays, Test
using LinearAlgebra
using ..BlockType

export Block2DFaces, print_connections, compare_blocks

# Face2D now stores all points along the face
struct Face2D{T}
    points::Vector{SVector{2,T}}
end

# Compare two faces by all points (with tolerance), allowing for reversed orientation
function faces_match(f1::Face2D, f2::Face2D; tol=1e-8, label1="", label2="", block1="", block2="")
    n = length(f1.points)
    n == length(f2.points) || return false
    println("Comparing $block1:$label1 to $block2:$label2")
    # Print coordinates being compared (forward)
    println("  Forward orientation:")
    for i in 1:n
        println("    f1: ", f1.points[i], "  f2: ", f2.points[i])
    end
    forward_match = all(norm(f1.points[i] - f2.points[i]) < tol for i in 1:n)
    if forward_match
        println("  => Match (forward orientation)")
        return true
    end
    # Print coordinates being compared (reverse)
    println("  Reverse orientation:")
    for i in 1:n
        println("    f1: ", f1.points[i], "  f2: ", f2.points[n-i+1])
    end
    reverse_match = all(norm(f1.points[i] - f2.points[n-i+1]) < tol for i in 1:n)
    if reverse_match
        println("  => Match (reverse orientation)")
        return true
    end
    println("  => No match")
    return false
end

# Block2DFaces now stores faces as Dict{Symbol, Face2D}
struct Block2DFaces
    faces::Dict{Symbol, Face2D}
end

# Construct Block2DFaces for a 2D block, extracting all points along each face
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
            if faces_match(f1, f2; label1=string(f1dim), label2=string(f2dim), block1=block1_label, block2=block2_label)
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
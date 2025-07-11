#DOESN'T WORK FOR ROUNDS - CORRECT FILE IS Match_Rounds.jl
module Faces

using StaticArrays, Test

using ..BlockType

export Block2DFaces, print_connections, compare_blocks

struct Face2D{T}
  x1::SVector{2,T}
  x2::SVector{2,T}
end

Base.:(==)(F1::Face2D, F2::Face2D) = isequal(F1, F2)

function Base.isequal(F1::Face2D, F2::Face2D)
  return (
    (isequal(F1.x1, F2.x1) && isequal(F1.x2, F2.x2)) ||
    (isequal(F1.x1, F2.x2) && isequal(F1.x2, F2.x1))
  )
end

struct Block2DFaces
  faces::Dict{Symbol, Face2D}
end

function Block2DFaces(b::Block2D)
  println("Block dimensions: ", size(b.X))
  faces = Dict{Symbol, Face2D}()
  faces[:ilo] = Face2D(
    SVector(b.X[begin, begin], b.Y[begin, begin]),
    SVector(b.X[end, begin], b.Y[end, begin]))
    println("ilo: ", faces[:ilo]) 
  faces[:ihi] = Face2D(
    SVector(b.X[begin, end], b.Y[begin, end]),
    SVector(b.X[end, end], b.Y[end, end])) 
    println("ihi: ", faces[:ihi])
  faces[:jlo] = Face2D(
    SVector(b.X[begin, begin], b.Y[begin, begin]),
    SVector(b.X[begin, end], b.Y[begin, end])) 
    println("jlo: ", faces[:jlo])
  faces[:jhi] = Face2D(
    SVector(b.X[end, begin], b.Y[end, begin]),
    SVector(b.X[end, end], b.Y[end, end]))
    println("jhi: ", faces[:jhi])
  return Block2DFaces(faces)
end

function compare_blocks(block1::Block2DFaces, block2::Block2DFaces)
  for (f1dim, f1) in block1.faces
    for (f2dim, f2) in block2.faces
      if f1 == f2
        return f1dim, f2dim
      end
    end
  end
  return nothing, nothing
end

function print_connections(block_faces)
  block_connections = []
  for i in 1:length(block_faces)
    for j in i+1:length(block_faces)
        result = compare_blocks(block_faces[i], block_faces[j])  # Store the result
        if result != (nothing, nothing)
          push!(block_connections, ((Symbol("blk$i"), result[1]) => (Symbol("blk$j"), result[2])))
        end
    end
  end
  return tuple(block_connections...)
end
struct Face3D{T}
  x1::SVector{3,T}
  x2::SVector{3,T}
end
Base.:(==)(F1::Face3D, F2::Face3D) = isequal(F1, F2)
function Base.isequal(F1::Face3D, F2::Face3D)
  return (
    (isequal(F1.x1, F2.x1) && isequal(F1.x2, F2.x2)) ||
    (isequal(F1.x1, F2.x2) && isequal(F1.x2, F2.x1))
  )
end
struct Block3DFaces
  faces::Dict{Symbol, Face3D}

end

function Block3DFaces(b::Block3D)
  println("Block dimensions: ", size(b.X))
  faces = Dict{Symbol, Face3D}()
  faces[:ilo] = Face3D(
    SVector(b.X[begin, begin, begin], b.Y[begin, begin, begin], b.Z[begin, begin, begin]),
    SVector(b.X[end, begin, begin], b.Y[end, begin, begin], b.Z[end, begin, begin]))
  faces[:ihi] = Face3D(
    SVector(b.X[begin, end, end], b.Y[begin, end, end], b.Z[begin, end, end]),
    SVector(b.X[end, end, end], b.Y[end, end, end], b.Z[end, end, end]))
  faces[:jlo] = Face3D(
    SVector(b.X[begin, begin, begin], b.Y[begin, begin, begin], b.Z[begin, begin, begin]),
    SVector(b.X[begin, end, begin], b.Y[begin, end, begin], b.Z[begin, end, begin]))
  faces[:jhi] = Face3D(
    SVector(b.X[end, begin, begin], b.Y[end, begin, begin], b.Z[end, begin, begin]),
    SVector(b.X[end, end, end], b.Y[end, end, end], b.Z[end,end,end]))
  faces[:klo] = Face3D(
    SVector(b.X[begin, begin, begin], b.Y[begin, begin, begin], b.Z[begin, begin, begin]),
    SVector(b.X[begin, begin, end], b.Y[begin, begin, end], b.Z[begin, begin, end]))
  faces[:khi] = Face3D(
    SVector(b.X[end, begin, begin], b.Y[end, begin, begin], b.Z[end, begin, begin]),
    SVector(b.X[end, end, end], b.Y[end, end, end], b.Z[end,end,end]))
  return Block3DFaces(faces)
end

function compare_3D_blocks(block1::Block3DFaces, block2::Block3DFaces)
  for (f1dim, f1) in block1.faces
    for (f2dim, f2) in block2.faces
      if f1 == f2
        return f1dim, f2dim
      end
    end
  end
  return nothing, nothing
end

function print_3D_connections(block_faces)
  block_connections = []
  for i in 1:length(block_faces)
    for j in i+1:length(block_faces)
        result = compare_3D_blocks(block_faces[i], block_faces[j])  # Store the result
        if result != (nothing, nothing)
          push!(block_connections, ((Symbol("blk$i"), result[1]) => (Symbol("blk$j"), result[2])))
        end
    end
  end
  return tuple(block_connections...)
end

end
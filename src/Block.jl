module BlockType

export Block2D, Block3D, read_blocks

abstract type AbstractBlock end

"""
    Block2D{T}(X, Y)

Represents a 2D computational block with coordinate fields `X` and `Y`, each a 2D array of type `T`.

# Fields
- `X::AbstractArray{T,2}`: x-coordinate values
- `Y::AbstractArray{T,2}`: y-coordinate values
"""
struct Block2D{T, AX<:AbstractArray{T,2}, AY<:AbstractArray{T,2}} <: AbstractBlock
    X::AX
    Y::AY
end


"""
    Block3D{T, AA<:AbstractArray{T,3}}

A 3D Plot3D block holding coordinate arrays `X`, `Y`, and `Z`.
"""
struct Block3D{T,AA<:AbstractArray{T,3}} <: AbstractBlock
  X::AA
  Y::AA
  Z::AA
end

Base.size(b::AbstractBlock) = size(b.X)

"""
    read_block(starting_line::Int, block_dims::NTuple{3,Int}, data::Vector{String}; T=Float64)

Read the coordinates of a single 3D block from the data vector. Returns a `Block3D`.
"""
function read_block(
  starting_line::Int, block_dims::NTuple{3,Int}, data::Vector{String}; T=Float64
)
  x = zeros(T, block_dims)
  y = zeros(T, block_dims)
  z = zeros(T, block_dims)

  ni, nj, nk = block_dims

  line_number = starting_line

  for k in 1:nk
    for j in 1:nj
      x1d = @view x[:, j, k]
      x1d .= parse.(T, split(data[line_number]))
      line_number += 1
    end
  end

  for k in 1:nk
    for j in 1:nj
      y1d = @view y[:, j, k]
      y1d .= parse.(T, split(data[line_number]))
      line_number += 1
    end
  end

  for k in 1:nk
    for j in 1:nj
      z1d = @view z[:, j, k]
      z1d .= parse.(T, split(data[line_number]))
      line_number += 1
    end
  end

  return Block3d(x, y, z)
end

"""
    read_block(starting_line::Int, block_dims::NTuple{2,Int}, data::Vector{String}; T=Float64)

Read the coordinates of a single 2D block from the data vector. Returns a `Block2D`.
"""
function read_block(
  starting_line::Int, block_dims::NTuple{2,Int}, data::Vector{String}; T=Float64
)
  x = zeros(T, block_dims)
  y = zeros(T, block_dims)

  ni, nj = block_dims

  line_number = starting_line

  for j in 1:nj
    x1d = @view x[:, j]
    x1d .= parse.(T, split(data[line_number]))
    line_number += 1
  end

  for j in 1:nj
    y1d = @view y[:, j]
    y1d .= parse.(T, split(data[line_number]))
    line_number += 1
  end

  return Block2D(x, y)
end

"""
    read_blocks(filename::String)

Read all blocks from a Plot3D file. Returns a vector of `Block2D` or `Block3D` depending on the file contents.
"""
function read_blocks(filename)
  data = readlines(filename)

  n_blocks = parse(Int, data[1])
  block_dims_str = split.(data[2:(2 + n_blocks - 1)])

  first_block_dim = parse.(Int, block_dims_str[1])

  if last(first_block_dim) == 1
    blocks = read_2d_blocks(data)
  else
    blocks = read_3d_blocks(data)
  end

  return blocks
end

"""
    read_3d_blocks(data::Vector{String})

Parse all 3D blocks from the given Plot3D data vector. Returns a vector of `Block3D`.
"""
function read_3d_blocks(data)
  n_blocks = parse(Int, data[1])
  block_dims_str = split.(data[2:(2 + n_blocks - 1)])
  block_dims = Vector{NTuple{3,Int}}(undef, n_blocks)

  for i in 1:n_blocks
    block_dims[i] = Tuple(parse.(Int, block_dims_str[i]))
  end

  # how many lines does each block have in the file? the 3 is for x, y, and z coords
  lines_per_block = [block_dim[2] * block_dim[3] * 3 for block_dim in block_dims]
  first_line = n_blocks + 2 # starting line of the block coordinates

  # starting line in 'data' for each block
  starting_lines = [first_line, first_line .+ cumsum(lines_per_block)...]

  blocks = Vector{Block3D}(undef, n_blocks)
  for b in 1:n_blocks
    blocks[b] = read_block(starting_lines[b], block_dims[b], data)
  end

  return blocks
end

"""
    read_2d_blocks(data::Vector{String})

Parse all 2D blocks from the given Plot3D data vector. Returns a vector of `Block2D`.
"""
function read_2d_blocks(data)
  n_blocks = parse(Int, data[1])
  block_dims_str = split.(data[2:(2 + n_blocks - 1)])
  block_dims = Vector{NTuple{2,Int}}(undef, n_blocks)

  for i in 1:n_blocks
    block_dims[i] = Tuple(parse.(Int, block_dims_str[i])[1:2])
  end

  # how many lines does each block have in the file? the 3 is for x, y, and z coords
  # even though this is 2D, the format still contains the z points
  lines_per_block = [block_dim[2] * 3 for block_dim in block_dims]
  first_line = n_blocks + 2 # starting line of the block coordinates

  # starting line in 'data' for each block
  starting_lines = [first_line, first_line .+ cumsum(lines_per_block)...]

  blocks = Vector{Block2D}(undef, n_blocks)
  for b in 1:n_blocks
    blocks[b] = read_block(starting_lines[b], block_dims[b], data)
  end

  return blocks
end

end
module BlockType

# Exported types and functions for use in other modules
export Block2D, Block3D, read_blocks

# Abstract type for blocks (for type hierarchy)
abstract type AbstractBlock end

# 3D block type, parameterized by element type and array type
struct Block3D{T,AA<:AbstractArray{T,3}} <: AbstractBlock
  X::AA  # X coordinates (3D array)
  Y::AA  # Y coordinates (3D array)
  Z::AA  # Z coordinates (3D array)
end

# 2D block type, parameterized by element type and array type
struct Block2D{T,AA<:AbstractArray{T,2}} <: AbstractBlock
  X::AA  # X coordinates (2D array)
  Y::AA  # Y coordinates (2D array)
end

# Generic size function for blocks (returns size of X array)
Base.size(b::AbstractBlock) = size(b.X)

"""
    read_block(starting_line::Int, block_dims::NTuple{N,Int}, data::Vector{String})

Read the coordinates of a single block from the data vector. The `starting_line` is the line
that the coordinates for this block start. The `data` vector is the raw string read from the file.
"""

# Read a single 3D block from the data vector
function read_block(
  starting_line::Int, block_dims::NTuple{3,Int}, data::Vector{String}, T=Float64
)
  x = zeros(T, block_dims)  # Allocate X array
  y = zeros(T, block_dims)  # Allocate Y array
  z = zeros(T, block_dims)  # Allocate Z array

  ni, nj, nk = block_dims  # Unpack dimensions

  line_number = starting_line  # Track current line in data

  # Read X coordinates
  for k in 1:nk
    for j in 1:nj
      x1d = @view x[:, j, k]
      x1d .= parse.(T, split(data[line_number]))
      line_number += 1
    end
  end

  # Read Y coordinates
  for k in 1:nk
    for j in 1:nj
      y1d = @view y[:, j, k]
      y1d .= parse.(T, split(data[line_number]))
      line_number += 1
    end
  end

  # Read Z coordinates
  for k in 1:nk
    for j in 1:nj
      z1d = @view z[:, j, k]
      z1d .= parse.(T, split(data[line_number]))
      line_number += 1
    end
  end

  return Block3d(x, y, z)  # Return a Block3D object
end

# Read a single 2D block from the data vector
function read_block(
  starting_line::Int, block_dims::NTuple{2,Int}, data::Vector{String}, T=Float64
)
  x = zeros(T, block_dims)  # Allocate X array
  y = zeros(T, block_dims)  # Allocate Y array

  ni, nj = block_dims  # Unpack dimensions

  line_number = starting_line  # Track current line in data

  # Read X coordinates
  for j in 1:nj
    x1d = @view x[:, j]
    x1d .= parse.(T, split(data[line_number]))
    line_number += 1
  end

  # Read Y coordinates
  for j in 1:nj
    y1d = @view y[:, j]
    y1d .= parse.(T, split(data[line_number]))
    line_number += 1
  end

  return Block2D(x, y)  # Return a Block2D object
end

"""
    read_blocks(filename)

Read all of the blocks from the given Plot3D file
"""
function read_blocks(filename)
  data = readlines(filename)  # Read all lines from the file

  n_blocks = parse(Int, data[1])  # Number of blocks
  block_dims_str = split.(data[2:(2 + n_blocks - 1)])  # Block dimensions as strings

  first_block_dim = parse.(Int, block_dims_str[1])  # Dimensions of the first block

  # Decide if the file contains 2D or 3D blocks based on the last dimension
  if last(first_block_dim) == 1
    blocks = read_2d_blocks(data)
  else
    blocks = read_3d_blocks(data)
  end

  return blocks  # Return all blocks
end

# Read all 3D blocks from the data vector
function read_3d_blocks(data)
  n_blocks = parse(Int, data[1])  # Number of blocks
  block_dims_str = split.(data[2:(2 + n_blocks - 1)])  # Block dimensions as strings
  block_dims = Vector{NTuple{3,Int}}(undef, n_blocks)  # Allocate block dimensions

  # Parse block dimensions for each block
  for i in 1:n_blocks
    block_dims[i] = Tuple(parse.(Int, block_dims_str[i]))
  end

  # Calculate number of lines per block (for X, Y, Z)
  lines_per_block = [block_dim[2] * block_dim[3] * 3 for block_dim in block_dims]
  first_line = n_blocks + 2 # Starting line of the block coordinates

  # Calculate starting line for each block in the data vector
  starting_lines = [first_line, first_line .+ cumsum(lines_per_block)...]

  blocks = Vector{Block3D}(undef, n_blocks)  # Allocate blocks
  for b in 1:n_blocks
    blocks[b] = read_block(starting_lines[b], block_dims[b], data)
  end

  return blocks  # Return all 3D blocks
end

# Read all 2D blocks from the data vector
function read_2d_blocks(data)
  n_blocks = parse(Int, data[1])  # Number of blocks
  block_dims_str = split.(data[2:(2 + n_blocks - 1)])  # Block dimensions as strings
  block_dims = Vector{NTuple{2,Int}}(undef, n_blocks)  # Allocate block dimensions

  # Parse block dimensions for each block (only first two dimensions for 2D)
  for i in 1:n_blocks
    block_dims[i] = Tuple(parse.(Int, block_dims_str[i])[1:2])
  end

  # Calculate number of lines per block (for X, Y, Z)
  # Even though this is 2D, the format still contains the z points
  lines_per_block = [block_dim[2] * 3 for block_dim in block_dims]
  first_line = n_blocks + 2 # Starting line of the block coordinates

  # Calculate starting line for each block in the data vector
  starting_lines = [first_line, first_line .+ cumsum(lines_per_block)...]

  blocks = Vector{Block2D}(undef, n_blocks)  # Allocate blocks
  for b in 1:n_blocks
    blocks[b] = read_block(starting_lines[b], block_dims[b], data)
  end

  return blocks  # Return all 2D blocks
end

end
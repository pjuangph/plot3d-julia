"""
    Run_File3D.jl

A generic script to read a Plot3D .p3d or .xyz mesh file, extract block face connectivity,
and print block connections. Usage:

    julia Run_File3D.jl path/to/mesh.p3d

Requires the Plot3D.jl package.
"""

using Plot3D

# Parse command-line arguments
if length(ARGS) < 1
    println("Usage: julia Run_File3D.jl path/to/mesh.p3d")
    return
end

meshfile = ARGS[1]
if !isfile(meshfile)
    println("ERROR: File not found: $meshfile")
    return
end

# Read mesh file (supports .p3d and .xyz)
blocks = endswith(meshfile, ".xyz") ?
    read_plot3D_ascii(meshfile) :
    read_blocks(meshfile)

# Use the correct face extraction for 2D or 3D blocks
if typeof(blocks[1]) <: Plot3D.Block
    block_faces = [get_faces(b) for b in blocks]
end

# Print connections
connections = print_connections(block_faces)
println("Block connections: ", connections)

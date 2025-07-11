include("/home/maiaross/.julia/dev/Plot3D/src/Plot3D.jl") 
using .Plot3D

# Read blocks from a test file
blocks = read_blocks(joinpath(@__DIR__, "/home/maiaross/.julia/dev/Cygnus/test/integration/2d/double_mach_reflection/DMR.p3d"))

# Convert blocks to Block2DFaces
block_faces = [Block2DFaces(b) for b in blocks]

# Print connections
connections = print_connections(block_faces)
println("2D Block connections: ", connections)
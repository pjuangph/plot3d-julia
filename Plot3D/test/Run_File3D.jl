include("/home/maiaross/JuliaGMSH/Plot3D/src/Plot3D.jl")
using .Plot3D

# Read structured XYZ file
# blocks = read_structured_xyz("/home/maiaross/JuliaGMSH/Plot3D/test/mesh.xyz")

# Normal - read .p3d file
blocks = read_blocks("/home/maiaross/.julia/dev/Plot3D/test/FewerBlocks.p3d")

# Convert blocks to Block2DFaces
block_faces = [Plot3D.Faces.Block2DFaces(b) for b in blocks]

# Print connections
connections = print_connections(block_faces)
println("2D Block connections: ", connections)
using Plot3D

# Path to your mesh file (.p3d or .xyz)
meshfile = "/home/maiaross/Plot3D/data/Copied_C_grid.p3d"  # Change to your actual file

# Read blocks from the mesh file
blocks = read_blocks(meshfile)  # For .xyz files, use read_structured_xyz(meshfile)

# Generate block faces (for 2D meshes; use Block3DFaces for 3D)
block_faces = [Plot3D.Faces.Block2DFaces(b) for b in blocks]

# Print block connections
connections = print_connections(block_faces)
println("Block connections: ", connections)

# Optionally, inspect a block
println("First block X coordinates:")
println(blocks[1].X[1])
# println("First block Y coordinates:")
# println(blocks[1].Y)

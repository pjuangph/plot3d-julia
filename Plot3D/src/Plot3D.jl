module Plot3D    
    using Parameters, IterTools, Printf
    using ProgressMeter
    using DataFrames
    import Base

    include("/home/maiaross/JuliaGMSH/Plot3D/src/smilBlocks.jl")
    using .BlockType
    export Block2D, Block3D, read_blocks

    #FOR 2D ONLY
    # include("2D_Match.jl")
    # using .Faces
    # export read_blocks, Block2DFaces, Block2D, Block3DFaces, Block3D, print_connections, block_faces

    #FOR .xyz files
    include("/home/maiaross/JuliaGMSH/Plot3D/src/xyz_readin.jl")
    using .xyzreader
    export Block3DXYZ, read_structured_xyz, read_xyz_block

    # FOR 3D
    include("3D_Match.jl")
    using .Faces
    export read_blocks, Block2DFaces, Block3DFaces, Block2D, Block3D, print_connections, block_faces,
    FaceND, forward_match, reverse_match, final_match_check, block_xyz_faces




end

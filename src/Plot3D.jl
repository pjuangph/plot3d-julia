module Plot3D    
    using Parameters, IterTools, Printf
    using ProgressMeter
    using DataFrames
    import Base

    # Load block structure
    include("Block.jl")
    
    # Load .xyz reader
    include("xyz_Reader.jl")
    using .xyzreader

    # Load face-matching utilities
    include("Face.jl")
    using .Faces



    export Block2D, Block3D, read_blocks, Block3DXYZ,
           read_structured_xyz, read_xyz_block,
           Block2DFaces, Block3DFaces, print_connections,
           block_faces, FaceND, forward_match,
           reverse_match, final_match_check, block_xyz_faces

end

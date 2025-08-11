module Plot3D

# =============== Explicit imports (alphabetical) ===============
import IterTools: product
import LinearAlgebra: dot, norm
import Parameters: @unpack
import Printf: @printf
import ProgressMeter: @showprogress
import Statistics: mean

# =============== Includes (load order matters) =================
include("block3D.jl")          # module Block3D
include("blockfunctions.jl")   # module BlockFunctions
include("face.jl")             # module Face3D
include("utils.jl")            # module Utils  <-- needed by facefunctions/connectivity
include("facefunctions.jl")    # functions live directly in Plot3D
include("connectivity.jl")     # functions live directly in Plot3D
include("readplot3d.jl")       # module ReadPlot3D
include("xyz_reader.jl")       # module XYZReader
include("writeplot3d.jl")      # module WritePlot3D

# =============== # Core types =================================
export Block
using .Block3D: Block

export reduce_blocks
using .BlockFunctions: reduce_blocks

# =============== # IO ==========================================
export read_plot3D_ascii, read_blocks, read_plot3D_binary, read_xyz_ascii
using .ReadPlot3D: read_plot3D_ascii, read_blocks, read_plot3D_binary
using .XYZReader:  read_xyz_ascii

export write_plot3D
using .WritePlot3D: write_plot3D

# =============== # Face core ===================================
export Face, add_vertex, vertices_equals, index_equals, match_indices,
       normal, to_dict, set_block_index, set_face_id, is_edge
using .Face3D: Face, add_vertex, vertices_equals, index_equals, match_indices,
               normal, to_dict, set_block_index, set_face_id, is_edge

# =============== # Face algorithms =============================
# (defined directly by include("facefunctions.jl"))
export get_faces, faces_match, find_matching_faces,
       get_outer_faces, create_face_from_diagonals, find_connected_faces,
       find_closest_block, find_bounding_faces, split_face,
       find_face_nearest_point, outer_face_dict_to_list,
       match_faces_dict_to_list, face_matches_to_dict

# =============== # Connectivity ================================
# (defined directly by include("connectivity.jl"))
export FaceMatchSet, point_match, select_multi_dimensional,
       find_matching_blocks, combinations_of_nearest_blocks,
       get_face_intersection

# =============== # Utils / Math ================================
export unique_pairs, ensure3d
using .Utils: unique_pairs, ensure3d

end # module

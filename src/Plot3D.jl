module Plot3D

# =============== Explicit imports (alphabetical) ===============
import IterTools: product
import LinearAlgebra: dot, norm, cross
import Parameters: @unpack
import Printf: @printf
import ProgressMeter: @showprogress
import Statistics: mean

# =============== Includes (load order matters) =================
include("block3D.jl")          # module Block3D
include("face.jl")             # module Face3D
include("utils.jl")            # module Utils  <-- needed by facefunctions/connectivity
include("block_face_functions.jl")     # reduce_blocks + face algorithms in Plot3D scope
include("connectivity.jl")     # functions live directly in Plot3D
include("readplot3d.jl")       # module ReadPlot3D
include("xyz_reader.jl")       # module XYZReader
include("writeplot3d.jl")      # module WritePlot3D
include("periodicity.jl")       # functions live directly in Plot3D
include("differencing.jl")
include("graph.jl")
include("split_block.jl")

include("gridpro.jl")               # module GridPro
include("glennht_classes.jl")       # module GlennHTClasses
include("glennht_export.jl")        # module glennht_export
include("glennht_import.jl")        # module GlennHTImport


# =============== # Core types =================================
export Block
using .Block3D: Block, recompute_centroid!


# =============== # IO ==========================================
export read_plot3D_ascii, read_blocks, read_plot3D_binary, read_xyz_ascii
using .ReadPlot3D: read_plot3D_ascii, read_blocks, read_plot3D_binary 
using .XYZReader:  read_xyz_ascii

# =============== # Face core ===================================
export Face, add_vertex, vertices_equals, index_equals, match_indices,
       normal, to_dict, set_block_index, set_face_id, is_edge
using .Face3D: Face, add_vertex, vertices_equals, index_equals, match_indices,
               normal, to_dict, set_block_index, set_face_id, is_edge

# =============== # Face algorithms =============================
# (defined directly by include("block_face_functions.jl"))
export get_faces, faces_match, find_matching_faces,
       get_outer_faces, create_face_from_diagonals, find_connected_faces,
       find_closest_block, find_bounding_faces, split_face,
       find_face_nearest_point, outer_face_dict_to_list,
       match_faces_dict_to_list, face_matches_to_dict

export reduce_blocks
# =============== # Connectivity ================================
# (defined directly by include("connectivity.jl"))
using .Connectivity: FaceMatchSet, point_match, select_multi_dimensional,
       find_matching_blocks, combinations_of_nearest_blocks,
       get_face_intersection, connectivity_fast, block_connection_matrix
export FaceMatchSet, point_match, select_multi_dimensional,
       find_matching_blocks, combinations_of_nearest_blocks,
       get_face_intersection, connectivity_fast, block_connection_matrix

# =============== # Periodicity ================================
export create_rotation_matrix, linear_real_transform,
       periodicity, periodicity_fast, rotated_periodicity, translational_periodicity

# =============== # Differencing ================================
export find_edges, find_face_edges

# =============== # Graph ================================
export add_connectivity_to_graph, block_to_graph, get_face_vertex_indices, get_starting_vertex


# =============== # Split Block ================================
export Direction, Direction_i, Direction_j, Direction_k, max_aspect_ratio, split_blocks


# =============== # GridPro ================================
export read_gridpro_to_blocks, read_gridpro_connectivity

# =============== # GlennHT I/O + helpers ================================
export BoundaryConditionType, InletBC_Subtype, InletBC_Direction, OutletSubtype,
       SymmetricSlipSubtype, WallSubtype, GIFCoordinate, GIFType, GIFOrder, TbModelType,
       BoundaryCondition, InletBC, OutletBC, SymmetricSlipBC, WallBC, GIF,
       BCGroup, JobFiles, JobControl, TurbModelInput, Plot3DParameters, InitialCond,
       TimeStpControl, SPDSchemeControl, RKSchemeControl, MGSchemeControl, GasPropertiesInput,
       ReferenceCond, ReferenceCondFull, Job

export to_pa, ideal_R, mach_from_p0_over_p, T_from_T0, mu_suth, a_sound, cp_from_gamma_R,
       pr_and_k, rpm_to_omegab, populate_reference_from_inputs,
       export_to_boundary_condition, export_to_job_file, summarize_contiguous,
       export_to_glennht_conn, read_ght_conn


# =============== # Utils / Math ================================
export unique_pairs, ensure3d
using .Utils: unique_pairs, ensure3d

end # module
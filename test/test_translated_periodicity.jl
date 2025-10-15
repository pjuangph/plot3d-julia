# test/manual/translated_periodicity_debug.jl
# A dead-simple, print-heavy debug script (no Test / no @testset)

using Downloads
using Dates
using Printf

# --- resolve repo paths -------------------------------------------------------
# This file lives in test/manual/, so src is two levels up.
SRC_DIR  = normpath(joinpath(@__DIR__, "..", "src"))
include(joinpath(SRC_DIR, "plot3d.jl"))
using .Plot3D: read_plot3D_binary, connectivity_fast, Block, match_faces_dict_to_list, outer_face_dict_to_list

# --- download helper ----------------------------------------------------------
function ensure_download(url::AbstractString, dest::AbstractString; force::Bool=false)
    if force || !isfile(dest)
        mkpath(dirname(dest))
        println("[$(Dates.now())] Downloading:\n  $url\n→ $dest")
        Downloads.download(url, dest)
    else
        println("[$(Dates.now())] File already present:\n  $dest")
    end
    isfile(dest) || error("File not found after download: $(abspath(dest))")
    sz = filesize(dest)
    sz > 0 || error("Downloaded file is empty: $(abspath(dest))")
    return sz
end

# --- main ---------------------------------------------------------------------
function main(; force_download::Bool=false)
    println("pwd = ", pwd())
    println("SRC_DIR = ", SRC_DIR)

    # put assets under test/manual/data to keep things tidy
    data_dir = joinpath(@__DIR__, "data")
    mkpath(data_dir)
    url  = "https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/iso65_64blocks.xyz"
    path = joinpath(data_dir, "iso65_64blocks.xyz")
    sz = ensure_download(url, path; force=force_download)
    println("exists? ", isfile(path), "   size(bytes) = ", sz)

    # 1) Read blocks (binary, float32)
    blocks = read_plot3D_binary(path; T=Float32, big_endian=false)
    println("blocks read = ", length(blocks))
    isempty(blocks) && error("No blocks were read.")

    # quick shape sanity for first block
    b1 = blocks[1]
    println("First block dims: IMAX=$(b1.IMAX) JMAX=$(b1.JMAX) KMAX=$(b1.KMAX)")
    size(b1.X) == (b1.IMAX, b1.JMAX, b1.KMAX) || error("X size mismatch")
    size(b1.Y) == (b1.IMAX, b1.JMAX, b1.KMAX) || error("Y size mismatch")
    size(b1.Z) == (b1.IMAX, b1.JMAX, b1.KMAX) || error("Z size mismatch")
    
    # 2) Connectivity: face_matches, outer_faces
    println("\nFinding connectivity …")
    methods(connectivity_fast) |> println
    typeof(blocks) |> println  # should be Vector{Plot3D.Block} (aka .Block3D.Block)

    face_matches, outer_faces = connectivity_fast(blocks)
    @printf("Got %d face_matches and %d outer_faces (dicts)\n",
            length(face_matches), length(outer_faces))

    # Some Python examples stored a 'match' table inside each match dict; drop if present
    for m in face_matches
        haskey(m, "match") && delete!(m, "match")
    end

    # 3) Build all_faces (as Face objects), then to dicts (like Python)
    println("Organizing split faces + outer faces into all_faces …")
    faces_from_matches = match_faces_dict_to_list(blocks, face_matches)   # -> Vector{Face}
    faces_from_outer   = outer_face_dict_to_list(blocks, outer_faces)     # -> Vector{Face}
    all_faces_faces    = vcat(faces_from_matches, faces_from_outer)       # Faces we’ll convert to dicts

    # (Optional) unique-ify faces to avoid duplicates
    let seen = Set{Tuple{Int,Int,Int,Int,Int,Int,Int}}()
        all_faces_faces = [f for f in all_faces_faces if begin
            key = (f.BlockIndex, f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX)
            if key in seen
                false
            else
                push!(seen, key); true
            end
        end]
    end

    # Convert to dicts so our APIs that expect dicts can use them
    all_faces = [to_dict(f) for f in all_faces_faces]
    @printf("all_faces: %d total (dicts)\n", length(all_faces))

    # 4) Block-to-block connection matrix (like Python)
    println("\nCreating block connection matrix …")
    C = block_connection_matrix(blocks, all_faces)  # returns a matrix/array
    @printf("Connection matrix size: %s\n", size(C))

    # 5) Translated periodicity
    println("\n=== Translational periodicity ===")
    # Python flow:
    # z_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks, all_faces, translational_direction='z')
    # x_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks, outer_faces, translational_direction='x')
    # y_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks, outer_faces, translational_direction='y')
    #
    # Julia version of translational_periodicity returns:
    # (periodic_faces_export::Vector{Dict}, periodic_pairs::Vector{Tuple{Face,Face,Dict}}, outer_faces_remaining::Vector{Dict})

    println("--> Z direction (using all_faces as candidates)")
    z_export, _z_pairs, outer_after_z = translational_periodicity(
        blocks, all_faces; translational_direction="z"
    )
    @printf("  z_export=%d   outer_after_z=%d\n", length(z_export), length(outer_after_z))

    println("--> X direction (use remaining outer faces from Z)")
    x_export, _x_pairs, outer_after_x = translational_periodicity(
        blocks, outer_after_z; translational_direction="x"
    )
    @printf("  x_export=%d   outer_after_x=%d\n", length(x_export), length(outer_after_x))

    println("--> Y direction (use remaining outer faces from X)")
    y_export, _y_pairs, outer_after_y = translational_periodicity(
        blocks, outer_after_x; translational_direction="y"
    )
    @printf("  y_export=%d   outer_after_y=%d\n", length(y_export), length(outer_after_y))

    # 6) Combine matched faces like Python did, then filter remaining outer faces
    println("\nCombining periodic matches (x+y+z) with original face_matches …")
    matched_faces = vcat(x_export, y_export, z_export, face_matches)   # dicts

    # Rebuild lists for set arithmetic on Face objects (unique, etc.)
    outer_list = outer_face_dict_to_list(blocks, outer_after_y)        # -> Faces
    matched_list = match_faces_dict_to_list(blocks, matched_faces)     # -> Faces

    # Unique sets
    matched_set = Set(matched_list)
    outer_filtered = [o for o in outer_list if !(o in matched_set)]

    @printf("Number of outer_faces remaining after filtering: %d (expected 0 in the tutorial case)\n",
            length(outer_filtered))

    # If you want the final outer_faces back as dicts:
    outer_faces_final = [to_dict(o) for o in outer_filtered]

    println("\n=== Summary ===")
    @printf("x_export: %d, y_export: %d, z_export: %d\n", length(x_export), length(y_export), length(z_export))
    @printf("all_periodic (export dicts) total: %d\n", length(x_export) + length(y_export) + length(z_export))
    @printf("outer_faces_final (dicts): %d\n", length(outer_faces_final))

    return (; all_faces, C, x_export, y_export, z_export, outer_faces_final)
end


main()  # set force_download=true to re-fetch

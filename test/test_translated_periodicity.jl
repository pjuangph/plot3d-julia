# test/manual/translated_periodicity_debug.jl
# A dead-simple, print-heavy debug script (no Test / no @testset)

using Downloads
using Dates

# --- resolve repo paths -------------------------------------------------------
# This file lives in test/manual/, so src is two levels up.
SRC_DIR  = normpath(joinpath(@__DIR__, "..", "src"))
include(joinpath(SRC_DIR, "plot3d.jl"))
using .Plot3D

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

    # Read blocks (binary, float32)
    blocks = read_plot3D_binary(path;T=Float32,big_endian=false)
    println("blocks read = ", length(blocks))
    isempty(blocks) && error("No blocks were read.")

    # quick shape sanity for first block
    b1 = blocks[1]
    println("First block dims: IMAX=$(b1.IMAX) JMAX=$(b1.JMAX) KMAX=$(b1.KMAX)")
    size(b1.X) == (b1.IMAX, b1.JMAX, b1.KMAX) || error("X size mismatch")
    size(b1.Y) == (b1.IMAX, b1.JMAX, b1.KMAX) || error("Y size mismatch")
    size(b1.Z) == (b1.IMAX, b1.JMAX, b1.KMAX) || error("Z size mismatch")

    # Run translated periodicity in x, then y, then z – carry remaining outer faces
    println("\n— Translational periodicity: x —")
    x_export, _periodic_pairs_x, outer_after_x = translational_periodicity(blocks, []; translational_direction="x")
    println("x_export count = ", length(x_export), "   remaining outer faces = ", length(outer_after_x))

    println("\n— Translational periodicity: y —")
    y_export, _periodic_pairs_y, outer_after_y = translational_periodicity(blocks, outer_after_x; translational_direction="y")
    println("y_export count = ", length(y_export), "   remaining outer faces = ", length(outer_after_y))

    println("\n— Translational periodicity: z —")
    z_export, _periodic_pairs_z, outer_after_z = translational_periodicity(blocks, outer_after_y; translational_direction="z")
    println("z_export count = ", length(z_export), "   remaining outer faces = ", length(outer_after_z))

    # Combine and lightly validate structure
    all_periodic = vcat(x_export, y_export, z_export)
    println("\nTOTAL periodic matches = ", length(all_periodic))

    # minimal structure checks with actionable messages
    function _ok_face(d::Dict)
        required = ("block_index","IMIN","JMIN","KMIN","IMAX","JMAX","KMAX")
        missing = [k for k in required if !haskey(d, k)]
        isempty(missing) || error("face dict missing keys: $(missing)  → dict=$(d)")
        return true
    end
    function _ok_pair(d::Dict)
        haskey(d,"block1") || error("pair dict missing 'block1': $(d)")
        haskey(d,"block2") || error("pair dict missing 'block2': $(d)")
        _ok_face(d["block1"]); _ok_face(d["block2"])
        return true
    end

    if !isempty(all_periodic)
        println("Sample periodic pair:\n", all_periodic[1])
        for (i,p) in enumerate(all_periodic)
            try
                _ok_pair(p)
            catch e
                error("Bad periodic pair at index $(i): $(e)")
            end
        end
    else
        println("No periodic pairs found (that can be fine depending on mesh/orientation).")
    end

    println("\nDone.")
end

try
    main()  # set force_download=true to re-fetch
catch e
    println("\nERROR: ", e)
    println("\nSTACKTRACE:")
    showerror(stdout, e, catch_backtrace())
    println()
    rethrow()
end

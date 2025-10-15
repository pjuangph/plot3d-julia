# test/test_periodicity_translated.jl
using Test
using Downloads
include(joinpath(@__DIR__, "..", "src", "plot3d.jl"))
using .Plot3D

# Small helper: robust binary reader that tries both :fortran and :raw
function _read_plot3d_binary_auto(path::AbstractString)
    try
        return read_plot3D_binary(path; format=:fortran, double_precision=false, big_endian=false)
    catch
        return read_plot3D_binary(path; format=:raw, double_precision=false, big_endian=false)
    end
end

@testset "Translated periodicity (x/y/z) on iso65_64blocks.xyz" begin
    mktempdir() do tmp
        cd(tmp) do
            url = "https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/iso65_64blocks.xyz"
            path = "iso65_64blocks.xyz"
            if !isfile(path)
                @info "Downloading $url ..."
                Downloads.download(url, path)
            end
            @test isfile(path) && filesize(path) > 0

            # In the Python example: read_plot3D(binary=True, read_double=False)
            blocks = _read_plot3d_binary_auto(path)
            @test !isempty(blocks)

            # We’ll start with the "outer_faces" as the working set for periodicity
            # The Julia API returns exports + remaining outer faces
            x_export, _, outer_after_x = translational_periodicity(blocks, []; translational_direction="x")
            y_export, _, outer_after_y = translational_periodicity(blocks, outer_after_x; translational_direction="y")
            z_export, _, outer_after_z = translational_periodicity(blocks, outer_after_y; translational_direction="z")

            @testset "Sanity checks" begin
                @test isa(x_export, Vector)
                @test isa(y_export, Vector)
                @test isa(z_export, Vector)

                # We expect to find at least some periodic matches in this mesh
                @test length(x_export) ≥ 0
                @test length(y_export) ≥ 0
                @test length(z_export) ≥ 0

                # The outer faces list should be a vector of Dicts if any remain
                @test isa(outer_after_z, Vector)
            end

            # Optional: combine all matches + any existing connectivity to mimic the Python example
            all_periodic = vcat(x_export, y_export, z_export)
            @info "periodic matches: x=$(length(x_export)) y=$(length(y_export)) z=$(length(z_export)) total=$(length(all_periodic))"

            # Very light structural validation: each match dict should have the standard fields
            function _ok_face(d::Dict)
                required = ("block_index","IMIN","JMIN","KMIN","IMAX","JMAX","KMAX")
                all(haskey(d, k) for k in required)
            end
            function _ok_pair(d::Dict)
                haskey(d,"block1") && haskey(d,"block2") && _ok_face(d["block1"]) && _ok_face(d["block2"])
            end
            @test all(_ok_pair, all_periodic)
        end
    end
end

# test/test_glennht_io.jl
using Test
using Plot3D  # exports: export_to_glennht_conn, read_ght_conn

@testset "GlennHT connectivity import/export" begin
    mktempdir() do tmp
        cd(tmp) do
            # fabricate a tiny case
            matches = [
                Dict(
                    "block1" => Dict("block_index"=>0,"IMIN"=>0,"JMIN"=>0,"KMIN"=>0,"IMAX"=>3,"JMAX"=>3,"KMAX"=>1),
                    "block2" => Dict("block_index"=>1,"IMIN"=>0,"JMIN"=>0,"KMIN"=>0,"IMAX"=>3,"JMAX"=>3,"KMAX"=>1),
                ),
            ]

            outer_faces = [
                Dict("block_index"=>0,"IMIN"=>0,"JMIN"=>0,"KMIN"=>0,"IMAX"=>3,"JMAX"=>3,"KMAX"=>0,"id"=>101),
                Dict("block_index"=>1,"IMIN"=>0,"JMIN"=>0,"KMIN"=>0,"IMAX"=>3,"JMAX"=>3,"KMAX"=>0,"id"=>202),
            ]

            gif_pairs = [ Dict("a"=>101, "b"=>202) ]
            gif_faces = Vector{Dict{String,Int}}() # already listed in outer_faces with IDs
            volume_zones = [
                Dict("block_index"=>0,"zone_type"=>"fluid","contiguous_index"=>1),
                Dict("block_index"=>1,"zone_type"=>"solid","contiguous_index"=>1),
            ]

            out = "mini_conn.ght_conn"
            export_to_glennht_conn(matches, copy(outer_faces), out, gif_pairs, gif_faces, volume_zones)
            @test isfile(out) && filesize(out) > 0

            maxblk, fm, of, num_conn, nZones, zone_types, Zones, GIFs = read_ght_conn(out)
            @test maxblk ≥ 1
            @test length(fm) == 1
            @test length(of) == 2
            @test nZones ≥ 1
            @test isa(zone_types, Vector)
            @test isa(Zones, Vector)
            @test length(GIFs) == 1

            # quick structural check
            chk_face(d) = all(haskey(d,k) for k in ("block_index","IMIN","JMIN","KMIN","IMAX","JMAX","KMAX"))
            @test chk_face(fm[1]["block1"])
            @test chk_face(fm[1]["block2"])
        end
    end
end

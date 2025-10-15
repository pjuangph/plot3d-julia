# test/test_gridpro_io.jl
using Test
using Plot3D  # exports: read_gridpro_connectivity

const GRIDPRO_SAMPLE = """
# Minimal GridPro-like connectivity
SB  1  something something  1   # pty near end → 1 (fluid)
SB  2  something something  2   # pty near end → 2 (solid)

# Patch lines:
# P pid sb1 sf1 sb2 sf2 fmap  L1i L1j L1k H1i H1j H1k  L2i L2j L2k H2i H2j H2k  pty lbid

# connection (pty=1), 0-based indices in file
P  10  1 0   2 0  012   0 0 0   3 3 0    0 0 0   3 3 0   1  999

# periodic (pty=3)
P  20  1 5   2 7  012   0 0 0   3 3 0    0 0 0   3 3 0   3  1000

# outer face (sb2=0 in file → -1 after normalization), keep pty not in [6,5,4,2]
P  30  1 1   0 0  012   0 0 0   3 3 0    0 0 0   0 0 0   7  42
"""

@testset "GridPro connectivity parser (smoke)" begin
    mktempdir() do tmp
        cd(tmp) do
            fn = "mini_gridpro.conn"
            open(fn, "w") do io
                write(io, GRIDPRO_SAMPLE)
            end

            result = read_gridpro_connectivity(fn; sb_zero_based_in_file=true, index_zero_based_in_file=true)
            @test isa(result, Dict)

            @test haskey(result, "face_matches")
            @test haskey(result, "outer_faces")
            @test haskey(result, "bc_group")
            @test haskey(result, "gif_faces")
            @test haskey(result, "periodic_faces")
            @test haskey(result, "volume_zones")

            conns = result["face_matches"]; @test length(conns) == 1
            per   = result["periodic_faces"]; @test length(per) == 1
            outer = result["outer_faces"];    @test length(outer) == 1

            bc = result["bc_group"]; @test haskey(bc, "inlet") && haskey(bc, "outlet") && haskey(bc, "symm_slip") && haskey(bc, "wall")

            # volume zones (SB lines → contiguous groups by odd/even pty)
            vz = result["volume_zones"]; @test !isempty(vz)
            @test all(haskey.(vz, "zone_type"))
        end
    end
end

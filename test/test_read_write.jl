using Test
using Downloads
using Plot3D

@testset "Plot3D I/O roundtrip (ASCII ⇄ Fortran binary, Float32)" begin
    # Work in an isolated temp dir so CI doesn't litter the repo
    mktempdir() do tmp
        cd(tmp) do
            # --- fetch test asset ---
            url = "https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz"
            ascii_path = "VSPT_ASCII.xyz"
            if !isfile(ascii_path)
                @info "Downloading $url ..."
                Downloads.download(url, ascii_path)
            end
            @test isfile(ascii_path)
            @test filesize(ascii_path) > 0

            # --- read ASCII ---
            blocks = read_plot3D_ascii(ascii_path)
            @info "Read $(length(blocks)) blocks from ASCII"
            @test !isempty(blocks)

            # basic integrity checks for each block
            @testset "ASCII block shapes" begin
                for b in blocks
                    @test size(b.X) == (b.IMAX, b.JMAX, b.KMAX)
                    @test size(b.Y) == (b.IMAX, b.JMAX, b.KMAX)
                    @test size(b.Z) == (b.IMAX, b.JMAX, b.KMAX)
                end
            end

            # --- write Fortran-unformatted Float32 binary ---
            bin_path = "VSPT_BINARY.xyzb"
            write_plot3D(
                bin_path, blocks;
                binary=true,
                format=:fortran,
                double_precision=false,
                big_endian=false,
            )
            @test isfile(bin_path)
            @info "Wrote binary: $(bin_path) size=$(filesize(bin_path)) bytes"

            # --- read back binary and verify shapes & values ---
            blocks_bin = read_plot3D_binary(
                bin_path; format=:fortran, double_precision=false, big_endian=false
            )
            @info "Read $(length(blocks_bin)) blocks from binary"
            @test length(blocks_bin) == length(blocks)

            @testset "Binary block shapes" begin
                for (a, b) in zip(blocks, blocks_bin)
                    @test (b.IMAX, b.JMAX, b.KMAX) == (a.IMAX, a.JMAX, a.KMAX)
                    @test size(b.X) == size(a.X)
                    @test size(b.Y) == size(a.Y)
                    @test size(b.Z) == size(a.Z)
                end
            end

            # Value equivalence (bitwise for Float32 I/O)
            @testset "Data equality after roundtrip" begin
                for (a, b) in zip(blocks, blocks_bin)
                    @test a.X == b.X
                    @test a.Y == b.Y
                    @test a.Z == b.Z
                end
            end

            # --- size sanity check for Fortran records (Float32) ---
            expected_bytes_fortran_float32(blocks) = begin
                nb = length(blocks)
                total = 0
                # record(nblocks: 1*UInt32)
                total += 4 + 4 + 4
                # dims per block (3*UInt32) as a record
                total += nb * (4 + 3*4 + 4)
                # payload per block: X,Y,Z records, each n*Float32 wrapped
                for b in blocks
                    n = b.IMAX*b.JMAX*b.KMAX
                    total += 3 * (4 + n*4 + 4)
                end
                total
            end

            expected = expected_bytes_fortran_float32(blocks)
            actual = filesize(bin_path)
            @info "Expected bytes (Float32 + Fortran records): $expected"
            @info "Actual bytes: $actual"
            @test actual == expected
        end
    end
end

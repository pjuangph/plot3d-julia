# using Test
# using Plot3D
using Downloads
include("../src/block3D.jl")
include("../src/readplot3D.jl")
include("../src/writeplot3d.jl")

# --- test starts here --------------------------------------------------------
url = "https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz"
ascii_path = "VSPT_ASCII.xyz"
if !isfile(ascii_path)
    @info "Downloading $url ..."
    Downloads.download(url, ascii_path)
end

# read ASCII using your package
blocks = ReadPlot3D.read_plot3D_ascii(ascii_path)
@info "Read $(length(blocks)) blocks from ASCII"

# basic integrity checks for each block
for b in blocks
    @assert size(b.X) == (b.IMAX, b.JMAX, b.KMAX)
    @assert size(b.Y) == (b.IMAX, b.JMAX, b.KMAX)
    @assert size(b.Z) == (b.IMAX, b.JMAX, b.KMAX)
end

# write binary plot3d (Fortran unformatted)
bin_path = "VSPT_BINARY.xyzb"
WritePlot3D.write_plot3D(bin_path, blocks;
             binary=true,
             format=:fortran,      # not ":="
             double_precision=false,
             big_endian=false)

blocks_binary = ReadPlot3D.read_plot3D_binary(bin_path; format=:fortran,
                             double_precision=false, big_endian=false)

@info "Wrote $(bin_path) size=$(filesize(bin_path)) bytes"
@info "Wrote $(length(blocks_binary)) blocks from binary"

# rough size sanity (not exact, but catches obvious mistakes)
# bytes = rec(nblocks) + sum[ rec(dims) + rec(X) + rec(Y) + rec(Z) ]
# where rec(payload) = 4 + sizeof(payload) + 4
function expected_bytes_fortran_float32(blocks)
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
    return total
end
expected = expected_bytes_fortran_float32(blocks)
@info "Expected bytes (Float32 + Fortran records): $expected"
@info "Actual bytes: $(filesize(bin_path))"
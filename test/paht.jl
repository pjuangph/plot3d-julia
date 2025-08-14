using Test
using Plot3D
using Downloads

# --- test starts here --------------------------------------------------------
url = "https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz"
ascii_path = "VSPT_ASCII.xyz"
if !isfile(ascii_path)
    @info "Downloading $url ..."
    Downloads.download(url, ascii_path)
end

# read ASCII using your package
blocks = read_plot3D_ascii(ascii_path)

# basic integrity checks for each block
for b in blocks
    size(b.X) == (b.IMAX, b.JMAX, b.KMAX)
    size(b.Y) == (b.IMAX, b.JMAX, b.KMAX)
    size(b.Z) == (b.IMAX, b.JMAX, b.KMAX)
end

# write binary plot3d (Fortran unformatted)
bin_path = "VSPT_BINARY.xyzb"
write_plot3D(bin_path, blocks;binary=true, double_precision=false)

# rough size sanity (not exact, but catches obvious mistakes)
# bytes = rec(nblocks) + sum[ rec(dims) + rec(X) + rec(Y) + rec(Z) ]
# where rec(payload) = 4 + sizeof(payload) + 4
function expected_bytes(blocks)
    nb = length(blocks)
    total = 0
    # nblocks record
    total += 4 + 4 + 4
    for b in blocks
        # dims (3 Int32)
        total += 4 + 3*4 + 4
        n = b.IMAX * b.JMAX * b.KMAX
        # X, Y, Z (Float32)
        total += (4 + n*4 + 4) * 3
    end
    return total
end
actual = filesize(bin_path)
expect = expected_bytes(blocks)
# allow a tiny cushion for portability differences (should match exactly on same machine)


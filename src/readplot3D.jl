module ReadPlot3D

using ..Block3D: Block

export read_plot3D_ascii, read_blocks, read_plot3D_binary

# ==============================================================
# Helpers
# ==============================================================

# Read one block (ASCII): dims line then X,Y,Z in that order
function _read_block_ascii!(io::IO, ::Type{T}) where {T<:Real}
    line = ""
    while !eof(io) && isempty(strip(line))
        line = readline(io)               # skip empties/comments space
    end
    eof(io) && error("Unexpected EOF while reading block dimensions")
    dims = split(strip(line))
    length(dims) == 3 || error("Expected 3 integers for I J K, got: $line")
    I = parse(Int, dims[1]); J = parse(Int, dims[2]); K = parse(Int, dims[3])

    n = I*J*K
    X = Array{T}(undef, I,J,K)
    Y = Array{T}(undef, I,J,K)
    Z = Array{T}(undef, I,J,K)

    # Read n numbers for each of X, Y, Z
    # Use a local buffer and fill! with Cartesian indexing
    function _fill3!(A::Array{T,3})
        filled = 0
        Ci = CartesianIndices(A)
        it = Iterators.take(Ci, n)
        for IJK in it
            # read numbers skipping blank lines
            val::Union{Nothing,T} = nothing
            while val === nothing
                eof(io) && error("Unexpected EOF reading coordinate values")
                line = strip(readline(io))
                isempty(line) && continue
                for tok in eachsplit(line)
                    A[IJK] = parse(T, tok)
                    filled += 1
                    IJK = iterate(Ci, IJK)[1]  # move to next IJK
                    if filled == n
                        return
                    end
                end
            end
        end
    end

    _fill3!(X); _fill3!(Y); _fill3!(Z)
    return Block(X,Y,Z)
end

# Unformatted Fortran read helpers (4-byte record markers)
# If record_markers=false, we read raw payload directly.
struct _Endian{B} end
const _BE = _Endian{:be}()
const _LE = _Endian{:le}()

@inline function _read_i32(io::IO, ::_Endian{:be})
    b = read(io, UInt32); return Int32(ntoh(b))   # ntoh = big-endian
end
@inline function _read_i32(io::IO, ::_Endian{:le})
    b = read(io, UInt32); return Int32(b)        # host assumed little; UInt32 read is LE
end

# read `N` reals of type T into a preallocated vector
function _read_reals!(io::IO, buf::Vector{T}) where {T<:Real}
    read!(io, reinterpret(UInt8, buf))  # raw bulk read
    return buf
end

# Read one block (binary): dims then X,Y,Z (each in its own record if markers=true)
function _read_block_binary!(io::IO, ::Type{T}; big_endian::Bool=true, record_markers::Bool=true) where {T<:Real}
    endian = big_endian ? _BE : _LE

    # dims: either in a record (with marker) or raw
    if record_markers
        reclen1 = _read_i32(io, endian)
        reclen1 == 12 || error("Expected 12-byte dims record, got $reclen1")
    end
    I = Int(read(io, Int32)); J = Int(read(io, Int32)); K = Int(read(io, Int32))
    if record_markers
        _ = _read_i32(io, endian)  # closing marker
    end

    n = I*J*K
    X = Array{T}(undef, I,J,K)
    Y = Array{T}(undef, I,J,K)
    Z = Array{T}(undef, I,J,K)

    # helper to read one record of n*T bytes (or raw if no markers)
    function _read_field!(A::Array{T,3})
        buf = Vector{T}(undef, length(A))
        if record_markers
            reclen = _read_i32(io, endian)
            reclen == sizeof(T)*length(A) || error("Record length mismatch: expected $(sizeof(T)*length(A)), got $reclen")
            _read_reals!(io, buf)
            _ = _read_i32(io, endian)
        else
            _read_reals!(io, buf)
        end
        # Fortran usually writes in the same linear order as Julia's column-major,
        # but some generators differ; if you observe transposed axes, swap here.
        A[:] .= buf
        return nothing
    end

    _read_field!(X); _read_field!(Y); _read_field!(Z)
    return Block(X,Y,Z)
end

# ==============================================================
# Public API
# ==============================================================

"""
    read_plot3D_ascii(path; T=Float64)

Read a multi-block Plot3D **grid** file in ASCII format.
File layout (common NASA variant):
- first non-empty line: integer `nb` (number of blocks)
- for each block: line with `I J K`
- then `I*J*K` reals for `X`, then same for `Y`, then same for `Z`.

Returns `Vector{Block}`.
"""
function read_plot3D_ascii(path::AbstractString; T::Type{<:Real}=Float64)
    open(path, "r") do io
        # skip empties/comments to find nb
        line = ""
        while !eof(io) && isempty(strip(line))
            line = readline(io)
        end
        eof(io) && error("Empty Plot3D ASCII file")
        nb = parse(Int, split(strip(line))[1])
        blocks = Vector{Block}(undef, nb)
        for b in 1:nb
            blocks[b] = _read_block_ascii!(io, T)
        end
        return blocks
    end
end

"""
    read_plot3D_binary(path; T=Float64, big_endian=true, record_markers=true)

Read a multi-block Plot3D **grid** file in unformatted Fortran binary.
Assumes:
- Leading record with `nb::Int32` (if `record_markers=true`)
- For each block:
    - dims record: `Int32 I, Int32 J, Int32 K`
    - record for `X` (I*J*K values of `T`)
    - record for `Y`
    - record for `Z`

If `record_markers=false`, we read raw payload in the same order without markers.
Returns `Vector{Block}`.
"""
function read_plot3D_binary(path::AbstractString; T::Type{<:Real}=Float64, big_endian::Bool=true, record_markers::Bool=true)
    open(path, "r") do io
        endian = big_endian ? _BE : _LE
        nb::Int
        if record_markers
            len1 = _read_i32(io, endian)
            len1 == 4 || error("Expected 4-byte record for nb, got $len1")
            nb = Int(read(io, Int32))
            _ = _read_i32(io, endian)
        else
            nb = Int(read(io, Int32))
        end
        blocks = Vector{Block}(undef, nb)
        for b in 1:nb
            blocks[b] = _read_block_binary!(io, T; big_endian=big_endian, record_markers=record_markers)
        end
        return blocks
    end
end

"""
    read_blocks(path; fmt=:ascii, kwargs...)

Convenience wrapper:
- `fmt = :ascii`  => `read_plot3D_ascii(path; kwargs...)`
- `fmt = :binary` => `read_plot3D_binary(path; kwargs...)`
"""
function read_blocks(path::AbstractString; fmt::Symbol=:ascii, kwargs...)
    if fmt === :ascii
        return read_plot3D_ascii(path; kwargs...)
    elseif fmt === :binary
        return read_plot3D_binary(path; kwargs...)
    else
        error("read_blocks: unknown fmt=$fmt (use :ascii or :binary)")
    end
end

end # module

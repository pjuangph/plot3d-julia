module WritePlot3D

import Printf: @printf
import ..Block3D: Block

# =========================
# Endianness / numeric I/O
# =========================

_unsigned(::Type{Int32})   = UInt32
_unsigned(::Type{UInt32})  = UInt32
_unsigned(::Type{Float32}) = UInt32
_unsigned(::Type{Float64}) = UInt64

const HOST_IS_LITTLE_ENDIAN = ENDIAN_BOM == UInt32(0x04030201)

"""
    write_num(io, x, big_endian)

Write primitive `x` to `io` with chosen endianness.
"""
function write_num(io::IO, x::T, big_endian::Bool) where {T<:Union{Int32,UInt32,Float32,Float64}}
    # if the requested endianness matches the host, write directly; else swap
    
    if big_endian == !HOST_IS_LITTLE_ENDIAN
        write(io, x)
    else
        u = reinterpret(_unsigned(T), x)
        write(io, bswap(u))
    end
end

# convenience for writing a whole vector of numbers with endianness
function write_vec(io::IO, v::AbstractVector{T}, big_endian::Bool) where {T<:Union{Int32,UInt32,Float32,Float64}}
    @inbounds for x in v
        write_num(io, x, big_endian)
    end
end

# =========================
# Block payload writers
# =========================

# i–j–k order (Plot3D convention), supports 2D via KMAX==1
function _write_block_binary_raw(io::IO, B::Block; double_precision::Bool=true, big_endian::Bool=false)
    Tout = double_precision ? Float64 : Float32

    @inline function write_var(V::Array{Float64,3})
        @inbounds for k in 1:B.KMAX
            for j in 1:B.JMAX
                for i in 1:B.IMAX
                    write_num(io, convert(Tout, V[i,j,k]), big_endian)
                end
            end
        end
        nothing
    end

    write_var(B.X); write_var(B.Y); write_var(B.Z)
    return nothing
end

# Fortran unformatted sequential: wrap each logical record with 4-byte record markers
# Header records:
#   [ nblocks::UInt32 ]   as one record
#   [ IMAX JMAX KMAX ]    one record per block (UInt32 each)
# Payload records:
#   X (i–j–k), Y (i–j–k), Z (i–j–k) — one record per array
function _write_record(io::IO, bytes::Vector{UInt8}, big_endian::Bool)
    len_u = UInt32(length(bytes))
    write_num(io, len_u, big_endian)
    write(io, bytes)
    write_num(io, len_u, big_endian)
end

function _pack_vec_bytes(v::AbstractVector{T}, big_endian::Bool) where {T<:Union{Int32,UInt32,Float32,Float64}}
    io = IOBuffer()
    write_vec(io, v, big_endian)
    take!(io)
end

function _write_block_binary_fortran(io::IO, B::Block; double_precision::Bool=true, big_endian::Bool=false)
    Tout = double_precision ? Float64 : Float32

    # flatten in i–j–k order
    function flatten_var(V::Array{Float64,3})
        out = Vector{Tout}(undef, B.IMAX*B.JMAX*B.KMAX)
        idx = 1
        @inbounds for k in 1:B.KMAX, j in 1:B.JMAX, i in 1:B.IMAX
            out[idx] = convert(Tout, V[i,j,k]); idx += 1
        end
        out
    end

    # X
    xb = _pack_vec_bytes(flatten_var(B.X), big_endian)
    _write_record(io, xb, big_endian)

    # Y
    yb = _pack_vec_bytes(flatten_var(B.Y), big_endian)
    _write_record(io, yb, big_endian)

    # Z
    zb = _pack_vec_bytes(flatten_var(B.Z), big_endian)
    _write_record(io, zb, big_endian)

    return nothing
end

# =========================
# Public API
# =========================

"""
    write_plot3D(filename, blocks;
                 binary::Bool=true,
                 format::Symbol = :fortran,   # :fortran (default) or :raw
                 double_precision::Bool=true,
                 big_endian::Bool=false,
                 columns::Int=6)

Write a Plot3D multi-block file from a vector of `Block`s.

**Binary** (when `binary=true`):
- `format = :fortran` (**default**): Fortran unformatted sequential with 4-byte record markers.
  *Header records*:
    1) `[ UInt32(nblocks) ]`
    2) for each block: `[ UInt32(IMAX), UInt32(JMAX), UInt32(KMAX) ]`
  *Payload*: three records per block: X, Y, Z in **i–j–k** order.
  `double_precision` selects `Float64` vs `Float32`. `big_endian` controls byte order.

- `format = :raw`: raw header then raw payload, no record markers (the previous behavior).

**ASCII** (when `binary=false`): same layout as before, `columns` controls wrapping.

2D is supported by setting `KMAX==1`.
"""
function write_plot3D(filename::AbstractString, blocks::Vector{Block};
                      binary::Bool=true,
                      format::Symbol = :fortran,
                      double_precision::Bool=true,
                      big_endian::Bool=false,
                      columns::Int=6)

    if binary
        open(filename, "w") do io
            if format === :raw
                # RAW HEADER
                write_num(io, UInt32(length(blocks)), big_endian)
                @inbounds for b in blocks
                    write_num(io, UInt32(b.IMAX), big_endian)
                    write_num(io, UInt32(b.JMAX), big_endian)
                    write_num(io, UInt32(b.KMAX), big_endian)
                end
                # RAW PAYLOAD
                @inbounds for b in blocks
                    _write_block_binary_raw(io, b; double_precision=double_precision, big_endian=big_endian)
                end

            elseif format === :fortran
                # FORTRAN HEADER
                nblocks_rec = Vector{UInt32}(undef, 1); nblocks_rec[1] = UInt32(length(blocks))
                _write_record(io, _pack_vec_bytes(nblocks_rec, big_endian), big_endian)

                @inbounds for b in blocks
                    dims = UInt32[b.IMAX, b.JMAX, b.KMAX]
                    _write_record(io, _pack_vec_bytes(dims, big_endian), big_endian)
                end

                # FORTRAN PAYLOAD
                @inbounds for b in blocks
                    _write_block_binary_fortran(io, b; double_precision=double_precision, big_endian=big_endian)
                end
            else
                error("write_plot3D: unknown binary format $(format). Use :fortran or :raw.")
            end
        end

    else
        # ASCII
        open(filename, "w") do io
            println(io, length(blocks))
            @inbounds for b in blocks
                @printf(io, "%d %d %d\n", b.IMAX, b.JMAX, b.KMAX)
            end
            @inline function write_var(V::Array{Float64,3})
                col = 0
                @inbounds for k in 1:blocks[1].KMAX
                    for j in 1:blocks[1].JMAX
                        for i in 1:blocks[1].IMAX
                            @printf(io, "%0.8f ", V[i,j,k])
                            col += 1
                            if (col % columns) == 0
                                write(io, '\n')
                            end
                        end
                    end
                end
                if (col % columns) != 0
                    write(io, '\n')
                end
            end
            @inbounds for b in blocks
                write_var(b.X); write_var(b.Y); write_var(b.Z)
            end
        end
    end
    return nothing
end

end # module

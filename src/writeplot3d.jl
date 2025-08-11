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

"""
    write_num(io, x, big_endian)

Write primitive `x` to `io` with chosen endianness.
"""
function write_num(io::IO, x::T, big_endian::Bool) where {T<:Union{Int32,UInt32,Float32,Float64}}
    if big_endian == Base.isbigendian()
        write(io, x)
    else
        u = reinterpret(_unsigned(T), x)
        write(io, bswap(u))
    end
end

# =========================
# Binary payload writer
# =========================
# i–j–k order (Plot3D convention), supports 2D via KMAX==1

function _write_plot3D_block_binary(io::IO, B::Block; double_precision::Bool=true, big_endian::Bool=false)
    Tout = double_precision ? Float64 : Float32

    @inline function write_var(V::Array{Float64,3})
        @inbounds for k in 1:B.KMAX
            for j in 1:B.JMAX
                for i in 1:B.IMAX
                    write_num(io, convert(Tout, V[i,j,k]), big_endian)
                end
            end
        end
        return nothing
    end

    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)
    return nothing
end

# =========================
# ASCII payload writer
# =========================

function _write_plot3D_block_ascii(io::IO, B::Block; columns::Int=6)
    @inline function write_var(V::Array{Float64,3})
        col = 0
        @inbounds for k in 1:B.KMAX
            for j in 1:B.JMAX
                for i in 1:B.IMAX
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
        return nothing
    end

    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)
    return nothing
end

# =========================
# Public API
# =========================

"""
    write_plot3D(filename, blocks;
                 binary::Bool=true,
                 double_precision::Bool=true,
                 big_endian::Bool=false,
                 columns::Int=6)

Write a Plot3D multi-block file from a vector of `Block`s.

- When `binary=true` (**default**): header is written as:
    - `UInt32(nblocks)`
    - then for each block: `UInt32(IMAX)`, `UInt32(JMAX)`, `UInt32(KMAX)`
  followed by X, Y, Z in **i–j–k** order. Use `double_precision=false` to write `Float32`.
  Use `big_endian=true` to force big-endian output.

- When `binary=false`: ASCII format with the same logical structure; `columns` controls line wrapping.

2D is supported by setting `KMAX==1`.
"""
function write_plot3D(filename::AbstractString, blocks::Vector{Block};
                      binary::Bool=true, double_precision::Bool=true,
                      big_endian::Bool=false, columns::Int=6)
    if binary
        open(filename, "w") do io
            # header
            write_num(io, UInt32(length(blocks)), big_endian)
            @inbounds for b in blocks
                write_num(io, UInt32(b.IMAX), big_endian)
                write_num(io, UInt32(b.JMAX), big_endian)
                write_num(io, UInt32(b.KMAX), big_endian)
            end
            # payload
            @inbounds for b in blocks
                _write_plot3D_block_binary(io, b; double_precision=double_precision, big_endian=big_endian)
            end
        end
    else
        open(filename, "w") do io
            # header
            println(io, length(blocks))
            @inbounds for b in blocks
                @printf(io, "%d %d %d\n", b.IMAX, b.JMAX, b.KMAX)
            end
            # payload
            @inbounds for b in blocks
                _write_plot3D_block_ascii(io, b; columns=columns)
            end
        end
    end
    return nothing
end

end # module

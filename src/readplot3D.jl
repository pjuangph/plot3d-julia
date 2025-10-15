module ReadPlot3D

using ..Block3D: Block

export read_plot3D_ascii, read_blocks, read_plot3D_binary

# ==============================================================
# Endian + binary helpers (RAW only; no Fortran record markers)
# ==============================================================

_unsigned(::Type{Int32})   = UInt32
_unsigned(::Type{UInt32})  = UInt32
_unsigned(::Type{Float32}) = UInt32
_unsigned(::Type{Float64}) = UInt64

const HOST_IS_LITTLE_ENDIAN = ENDIAN_BOM == UInt32(0x04030201)

@inline function _read_u32(io::IO, big_endian::Bool)::UInt32
    u = read(io, UInt32)
    return (big_endian == HOST_IS_LITTLE_ENDIAN) ? bswap(u) : u
end
@inline _read_i32(io::IO, big_endian::Bool)::Int32 = reinterpret(Int32, _read_u32(io, big_endian))

# Interpret a byte buffer as a vector of T, swapping endianness in-place if needed
function _bytes_as_vector_T!(buf::Vector{UInt8}, ::Type{T}, big_endian::Bool) where {T<:Union{Float32,Float64}}
    (length(buf) % sizeof(T) == 0) || error("Byte buffer length $(length(buf)) not multiple of sizeof($(T))=$(sizeof(T))")
    vecT = reinterpret(T, buf)  # view onto buf
    if big_endian == HOST_IS_LITTLE_ENDIAN
        uw = reinterpret(_unsigned(T), vecT)
        @inbounds for i in eachindex(uw)
            uw[i] = bswap(uw[i])
        end
    end
    return vecT
end

# Read n values of T (RAW; no record markers)
function _read_values_T(io::IO, ::Type{T}, n::Int, big_endian::Bool) where {T<:Union{Float32,Float64}}
    nbytes = n * sizeof(T)
    buf = Vector{UInt8}(undef, nbytes)
    read!(io, buf)
    return _bytes_as_vector_T!(buf, T, big_endian)
end

# ==============================================================
# ASCII reader
# ==============================================================

# Consume and return next non-empty line (keeps simple ASCII variant)
function _next_data_line(io::IO)
    while !eof(io)
        ln = strip(readline(io))
        if !isempty(ln)
            return ln
        end
    end
    return ""
end

# Read exactly n reals of type T from ASCII (can span multiple lines)
function _read_n_reals_ascii!(io::IO, n::Int, ::Type{T}) where {T<:Real}
    out = Vector{T}(undef, n)
    i = 1
    while i <= n
        eof(io) && error("Unexpected EOF while reading ASCII values (needed $n, got $(i-1))")
        ln = strip(readline(io))
        isempty(ln) && continue
        for tok in eachsplit(ln)
            out[i] = parse(T, tok)
            i += 1
            i > n && break
        end
    end
    return out
end

# Read one block (ASCII): dims line then X,Y,Z in i–j–k order
function _read_block_ascii!(io::IO, ::Type{T}) where {T<:Real}
    dimline = strip(_next_data_line(io))
    isempty(dimline) && error("Unexpected EOF while reading block dimensions")
    parts = split(dimline)
    length(parts) >= 3 || error("Expected 3 integers for I J K, got: '$dimline'")
    I = parse(Int, parts[1]); J = parse(Int, parts[2]); K = parse(Int, parts[3])
    (I>0 && J>0 && K>0) || error("Non-positive block dims: ($I,$J,$K)")

    n = I*J*K
    vx = _read_n_reals_ascii!(io, n, T)
    vy = _read_n_reals_ascii!(io, n, T)
    vz = _read_n_reals_ascii!(io, n, T)

    X = Array{Float64,3}(undef, I,J,K)
    Y = Array{Float64,3}(undef, I,J,K)
    Z = Array{Float64,3}(undef, I,J,K)
    idx = 1
    @inbounds for k in 1:K, j in 1:J, i in 1:I
        X[i,j,k] = Float64(vx[idx])
        Y[i,j,k] = Float64(vy[idx])
        Z[i,j,k] = Float64(vz[idx])
        idx += 1
    end
    return Block(X,Y,Z)
end

"""
    read_plot3D_ascii(path; T=Float64)

Read a multi-block Plot3D **grid** file in ASCII format:
- first non-empty line: integer `nb` (number of blocks)
- for each block: `I J K`, then `I*J*K` reals for X, then Y, then Z (i–j–k order).
Returns `Vector{Block}`.
"""
function read_plot3D_ascii(path::AbstractString; T::Type{<:Real}=Float64)
    open(path, "r") do io
        first = strip(_next_data_line(io))
        isempty(first) && error("Empty Plot3D ASCII file")
        nb = parse(Int, split(first)[1])
        (nb > 0) || error("Invalid number of blocks: $nb")
        blocks = Vector{Block}(undef, nb)
        @inbounds for b in 1:nb
            blocks[b] = _read_block_ascii!(io, T)
        end
        return blocks
    end
end

# ==============================================================
# RAW binary reader (no Fortran record markers)
# ==============================================================

# Read one block (RAW): dims then X,Y,Z (all contiguous)
function _read_block_binary_raw!(io::IO, ::Type{T}; big_endian::Bool=false) where {T<:Union{Float32,Float64}}
    I = Int(_read_i32(io, big_endian))
    J = Int(_read_i32(io, big_endian))
    K = Int(_read_i32(io, big_endian))
    (I>0 && J>0 && K>0) || error("Non-positive block dims: ($I,$J,$K)")

    n = I*J*K
    vX = _read_values_T(io, T, n, big_endian)
    vY = _read_values_T(io, T, n, big_endian)
    vZ = _read_values_T(io, T, n, big_endian)

    X = Array{Float64,3}(undef, I,J,K)
    Y = Array{Float64,3}(undef, I,J,K)
    Z = Array{Float64,3}(undef, I,J,K)
    idx = 1
    @inbounds for k in 1:K, j in 1:J, i in 1:I
        X[i,j,k] = Float64(vX[idx])
        Y[i,j,k] = Float64(vY[idx])
        Z[i,j,k] = Float64(vZ[idx])
        idx += 1
    end
    return Block(X,Y,Z)
end

"""
    read_plot3D_binary(path; T=Float32, big_endian=false)

Read a multi-block Plot3D **grid** file in RAW binary (no record markers).
Layout:
- `Int32 nb`
- For each block: `Int32 I, Int32 J, Int32 K`
- Then for each block: X, Y, Z — each `I*J*K` reals of type `T` in i–j–k order.

Returns `Vector{Block}`.
"""
function read_plot3D_binary(path::AbstractString;
    T::Type{<:Union{Float32,Float64}} = Float32,
    big_endian::Bool = false
)
    open(path, "r") do io
        # --- header: number of blocks ---
        nb = Int(_read_i32(io, big_endian))
        (nb > 0) || error("Invalid number of blocks: $nb")

        # --- header: per-block dimensions (peek+store so payload can be read later) ---
        dims = Vector{NTuple{3,Int}}(undef, nb)
        @inbounds for b in 1:nb
            I = Int(_read_i32(io, big_endian))
            J = Int(_read_i32(io, big_endian))
            K = Int(_read_i32(io, big_endian))
            (I>0 && J>0 && K>0) || error("Non-positive dims for block $b: ($I,$J,$K)")
            dims[b] = (I, J, K)
        end

        # --- payload: X,Y,Z per block ---
        blocks = Vector{Block}(undef, nb)
        @inbounds for b in 1:nb
            I, J, K = dims[b]
            n = I*J*K

            vX = _read_values_T(io, T, n, big_endian)
            vY = _read_values_T(io, T, n, big_endian)
            vZ = _read_values_T(io, T, n, big_endian)

            X = Array{Float64,3}(undef, I,J,K)
            Y = Array{Float64,3}(undef, I,J,K)
            Z = Array{Float64,3}(undef, I,J,K)

            idx = 1
            @inbounds for k in 1:K, j in 1:J, i in 1:I
                X[i,j,k] = Float64(vX[idx])
                Y[i,j,k] = Float64(vY[idx])
                Z[i,j,k] = Float64(vZ[idx])
                idx += 1
            end

            blocks[b] = Block(X, Y, Z)
        end

        return blocks
    end
end

# ==============================================================
# Unified front-end
# ==============================================================

"""
    read_blocks(path; fmt=:ascii, kwargs...)

Convenience wrapper:
- `fmt = :ascii`  => `read_plot3D_ascii(path; kwargs...)`
- `fmt = :binary` => `read_plot3D_binary(path; kwargs...)`  (RAW only)
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

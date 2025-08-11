module ReadPlot3D

import ..Block3D: Block

export read_plot3D_ascii, read_blocks, read_plot3D_binary

# ------------------------------------------------------------
# Internal: read exactly N Float64s from an IO of ASCII numbers
# ------------------------------------------------------------
function _read_n_ascii_floats!(io::IO, N::Int, out::Vector{Float64})
    empty!(out)
    sizehint!(out, N)
    for ln in eachline(io)
        isempty(ln) && continue
        for tok in split(ln)
            isempty(tok) && continue
            push!(out, parse(Float64, tok))
            length(out) == N && return nothing
        end
    end
    error("Unexpected EOF while reading ASCII Plot3D payload (wanted $N values).")
end

# ------------------------------------------------------------
# Public: read ASCII Plot3D multi-block file
# ------------------------------------------------------------
"""
    read_plot3D_ascii(path::AbstractString) -> Vector{Block}

Read an **ASCII** Plot3D file with layout:

nblocks
IMAX JMAX KMAX # per block (one line each)
... # Block 1 X (IMAXJMAXKMAX numbers, i–j–k order)
... # Block 1 Y
... # Block 1 Z
(repeat for each block)


Supports 2D blocks by setting `KMAX == 1`.
"""
function read_plot3D_ascii(path::AbstractString)
    open(path, "r") do io
        # -------- header: number of blocks --------
        firstline = nothing
        for ln in eachline(io)
            if !isempty(strip(ln))
                firstline = ln
                break
            end
        end
        firstline === nothing && error("Empty file: $path")
        nblocks = parse(Int, split(firstline)[1])

        # -------- per-block sizes --------
        dims = Vector{NTuple{3,Int}}(undef, nblocks)
        for b in 1:nblocks
            ln = ""
            while isempty(strip(ln))
                ln = readline(io)
            end
            toks = split(ln)
            length(toks) >= 3 || error("Bad size line for block $b: '$ln'")
            IMAX = parse(Int, toks[1]); JMAX = parse(Int, toks[2]); KMAX = parse(Int, toks[3])
            dims[b] = (IMAX, JMAX, KMAX)
        end

        # -------- payload: X, Y, Z for each block --------
        tmp = Float64[]
        blocks = Vector{Block}(undef, nblocks)
        for b in 1:nblocks
            IMAX, JMAX, KMAX = dims[b]
            n = IMAX*JMAX*KMAX

            # X
            _read_n_ascii_floats!(io, n, tmp)
            X = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
            idx = 1
            @inbounds for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX   # i–j–k order
                X[i,j,k] = tmp[idx]; idx += 1
            end

            # Y
            _read_n_ascii_floats!(io, n, tmp)
            Y = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
            idx = 1
            @inbounds for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                Y[i,j,k] = tmp[idx]; idx += 1
            end

            # Z
            _read_n_ascii_floats!(io, n, tmp)
            Z = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
            idx = 1
            @inbounds for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                Z[i,j,k] = tmp[idx]; idx += 1
            end

            blocks[b] = Block(X, Y, Z; index=b-1)  # store 0-based index for consistency
        end

        return blocks
    end
end

# Alias aligned with your earlier examples
read_blocks(path::AbstractString) = read_plot3D_ascii(path)

# ------------------------------------------------------------
# Binary helpers (match our writer: raw header + raw payload)
# ------------------------------------------------------------
_unsigned(::Type{Int32})   = UInt32
_unsigned(::Type{UInt32})  = UInt32
_unsigned(::Type{Float32}) = UInt32
_unsigned(::Type{Float64}) = UInt64

# read a single primitive with chosen endianness
function _read_num(io::IO, ::Type{T}, big_endian::Bool) where {T<:Union{Int32,UInt32,Float32,Float64}}
    nread = read(io, T)
    if big_endian == Base.isbigendian()
        return nread
    else
        u = reinterpret(_unsigned(T), nread)
        return reinterpret(T, bswap(u))
    end
end

# read N values of type T with chosen endianness
function _read_vector(io::IO, ::Type{T}, N::Int, big_endian::Bool) where {T<:Union{Float32,Float64,UInt32}}
    out = Vector{T}(undef, N)
    @inbounds for i in 1:N
        out[i] = _read_num(io, T, big_endian)
    end
    return out
end

# ------------------------------------------------------------
# Public: read BINARY Plot3D multi-block file
# ------------------------------------------------------------
"""
    read_plot3D_binary(path; double_precision=true, big_endian=false) -> Vector{Block}

Read a **binary** Plot3D file that matches `write_plot3D` in this library
(i.e., raw header followed by raw payload, **no** Fortran record markers).

Header layout:
- `UInt32(nblocks)`
- then for each block: `UInt32(IMAX)`, `UInt32(JMAX)`, `UInt32(KMAX)`

Payload:
- For each block, X then Y then Z arrays in **i–j–k** order (i fastest).
- `double_precision=true` expects `Float64`; set `false` for `Float32`.
- `big_endian` lets you request big-endian decoding (default assumes little-endian files on little-endian hosts).
"""
function read_plot3D_binary(path::AbstractString; double_precision::Bool=true, big_endian::Bool=false)
    open(path, "r") do io
        # header count
        nblocks_u = _read_num(io, UInt32, big_endian)
        nblocks = Int(nblocks_u)

        # per-block dims
        dims = Vector{NTuple{3,Int}}(undef, nblocks)
        @inbounds for b in 1:nblocks
            IMAX = Int(_read_num(io, UInt32, big_endian))
            JMAX = Int(_read_num(io, UInt32, big_endian))
            KMAX = Int(_read_num(io, UInt32, big_endian))
            dims[b] = (IMAX, JMAX, KMAX)
        end

        Tout = double_precision ? Float64 : Float32
        blocks = Vector{Block}(undef, nblocks)

        # payload: i–j–k order X, then Y, then Z
        @inbounds for b in 1:nblocks
            IMAX, JMAX, KMAX = dims[b]
            n = IMAX*JMAX*KMAX

            # X
            vx = _read_vector(io, Tout, n, big_endian)
            X = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
            idx = 1
            for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                X[i,j,k] = Float64(vx[idx]); idx += 1
            end

            # Y
            vy = _read_vector(io, Tout, n, big_endian)
            Y = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
            idx = 1
            for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                Y[i,j,k] = Float64(vy[idx]); idx += 1
            end

            # Z
            vz = _read_vector(io, Tout, n, big_endian)
            Z = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
            idx = 1
            for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                Z[i,j,k] = Float64(vz[idx]); idx += 1
            end

            blocks[b] = Block(X, Y, Z; index=b-1)
        end

        return blocks
    end
end

end # module

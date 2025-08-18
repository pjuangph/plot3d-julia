module ReadPlot3D

import ..Block3D: Block

export read_plot3D_ascii, read_blocks, read_plot3D_binary

const HOST_IS_LITTLE_ENDIAN = ENDIAN_BOM == UInt32(0x04030201)

# ---------------- ASCII reader (unchanged) ----------------
function _read_n_ascii_floats!(io::IO, N::Int, out::Vector{Float64})
    empty!(out); sizehint!(out, N)
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

"""
    read_plot3D_ascii(path::AbstractString) -> Vector{Block}
"""
function read_plot3D_ascii(path::AbstractString)
    open(path, "r") do io
        firstline = nothing
        for ln in eachline(io)
            if !isempty(strip(ln)); firstline = ln; break; end
        end
        firstline === nothing && error("Empty file: $path")
        nblocks = parse(Int, split(firstline)[1])

        dims = Vector{NTuple{3,Int}}(undef, nblocks)
        for b in 1:nblocks
            ln = ""
            while isempty(strip(ln)); ln = readline(io); end
            toks = split(ln); length(toks) >= 3 || error("Bad size line for block $b: '$ln'")
            IMAX = parse(Int, toks[1]); JMAX = parse(Int, toks[2]); KMAX = parse(Int, toks[3])
            dims[b] = (IMAX, JMAX, KMAX)
        end

        tmp = Float64[]
        blocks = Vector{Block}(undef, nblocks)
        for b in 1:nblocks
            IMAX, JMAX, KMAX = dims[b]; n = IMAX*JMAX*KMAX
            _read_n_ascii_floats!(io, n, tmp)
            X = reshape(copy(tmp), IMAX, JMAX, KMAX)
            _read_n_ascii_floats!(io, n, tmp)
            Y = reshape(copy(tmp), IMAX, JMAX, KMAX)
            _read_n_ascii_floats!(io, n, tmp)
            Z = reshape(copy(tmp), IMAX, JMAX, KMAX)
            blocks[b] = Block(X, Y, Z; index=b-1)
        end
        return blocks
    end
end

read_blocks(path::AbstractString) = read_plot3D_ascii(path)

# ---------------- Binary helpers ----------------
_unsigned(::Type{Int32})   = UInt32
_unsigned(::Type{UInt32})  = UInt32
_unsigned(::Type{Float32}) = UInt32
_unsigned(::Type{Float64}) = UInt64

function _read_num(io::IO, ::Type{T}, big_endian::Bool) where {T<:Union{Int32,UInt32,Float32,Float64}}
    x = read(io, T)
    

    if big_endian == !HOST_IS_LITTLE_ENDIAN
        return x
    else
        u = reinterpret(_unsigned(T), x)
        return reinterpret(T, bswap(u))
    end
end

function _read_vec(io::IO, ::Type{T}, N::Int, big_endian::Bool) where {T<:Union{Float32,Float64,UInt32}}
    out = Vector{T}(undef, N)
    @inbounds for i in 1:N
        out[i] = _read_num(io, T, big_endian)
    end
    out
end

# Fortran unformatted record: [len::UInt32] payload [len::UInt32]
function _read_record_bytes(io::IO, big_endian::Bool)
    len = _read_num(io, UInt32, big_endian)
    bytes = read(io, len)
    len2 = _read_num(io, UInt32, big_endian)
    len == len2 || error("Fortran record length mismatch ($len vs $len2).")
    bytes
end

# ---------------- Binary reader (Fortran or raw) ----------------
"""
    read_plot3D_binary(path; format=:fortran, double_precision=true, big_endian=false) -> Vector{Block}

Read a **binary** Plot3D file.

- `format = :fortran` (**default**): Fortran unformatted sequential with 4-byte record markers.
- `format = :raw`: raw header + payload (no record markers), matching the alternate writer path.

Header:
- nblocks (UInt32), then IMAX/JMAX/KMAX per block (UInt32 each).

Payload:
- For each block, X then Y then Z in **i–j–k** order.
"""
function read_plot3D_binary(path::AbstractString; format::Symbol=:fortran, double_precision::Bool=true, big_endian::Bool=false)
    open(path, "r") do io
        Tout = double_precision ? Float64 : Float32

        # ---- read header ----
        nblocks::Int = 0
        dims::Vector{NTuple{3,Int}} = NTuple{3,Int}[]

        if format === :raw
            nblocks = Int(_read_num(io, UInt32, big_endian))
            dims = Vector{NTuple{3,Int}}(undef, nblocks)
            @inbounds for b in 1:nblocks
                IMAX = Int(_read_num(io, UInt32, big_endian))
                JMAX = Int(_read_num(io, UInt32, big_endian))
                KMAX = Int(_read_num(io, UInt32, big_endian))
                dims[b] = (IMAX, JMAX, KMAX)
            end
        elseif format === :fortran
            # nblocks record
            nb_rec = _read_record_bytes(io, big_endian)
            nblocks = Int(reinterpret(UInt32, nb_rec)[1])
            @info "read plot3d binary reading $(nblocks) blocks"

            dims = Vector{NTuple{3,Int}}(undef, nblocks)
            @inbounds for b in 1:nblocks
                drec = _read_record_bytes(io, big_endian)
                u = reinterpret(UInt32, drec)
                length(u) >= 3 || error("Dimension record too short for block $b")
                dims[b] = (Int(u[1]), Int(u[2]), Int(u[3]))
            end
        else
            error("read_plot3D_binary: unknown format $(format). Use :fortran or :raw.")
        end

        # ---- payload ----
        blocks = Vector{Block}(undef, nblocks)

        if format === :raw
            @inbounds for b in 1:nblocks
                IMAX, JMAX, KMAX = dims[b]; n = IMAX*JMAX*KMAX

                vx = _read_vec(io, Tout, n, big_endian)
                vy = _read_vec(io, Tout, n, big_endian)
                vz = _read_vec(io, Tout, n, big_endian)

                X = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
                Y = similar(X); Z = similar(X)
                idx = 1
                for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                    X[i,j,k] = Float64(vx[idx])
                    Y[i,j,k] = Float64(vy[idx])
                    Z[i,j,k] = Float64(vz[idx])
                    idx += 1
                end
                blocks[b] = Block(X, Y, Z; index=b-1)
            end

        else # :fortran
            @inbounds for b in 1:nblocks
                IMAX, JMAX, KMAX = dims[b]; n = IMAX*JMAX*KMAX

                # X record
                xb = _read_record_bytes(io, big_endian)
                vx = reinterpret(Tout, xb)
                length(vx) == n || error("X record size mismatch on block $b")
                # Y record
                yb = _read_record_bytes(io, big_endian)
                vy = reinterpret(Tout, yb)
                length(vy) == n || error("Y record size mismatch on block $b")
                # Z record
                zb = _read_record_bytes(io, big_endian)
                vz = reinterpret(Tout, zb)
                length(vz) == n || error("Z record size mismatch on block $b")

                X = Array{Float64,3}(undef, IMAX, JMAX, KMAX)
                Y = similar(X); Z = similar(X)
                idx = 1
                for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                    X[i,j,k] = Float64(vx[idx])
                    Y[i,j,k] = Float64(vy[idx])
                    Z[i,j,k] = Float64(vz[idx])
                    idx += 1
                end
                blocks[b] = Block(X, Y, Z; index=b-1)
            end
        end

        return blocks
    end
end

end # module

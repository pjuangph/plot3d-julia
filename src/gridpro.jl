module GridPro

using ..Block3D: Block

export read_gridpro_to_blocks, read_gridpro_connectivity

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
_parse_header(st::AbstractString) = try
    parts = split(strip(st))
    length(parts) < 3 && return nothing
    i = Int(floor(parse(Float64, parts[1])))
    j = Int(floor(parse(Float64, parts[2])))
    k = Int(floor(parse(Float64, parts[3])))
    (i > 0 && j > 0 && k > 0) ? (i, j, k) : nothing
catch
    nothing
end

# ------------------------------------------------------------------
# read_gridpro_to_blocks
# ------------------------------------------------------------------
"""
    read_gridpro_to_blocks(filename; encoding="utf-8", comment_prefixes=("#","//"))

Reads a GridPro-style ASCII file that lists a header line `IMAX JMAX KMAX`, followed by
`IMAX*JMAX*KMAX` triplets `(x y z)` (whitespace-separated). Repeats for each block.

Returns `Vector{Block}`.
"""
function read_gridpro_to_blocks(filename::AbstractString;
                                encoding::AbstractString="utf-8",
                                comment_prefixes::Tuple{Vararg{AbstractString}}=("#","//"))
    blocks = Block[]
    open(filename, "r"; encoding=encoding) do io
        # seek first header
        line = readline(io; keep=true)
        while !eof(io) || !isempty(line)
            st = strip(line)
            if !isempty(st) && all(!startswith(st, p) for p in comment_prefixes)
                hdr = _parse_header(st)
                if hdr !== nothing
                    # main block loop from this header
                    while true
                        if hdr === nothing
                            break
                        end
                        IMAX, JMAX, KMAX = hdr
                        total_pts = IMAX*JMAX*KMAX
                        expected = 3*total_pts
                        # Slurp floats for this block
                        vals = Float64[]
                        sizehint!(vals, expected)
                        while length(vals) < expected && !eof(io)
                            s = strip(readline(io))
                            isempty(s) && continue
                            any(startswith(s, p) for p in comment_prefixes) && continue
                            # replace commas just in case
                            s = replace(s, ',' => ' ')
                            for tok in eachsplit(s)
                                isempty(tok) && continue
                                push!(vals, parse(Float64, tok))
                                length(vals) == expected && break
                            end
                        end
                        length(vals) == expected || error("Expected $expected floats for block ($IMAX,$JMAX,$KMAX), got $(length(vals)).")
                        # reshape
                        X = Array{Float64}(undef, IMAX,JMAX,KMAX)
                        Y = similar(X); Z = similar(X)
                        # vals is row of triplets
                        idx = 1
                        @inbounds for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
                            X[i,j,k] = vals[idx];   Y[i,j,k] = vals[idx+1]; Z[i,j,k] = vals[idx+2]
                            idx += 3
                        end
                        push!(blocks, Block(X,Y,Z))
                        # advance to next header (skip blanks/comments)
                        hdr = nothing
                        while !eof(io)
                            line2 = readline(io)
                            st2 = strip(line2)
                            isempty(st2) && continue
                            any(startswith(st2, p) for p in comment_prefixes) && continue
                            hdr = _parse_header(st2)
                            break
                        end
                        if hdr === nothing
                            return blocks
                        end
                    end
                end
            end
            eof(io) && break
            line = readline(io; keep=true)
        end
    end
    return blocks
end

# ------------------------------------------------------------------
# read_gridpro_connectivity
# ------------------------------------------------------------------
"""
    read_gridpro_connectivity(file_path; sb_zero_based_in_file=true, index_zero_based_in_file=true)

Parse a GridPro connectivity file (with lines starting `SB` and `P`) and return a Dict
mirroring the Python output:

Keys:
- `"face_matches"   => Vector{Dict}`
- `"outer_faces"    => Vector{Dict}`
- `"bc_group"       => Dict{String,Vector{Dict}}` with keys inlet/outlet/symm_slip/wall
- `"gif_faces"      => Vector{Dict}`
- `"periodic_faces" => Vector{Dict}`
- `"volume_zones"   => Vector{Dict}`
"""
function read_gridpro_connectivity(file_path::AbstractString;
                                   sb_zero_based_in_file::Bool=true,
                                   index_zero_based_in_file::Bool=true)
    superblock_ptys = Int[]
    patches = Vector{Dict{String,Any}}()

    sb_offset = sb_zero_based_in_file ? 0 : -1
    idx_offset = index_zero_based_in_file ? 0 : -1

    function parse_patch(tokens::Vector{SubString{String}})
        length(tokens) ≥ 21 || error("Patch line too short: $(join(tokens, ' '))")
        pid  = parse(Int, tokens[2])
        sb1_raw = parse(Int, tokens[3]); sf1 = parse(Int, tokens[4])
        sb2_raw = parse(Int, tokens[5]); sf2 = parse(Int, tokens[6])
        sb1 = sb1_raw - 1
        sb2 = (sb2_raw - 1 > -1) ? (sb2_raw - 1) : -1
        fmap = String(tokens[7])

        L1i = parse(Int, tokens[8 ]) + idx_offset
        L1j = parse(Int, tokens[9 ]) + idx_offset
        L1k = parse(Int, tokens[10]) + idx_offset
        H1i = parse(Int, tokens[11]) + idx_offset
        H1j = parse(Int, tokens[12]) + idx_offset
        H1k = parse(Int, tokens[13]) + idx_offset
        L2i = parse(Int, tokens[14]) + idx_offset
        L2j = parse(Int, tokens[15]) + idx_offset
        L2k = parse(Int, tokens[16]) + idx_offset
        H2i = parse(Int, tokens[17]) + idx_offset
        H2j = parse(Int, tokens[18]) + idx_offset
        H2k = parse(Int, tokens[19]) + idx_offset

        pty  = parse(Int, tokens[20])
        lbid = parse(Int, tokens[21])

        push!(patches, Dict(
            "pid"=>pid, "sb1"=>sb1, "sf1"=>sf1, "sb2"=>sb2, "sf2"=>sf2,
            "fmap"=>fmap,
            "L1i"=>L1i, "L1j"=>L1j, "L1k"=>L1k, "H1i"=>H1i, "H1j"=>H1j, "H1k"=>H1k,
            "L2i"=>L2i, "L2j"=>L2j, "L2k"=>L2k, "H2i"=>H2i, "H2j"=>H2j, "H2k"=>H2k,
            "pty"=>pty, "lbid"=>lbid
        ))
        return nothing
    end

    open(file_path, "r"; encoding="utf-8") do io
        for raw in eachline(io)
            line = strip(raw)
            isempty(line) && continue
            (startswith(line, "#") || startswith(line, "//")) && continue
            toks = split(line)
            tag = toks[1]
            if tag == "SB"
                # Last-but-one token is the pty in the provided sources
                push!(superblock_ptys, parse(Int, toks[end-1]))
            elseif tag == "P"
                parse_patch(toks)
            end
        end
    end

    face_dict(sb, imin,jmin,kmin, imax,jmax,kmax; pty::Int=-1) = Dict(
        "block_index"=>sb,
        "IMIN"=>imin, "JMIN"=>jmin, "KMIN"=>kmin,
        "IMAX"=>imax, "JMAX"=>jmax, "KMAX"=>kmax,
        "id"=>pty
    )
    function pair_dict(sb1, r1::NTuple{6,Int}, sb2, r2::NTuple{6,Int})
        L1i,L1j,L1k,H1i,H1j,H1k = r1
        L2i,L2j,L2k,H2i,H2j,H2k = r2
        return Dict(
            "block1"=>face_dict(sb1, L1i,L1j,L1k, H1i,H1j,H1k),
            "block2"=>face_dict(sb2, L2i,L2j,L2k, H2i,H2j,H2k),
        )
    end

    # connections: pty in {1,3} and sb2 != -1
    connections = Dict{String,Any}[]
    for p in patches
        if (p["sb2"] != -1) && (p["pty"] in (1,3))
            push!(connections, pair_dict(
                p["sb1"], (p["L1i"],p["L1j"],p["L1k"],p["H1i"],p["H1j"],p["H1k"]),
                p["sb2"], (p["L2i"],p["L2j"],p["L2k"],p["H2i"],p["H2j"],p["H2k"]),
            ))
        end
    end

    # outer faces: pty not in [6,5,4,2] and sb2 == -1
    outer_faces = Dict{String,Int}[]
    pty_exclude = Set([6,5,4,2])
    for p in patches
        if p["sb2"] == -1 && !(p["pty"] in pty_exclude)
            push!(outer_faces, face_dict(
                p["sb1"], p["L1i"],p["L1j"],p["L1k"], p["H1i"],p["H1j"],p["H1k"]; pty=p["pty"])
            )
        end
    end

    # bc groups
    function bc_faces_for(pty_val::Int)
        v = Dict{String,Int}[]
        for p in patches
            if p["pty"] == pty_val
                push!(v, face_dict(p["sb1"], p["L1i"],p["L1j"],p["L1k"], p["H1i"],p["H1j"],p["H1k"]; pty=p["pty"]))
            end
        end
        return v
    end
    bc_group = Dict(
        "inlet"     => bc_faces_for(5),
        "outlet"    => bc_faces_for(6),
        "symm_slip" => bc_faces_for(4),
        "wall"      => bc_faces_for(2),
    )

    # periodic faces (pty == 3)
    periodic_faces = Dict{String,Any}[]
    for p in patches
        if p["pty"] == 3 && p["sb2"] != -1
            push!(periodic_faces, pair_dict(
                p["sb1"], (p["L1i"],p["L1j"],p["L1k"],p["H1i"],p["H1j"],p["H1k"]),
                p["sb2"], (p["L2i"],p["L2j"],p["L2k"],p["H2i"],p["H2j"],p["H2k"]),
            ))
        end
    end

    # GIF faces (pty in 12..21 or == 1000)
    gif_faces = Dict{String,Int}[]
    for p in patches
        pty = p["pty"]
        if (12 ≤ pty ≤ 21) || (pty == 1000)
            d = face_dict(p["sb1"], p["L1i"],p["L1j"],p["L1k"], p["H1i"],p["H1j"],p["H1k"]; pty=pty)
            d["id"] = pty
            push!(gif_faces, d)
        end
    end

    # Volume zones: from superblock_ptys (odd→fluid, even→solid), contiguous groups
    volume_zones = Vector{Dict{String,Any}}()
    prev = ""
    cid = 0
    for (id,v) in enumerate(superblock_ptys)
        zone = (v % 2 != 0) ? "fluid" : "solid"
        if zone != prev
            cid += 1
            prev = zone
        end
        push!(volume_zones, Dict("block_index"=>id-1, "zone_type"=>zone, "contiguous_index"=>cid))
    end

    return Dict(
        "face_matches"=>connections,
        "outer_faces"=>outer_faces,
        "bc_group"=>bc_group,
        "gif_faces"=>gif_faces,
        "periodic_faces"=>periodic_faces,
        "volume_zones"=>volume_zones,
    )
end

end # module

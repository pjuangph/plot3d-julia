# connectivity.jl — lives directly in Plot3D (no module block)

import LinearAlgebra: norm
import Statistics: mean

using .Block3D: Block
using .Face3D: Face, set_block_index
# Pull helpers that live directly in the Plot3D module (from face algorithms)
using .Plot3D: faces_match, create_face_from_diagonals, face_matches_to_dict

# -----------------------------------------------------------------------------
# Types
# -----------------------------------------------------------------------------
"""
Light container for a set of face matches.

Use `pairs` when you only need (block, face-name) tuples,
or use the dicts returned by `find_matching_blocks` / `connectivity_fast`
if you need index ranges.
"""
struct FaceMatchSet
    pairs::Vector{Tuple{Tuple{Int,Symbol},Tuple{Int,Symbol}}}
end

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
"""
    point_match(p::NTuple{3,Real}, q::NTuple{3,Real}; tol=1e-8)

Approximate equality for 3D points.
"""
point_match(p::NTuple{3,Real}, q::NTuple{3,Real}; tol::Real=1e-8) =
    isapprox(p[1], q[1]; atol=tol) &&
    isapprox(p[2], q[2]; atol=tol) &&
    isapprox(p[3], q[3]; atol=tol)

"""
    select_multi_dimensional(A, ir::UnitRange, jr::UnitRange, kr::UnitRange)

Return a `@view` into a 3D array.
"""
select_multi_dimensional(A::AbstractArray, ir::UnitRange, jr::UnitRange, kr::UnitRange) = @view A[ir, jr, kr]

# Small helper: centroid of a block. Uses b.cx/b.cy/b.cz if present; otherwise computes mean.
function _block_centroid(b::Block)
    try
        return (getfield(b, :cx), getfield(b, :cy), getfield(b, :cz))
    catch
        return (mean(vec(b.X)), mean(vec(b.Y)), mean(vec(b.Z)))
    end
end

# Internal: get canonical 2D (X,Y,Z) face arrays by name
# face_name: "imin"|"imax"|"jmin"|"jmax"|"kmin"|"kmax"
function _get_face_arrays_named(b::Block, face_name::String)
    if face_name === "imin"
        return (view(b.X, 1, :, :),           view(b.Y, 1, :, :),           view(b.Z, 1, :, :))
    elseif face_name === "imax"
        return (view(b.X, b.IMAX, :, :),      view(b.Y, b.IMAX, :, :),      view(b.Z, b.IMAX, :, :))
    elseif face_name === "jmin"
        return (view(b.X, :, 1, :),           view(b.Y, :, 1, :),           view(b.Z, :, 1, :))
    elseif face_name === "jmax"
        return (view(b.X, :, b.JMAX, :),      view(b.Y, :, b.JMAX, :),      view(b.Z, :, b.JMAX, :))
    elseif face_name === "kmin"
        return (view(b.X, :, :, 1),           view(b.Y, :, :, 1),           view(b.Z, :, :, 1))
    elseif face_name === "kmax"
        return (view(b.X, :, :, b.KMAX),      view(b.Y, :, :, b.KMAX),      view(b.Z, :, :, b.KMAX))
    else
        error("Unknown face name: $face_name")
    end
end

# Provide all canonical names a block actually has (skip k-faces if KMAX==1)
function _all_face_names(b::Block)
    names = String["imin","imax","jmin","jmax"]
    if b.KMAX > 1
        push!(names, "kmin"); push!(names, "kmax")
    end
    return names
end

# Build a Face object that spans the full canonical face (index-space) for a block
function _make_full_face(b::Block, face_name::String)
    if face_name === "imin"
        return create_face_from_diagonals(b, 0, 0, 0, 0, b.JMAX-1, b.KMAX-1)
    elseif face_name === "imax"
        return create_face_from_diagonals(b, b.IMAX-1, 0, 0, b.IMAX-1, b.JMAX-1, b.KMAX-1)
    elseif face_name === "jmin"
        return create_face_from_diagonals(b, 0, 0, 0, b.IMAX-1, 0, b.KMAX-1)
    elseif face_name === "jmax"
        return create_face_from_diagonals(b, 0, b.JMAX-1, 0, b.IMAX-1, b.JMAX-1, b.KMAX-1)
    elseif face_name === "kmin"
        return create_face_from_diagonals(b, 0, 0, 0, b.IMAX-1, b.JMAX-1, 0)
    else # "kmax"
        return create_face_from_diagonals(b, 0, 0, b.KMAX-1, b.IMAX-1, b.JMAX-1, b.KMAX-1)
    end
end

# -----------------------------------------------------------------------------
# find_matching_blocks
# -----------------------------------------------------------------------------
"""
    find_matching_blocks(blocks; tol=1e-8)

Find matching faces across all blocks by comparing canonical faces
(`imin`, `imax`, `jmin`, `jmax`, `kmin`, `kmax` when present).

Returns a vector of **match dictionaries** compatible with your Python output,
using the same schema as produced by `face_matches_to_dict`.

Each element looks like:
    Dict(
    "block1" => Dict("block_index"=>i, "IMIN"=>..., "JMIN"=>..., "KMIN"=>..., "IMAX"=>..., "JMAX"=>..., "KMAX"=>..., "id"=>...),
    "block2" => Dict( … same keys … )
    )

    The index ranges are filled so that lower/upper corners correspond between faces.
"""
function find_matching_blocks(blocks::Vector{Block}; tol::Real=1e-8)
    n = length(blocks)
    matches = Dict{String,Any}[]

    # Precompute centroids for cheap pruning
    cents = [_block_centroid(b) for b in blocks]

    # Rough global scale for centroid threshold
    xs = vcat([vec(b.X) for b in blocks]...)
    ys = vcat([vec(b.Y) for b in blocks]...)
    zs = vcat([vec(b.Z) for b in blocks]...)
    bbox = (maximum(xs)-minimum(xs)) + (maximum(ys)-minimum(ys)) + (maximum(zs)-minimum(zs))
    centroid_thresh = max(1e-12, 1e-6 * bbox)

    face_names = [ _all_face_names(b) for b in blocks ]

    for i in 1:n
        bi = blocks[i]
        for j in i+1:n
            bj = blocks[j]

            # centroid gating
            dx = cents[i][1]-cents[j][1]; dy = cents[i][2]-cents[j][2]; dz = cents[i][3]-cents[j][3]
            dcent = sqrt(dx*dx + dy*dy + dz*dz)
            dcent <= centroid_thresh || continue

            # try all face name pairs
            for fi in face_names[i]
                Ai = _get_face_arrays_named(bi, fi)
                for fj in face_names[j]
                    Aj = _get_face_arrays_named(bj, fj)
                    ok, _ = faces_match(Ai, Aj; tol=tol)
                    if ok
                        f1 = _make_full_face(bi, fi)
                        f2 = _make_full_face(bj, fj)
                        set_block_index(f1, i-1); set_block_index(f2, j-1)
                        push!(matches, face_matches_to_dict(f1, f2, bi, bj))
                    end
                end
            end
        end
    end

    return matches
end

# -----------------------------------------------------------------------------
# combinations_of_nearest_blocks
# -----------------------------------------------------------------------------
"""
    combinations_of_nearest_blocks(blocks; nearest_nblocks=4)

Heuristic block pairing by nearest **block centroids**.
Returns a vector of unique `(i,j)` index pairs (1-based), with i < j.
"""
function combinations_of_nearest_blocks(blocks::Vector{Block}; nearest_nblocks::Int=4)
    n = length(blocks)
    if n ≤ 1
        return Tuple{Int,Int}[]
    end

    cents = [_block_centroid(b) for b in blocks]
    pairs = Set{Tuple{Int,Int}}()

    for i in 1:n
        dists = Vector{Tuple{Int,Float64}}()
        ci = cents[i]
        for j in 1:n
            j == i && continue
            cj = cents[j]
            dx = ci[1]-cj[1]; dy = ci[2]-cj[2]; dz = ci[3]-cj[3]
            push!(dists, (j, sqrt(dx*dx + dy*dy + dz*dz)))
        end
        sort!(dists, by = x->x[2])
        m = min(nearest_nblocks, length(dists))
        for k in 1:m
            j = dists[k][1]
            i < j ? push!(pairs, (i,j)) : push!(pairs, (j,i))
        end
    end

    return collect(pairs)
end

# -----------------------------------------------------------------------------
# get_face_intersection (simple full-face variant)
# -----------------------------------------------------------------------------
"""
    get_face_intersection(face1::Face, face2::Face, block1::Block, block2::Block; tol=1e-8)

Return a vector with **one** match dictionary describing how `face1` maps to `face2`
(using the same schema as `face_matches_to_dict`). This **first pass** assumes the
faces fully match (typical CFD abutting faces). If you need **partial overlaps**,
extend this to compute tight index windows by projecting and scanning parameter lines.

Returns: `Vector{Dict{String,Any}}` with length 1 when matched, or `Vector{Dict{String,Any}}()` if not matched.
"""
function get_face_intersection(face1::Face, face2::Face, block1::Block, block2::Block; tol::Real=1e-8)
    # Build 2D arrays for each face to test a match (fast corner test with reversals)
    function _arrays_from_face(b::Block, f::Face)
        if f.IMIN == f.IMAX
            a = f.IMIN + 1
            return (view(b.X, a, f.JMIN+1:f.JMAX+1, f.KMIN+1:f.KMAX+1),
                    view(b.Y, a, f.JMIN+1:f.JMAX+1, f.KMIN+1:f.KMAX+1),
                    view(b.Z, a, f.JMIN+1:f.JMAX+1, f.KMIN+1:f.KMAX+1))
        elseif f.JMIN == f.JMAX
            b0 = f.JMIN + 1
            return (view(b.X, f.IMIN+1:f.IMAX+1, b0, f.KMIN+1:f.KMAX+1),
                    view(b.Y, f.IMIN+1:f.IMAX+1, b0, f.KMIN+1:f.KMAX+1),
                    view(b.Z, f.IMIN+1:f.IMAX+1, b0, f.KMIN+1:f.KMAX+1))
        else
            c = f.KMIN + 1
            return (view(b.X, f.IMIN+1:f.IMAX+1, f.JMIN+1:f.JMAX+1, c),
                    view(b.Y, f.IMIN+1:f.IMAX+1, f.JMIN+1:f.JMAX+1, c),
                    view(b.Z, f.IMIN+1:f.IMAX+1, f.JMIN+1:f.JMAX+1, c))
        end
    end

    A1 = _arrays_from_face(block1, face1)
    A2 = _arrays_from_face(block2, face2)
    ok, _ = faces_match(A1, A2; tol=tol)
    if !ok
        return Dict{String,Any}[]  # not intersecting (within tolerance)
    end

    # For full-face matches, a single dict maps the lower/upper corners.
    return [face_matches_to_dict(face1, face2, block1, block2)]
end

# -----------------------------------------------------------------------------
# connectivity_fast (Python-compatible surface)
# -----------------------------------------------------------------------------
"""
    connectivity_fast(blocks; tol=1e-8)

Lightweight, pragmatic connectivity finder.

Returns:
- `face_matches::Vector{Dict{String,Any}}` (same schema as Python’s `connectivity_fast`)
- `outer_faces::Vector{Dict{String,Int}}`   (canonical faces that did not match)

Notes:
- This variant compares **canonical faces** only (imin/imax/jmin/jmax/kmin/kmax) and
  does **not** split partial overlaps. If you need splitting, we can extend this.
"""
function connectivity_fast(blocks::Vector{Block}; tol::Real=1e-8)
    n = length(blocks)
    face_names = [ _all_face_names(b) for b in blocks ]

    # Build (block, face_name, Face) inventory
    faces = Vector{Tuple{Int,String,Face}}()
    for i in 1:n
        bi = blocks[i]
        for fn in face_names[i]
            f = _make_full_face(bi, fn)
            set_block_index(f, i-1)  # 0-based in dicts, like Python
            push!(faces, (i, fn, f))
        end
    end

    # Try to match every pair across different blocks
    matched_pairs = Vector{Tuple{Int,Int}}()  # indices into `faces`
    match_dicts   = Dict{String,Any}[]

    for a in 1:length(faces)-1
        ia, fna, fa = faces[a]
        Aa = _get_face_arrays_named(blocks[ia], fna)
        for b in a+1:length(faces)
            ib, fnb, fb = faces[b]
            ia == ib && continue # skip same-block duplicates
            Ab = _get_face_arrays_named(blocks[ib], fnb)
            ok, _ = faces_match(Aa, Ab; tol=tol)
            if ok
                push!(matched_pairs, (a,b))
                push!(match_dicts, face_matches_to_dict(fa, fb, blocks[ia], blocks[ib]))
            end
        end
    end

    # Any face not participating in a match is considered an "outer face"
    used = Set{Int}()
    for (a,b) in matched_pairs
        push!(used, a); push!(used, b)
    end

    outer_faces = Dict{String,Int}[]
    for idx in 1:length(faces)
        if !(idx in used)
            _, _, f = faces[idx]
            push!(outer_faces, Dict(
                "block_index" => f.BlockIndex,
                "IMIN" => f.IMIN, "JMIN" => f.JMIN, "KMIN" => f.KMIN,
                "IMAX" => f.IMAX, "JMAX" => f.JMAX, "KMAX" => f.KMAX,
                "id"   => f.id,
            ))
        end
    end

    return match_dicts, outer_faces
end

# -----------------------------------------------------------------------------
# block_connection_matrix
# -----------------------------------------------------------------------------
"""
    block_connection_matrix(blocks, all_faces_dicts) -> Matrix{Int8}

Return an `n×n` symmetric adjacency matrix where `C[i,j]=1` if any face between
block `i-1` and `j-1` was found to match (via `find_matching_blocks` here),
else 0. Diagonal is set to 1.

`all_faces_dicts` is accepted to mirror the Python signature, but this pragmatic
implementation recomputes matches from `blocks` directly (fine for now).
"""
function block_connection_matrix(blocks::Vector{Block}, all_faces_dicts)::Matrix{Int8}
    n = length(blocks)
    C = fill(Int8(0), n, n)
    for i in 1:n
        C[i,i] = 1
    end

    matches = find_matching_blocks(blocks)
    for m in matches
        b1 = m["block1"]; b2 = m["block2"]
        i = Int(b1["block_index"]) + 1
        j = Int(b2["block_index"]) + 1
        C[i,j] = 1
        C[j,i] = 1
    end
    return C
end

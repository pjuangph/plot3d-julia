# block_face_functions.jl — merged Block+Face helpers (no circular imports)

import LinearAlgebra: norm, dot
import Statistics: mean

# Core types we depend on
import .Block3D: Block
import .Face3D: Face, add_vertex, vertices_equals, index_equals, normal, match_indices, to_dict, set_block_index, set_face_id, is_edge
import .Utils: unique_pairs

# Face + block helpers (from block_face_functions.jl)
export reduce_blocks,
       get_faces,
       faces_match,
       find_matching_faces,
       get_outer_faces,          # both single-block and vector-of-blocks methods
       get_outer_face_dicts,     # both single-block and vector-of-blocks methods
       create_face_from_diagonals,
       find_connected_faces,
       find_closest_block,
       find_bounding_faces,
       split_face,
       find_face_nearest_point,
       outer_face_dict_to_list,
       match_faces_dict_to_list,
       face_matches_to_dict,
       touches_by_nodes,
       shared_point_fraction
# =============================================================================
# Block reduction (moved from BlockFunctions)
# =============================================================================
"""
    reduce_blocks(blocks::Vector{Block}, stride::Integer)

Downsample each block by taking every `stride` index in i, j, k (keeps endpoints).
Returns a new vector of `Block`s with reduced resolution. `stride` must be ≥1.
"""
function reduce_blocks(blocks::Vector{Block}, stride::Integer)
    stride ≥ 1 || throw(ArgumentError("stride must be ≥ 1 (got $stride)"))
    out = Block[]
    @inbounds for b in blocks
        # make index ranges that always include the last index
        is = unique!(sort!(vcat(1:stride:b.IMAX, b.IMAX)))
        js = unique!(sort!(vcat(1:stride:b.JMAX, b.JMAX)))
        ks = unique!(sort!(vcat(1:stride:max(b.KMAX,1), max(b.KMAX,1))))
        Xr = b.X[is, js, ks]
        Yr = b.Y[is, js, ks]
        Zr = b.Z[is, js, ks]
        push!(out, Block(Xr, Yr, Zr))
    end
    return out
end

# =============================================================================
# Face functions (your current facefunctions.jl content, unchanged)
# =============================================================================

# -----------------------------------------------------------------------------
# get_faces: dictionary of face slices by name
# -----------------------------------------------------------------------------
function get_faces(block::Block)
    faces = Dict{String,Tuple{AbstractArray,AbstractArray,AbstractArray}}()
    faces["imin"] = (view(block.X, 1, :, :),  view(block.Y, 1, :, :),  view(block.Z, 1, :, :))
    faces["imax"] = (view(block.X, block.IMAX, :, :), view(block.Y, block.IMAX, :, :), view(block.Z, block.IMAX, :, :))
    faces["jmin"] = (view(block.X, :, 1, :),  view(block.Y, :, 1, :),  view(block.Z, :, 1, :))
    faces["jmax"] = (view(block.X, :, block.JMAX, :), view(block.Y, :, block.JMAX, :), view(block.Z, :, block.JMAX, :))
    if block.KMAX > 1
        faces["kmin"] = (view(block.X, :, :, 1),  view(block.Y, :, :, 1),  view(block.Z, :, :, 1))
        faces["kmax"] = (view(block.X, :, :, block.KMAX), view(block.Y, :, :, block.KMAX), view(block.Z, :, :, block.KMAX))
    end
    return faces
end

# -----------------------------------------------------------------------------
# faces_match: compare corners only (fast, no copies)
# -----------------------------------------------------------------------------
function faces_match(
    face1::Tuple{AbstractArray,AbstractArray,AbstractArray},
    face2::Tuple{AbstractArray,AbstractArray,AbstractArray};
    tol::Real = 1e-12,
)
    X1, Y1, Z1 = face1
    X2, Y2, Z2 = face2
    size(X1) == size(X2) || return (false, nothing)

    nr, nc = size(X1)
    tol2 = tol*tol

    @inline function corner_ok(flip_ud::Bool, flip_lr::Bool)
        r1 = flip_ud ? nr : 1
        rN = flip_ud ? 1  : nr
        c1 = flip_lr ? nc : 1
        cN = flip_lr ? 1  : nc

        dx = X1[1,1]   - X2[r1,c1]; dy = Y1[1,1]   - Y2[r1,c1]; dz = Z1[1,1]   - Z2[r1,c1]
        (dx*dx + dy*dy + dz*dz <= tol2) || return false
        dx = X1[1,end] - X2[r1,cN]; dy = Y1[1,end] - Y2[r1,cN]; dz = Z1[1,end] - Z2[r1,cN]
        (dx*dx + dy*dy + dz*dz <= tol2) || return false
        dx = X1[end,1] - X2[rN,c1]; dy = Y1[end,1] - Y2[rN,c1]; dz = Z1[end,1] - Z2[rN,c1]
        (dx*dx + dy*dy + dz*dz <= tol2) || return false
        dx = X1[end,end] - X2[rN,cN]; dy = Y1[end,end] - Y2[rN,cN]; dz = Z1[end,end] - Z2[rN,cN]
        return dx*dx + dy*dy + dz*dz <= tol2
    end

    corner_ok(false,false) && return (true, (false,false))
    corner_ok(true, false) && return (true, (true, false))
    corner_ok(false,true ) && return (true, (false,true ))
    corner_ok(true, true ) && return (true, (true, true ))
    return (false, nothing)
end

# -----------------------------------------------------------------------------
# find_matching_faces
# -----------------------------------------------------------------------------
function find_matching_faces(block1::Block, block2::Block; tol::Real = 1e-8)
    faces1 = get_faces(block1)
    faces2 = get_faces(block2)
    for (name1, data1) in faces1
        for (name2, data2) in faces2
            match, flips = faces_match(data1, data2; tol=tol)
            match && return (name1, name2, flips)
        end
    end
    return (nothing, nothing, nothing)
end

# -----------------------------------------------------------------------------
# get_outer_faces for multiple blocks (Python-style convenience)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# get_outer_faces — single block
# -----------------------------------------------------------------------------
"""
    get_outer_faces(b::Block; block_index::Int=0)
        -> (non_matching::Vector{Face}, matching::Vector{Tuple{Face,Face}})

Return the outer faces of a block along with any self-matching face pairs
(faces that coincide within the same block). Mirrors Python's behaviour.
"""
function get_outer_faces(b::Block; block_index::Int=0)
    faces = Face[]
    push!(faces, create_face_from_diagonals(b, 0, 0, 0, 0, b.JMAX-1, b.KMAX-1))                 # imin
    push!(faces, create_face_from_diagonals(b, b.IMAX-1, 0, 0, b.IMAX-1, b.JMAX-1, b.KMAX-1))   # imax
    push!(faces, create_face_from_diagonals(b, 0, 0, 0, b.IMAX-1, 0, b.KMAX-1))                 # jmin
    push!(faces, create_face_from_diagonals(b, 0, b.JMAX-1, 0, b.IMAX-1, b.JMAX-1, b.KMAX-1))   # jmax
    if b.KMAX > 1
        push!(faces, create_face_from_diagonals(b, 0, 0, 0, b.IMAX-1, b.JMAX-1, 0))             # kmin
        push!(faces, create_face_from_diagonals(b, 0, 0, b.KMAX-1, b.IMAX-1, b.JMAX-1, b.KMAX-1)) # kmax
    end
    for f in faces
        set_block_index(f, block_index)
    end

    matching_pairs = Tuple{Int,Int}[]
    non_matching = Face[]
    for i in eachindex(faces)
        match_found = false
        for j in eachindex(faces)
            i == j && continue
            if vertices_equals(faces[i], faces[j])
                push!(matching_pairs, (i, j))
                match_found = true
            end
        end
        if !match_found
            push!(non_matching, faces[i])
        end
    end

    uniq_pairs = unique_pairs(matching_pairs)
    matching = [(faces[i], faces[j]) for (i, j) in uniq_pairs]
    return non_matching, matching
end

"""
    get_outer_faces(blocks::Vector{Block})
        -> (non_matching::Vector{Face}, matching::Vector{Tuple{Face,Face}})

Return outer faces and self-matched face pairs for every block.
"""
function get_outer_faces(blocks::Vector{Block})
    non_matching_all = Face[]
    matching_all = Tuple{Face,Face}[]
    for (bi, b) in enumerate(blocks)
        outer, matching = get_outer_faces(b; block_index=bi-1)
        append!(non_matching_all, outer)
        append!(matching_all, matching)
    end
    return non_matching_all, matching_all
end

# -----------------------------------------------------------------------------
# Dict helpers (if you need Python-style dict outputs)
# -----------------------------------------------------------------------------
"""
    get_outer_face_dicts(b::Block; block_index::Int=0) -> Vector{Dict{String,Int}}
"""
function get_outer_face_dicts(b::Block; block_index::Int=0)
    outer, _ = get_outer_faces(b; block_index=block_index)
    return [Dict(
        "block_index"=>f.BlockIndex,
        "IMIN"=>f.IMIN, "JMIN"=>f.JMIN, "KMIN"=>f.KMIN,
        "IMAX"=>f.IMAX, "JMAX"=>f.JMAX, "KMAX"=>f.KMAX,
        "id"=>f.id,
    ) for f in outer]
end

"""
    get_outer_face_dicts(blocks::Vector{Block}) -> Vector{Dict{String,Int}}
"""
function get_outer_face_dicts(blocks::Vector{Block})
    outer, _ = get_outer_faces(blocks)
    return [Dict(
        "block_index"=>f.BlockIndex,
        "IMIN"=>f.IMIN, "JMIN"=>f.JMIN, "KMIN"=>f.KMIN,
        "IMAX"=>f.IMAX, "JMAX"=>f.JMAX, "KMAX"=>f.KMAX,
        "id"=>f.id,
    ) for f in outer]
end
# -----------------------------------------------------------------------------
# create_face_from_diagonals
# -----------------------------------------------------------------------------
function create_face_from_diagonals(block::Block,
    imin::Int, jmin::Int, kmin::Int, imax::Int, jmax::Int, kmax::Int)

    newFace = Face(4)
    @inbounds begin
        if imin == imax
            i = imin
            for j in (jmin, jmax), k in (kmin, kmax)
                x = block.X[i+1, j+1, k+1]; y = block.Y[i+1, j+1, k+1]; z = block.Z[i+1, j+1, k+1]
                add_vertex(newFace, x, y, z, i, j, k)
            end
        elseif jmin == jmax
            j = jmin
            for i in (imin, imax), k in (kmin, kmax)
                x = block.X[i+1, j+1, k+1]; y = block.Y[i+1, j+1, k+1]; z = block.Z[i+1, j+1, k+1]
                add_vertex(newFace, x, y, z, i, j, k)
            end
        elseif kmin == kmax
            k = kmin
            for i in (imin, imax), j in (jmin, jmax)
                x = block.X[i+1, j+1, k+1]; y = block.Y[i+1, j+1, k+1]; z = block.Z[i+1, j+1, k+1]
                add_vertex(newFace, x, y, z, i, j, k)
            end
        end
    end
    return newFace
end

# -----------------------------------------------------------------------------
# Face index/geometry helpers (parity with Python Face methods)
# -----------------------------------------------------------------------------
@inline function _face_const_type(f::Face)
    if f.IMIN == f.IMAX
        return 0  # I-constant
    elseif f.JMIN == f.JMAX
        return 1  # J-constant
    elseif f.KMIN == f.KMAX
        return 2  # K-constant
    else
        return -1
    end
end

@inline function _face_index_ranges(f::Face)
    return ((f.IMIN, f.IMAX), (f.JMIN, f.JMAX), (f.KMIN, f.KMAX))
end

function _face_axis_extreme(f::Face, axis::AbstractString)
    axis_l = lowercase(axis)
    if isempty(f.vertices)
        return (0.0, 0.0)
    end
    xs = getfield.(f.vertices, Val(1))
    ys = getfield.(f.vertices, Val(2))
    zs = getfield.(f.vertices, Val(3))
    if axis_l == "x"
        return (minimum(xs), maximum(xs))
    elseif axis_l == "y"
        return (minimum(ys), maximum(ys))
    elseif axis_l == "z"
        return (minimum(zs), maximum(zs))
    else
        throw(ArgumentError("axis must be \"x\", \"y\", or \"z\" (got $axis)"))
    end
end

function _global_axis_extreme(blocks::AbstractVector{<:Block}, axis::AbstractString)
    axis_l = lowercase(axis)
    mins = Float64[]; maxs = Float64[]
    if axis_l == "x"
        for b in blocks
            push!(mins, minimum(b.X)); push!(maxs, maximum(b.X))
        end
    elseif axis_l == "y"
        for b in blocks
            push!(mins, minimum(b.Y)); push!(maxs, maximum(b.Y))
        end
    elseif axis_l == "z"
        for b in blocks
            push!(mins, minimum(b.Z)); push!(maxs, maximum(b.Z))
        end
    else
        throw(ArgumentError("axis must be \"x\", \"y\", or \"z\" (got $axis)"))
    end
    return (minimum(mins), maximum(maxs))
end

function _face_grid_points(block::Block, f::Face; stride_u::Integer=1, stride_v::Integer=1)
    const_type = _face_const_type(f)
    if const_type == -1
        # fall back to stored vertices (may be partial faces)
        return [(v[1], v[2], v[3]) for v in f.vertices]
    end

    su = max(1, Int(stride_u))
    sv = max(1, Int(stride_v))

    points = NTuple{3,Float64}[]
    (ir, jr, kr) = _face_index_ranges(f)

    if const_type == 0
        i = ir[1] + 1
        for j in jr[1]:su:jr[2], k in kr[1]:sv:kr[2]
            jj = j + 1; kk = k + 1
            push!(points, (block.X[i, jj, kk], block.Y[i, jj, kk], block.Z[i, jj, kk]))
        end
    elseif const_type == 1
        j = jr[1] + 1
        for i in ir[1]:su:ir[2], k in kr[1]:sv:kr[2]
            ii = i + 1; kk = k + 1
            push!(points, (block.X[ii, j, kk], block.Y[ii, j, kk], block.Z[ii, j, kk]))
        end
    else
        k = kr[1] + 1
        for i in ir[1]:su:ir[2], j in jr[1]:sv:jr[2]
            ii = i + 1; jj = j + 1
            push!(points, (block.X[ii, jj, k], block.Y[ii, jj, k], block.Z[ii, jj, k]))
        end
    end
    return points
end

@inline function _quantize_points(points::AbstractVector{<:NTuple{3,Float64}}, tol::Real)
    s = tol > 0 ? float(tol) : 1e-12
    return Set{NTuple{3,Int}}( (round(Int, p[1] / s), round(Int, p[2] / s), round(Int, p[3] / s)) for p in points )
end

function shared_point_fraction(
    face1::Face, face2::Face,
    block1::Block, block2::Block;
    tol_xyz::Real = 1e-8,
    stride_u::Integer = 1,
    stride_v::Integer = 1,
)
    pts1 = _face_grid_points(block1, face1; stride_u=stride_u, stride_v=stride_v)
    pts2 = _face_grid_points(block2, face2; stride_u=stride_u, stride_v=stride_v)
    if isempty(pts1) || isempty(pts2)
        return 0.0
    end

    Q1 = _quantize_points(pts1, tol_xyz)
    Q2 = _quantize_points(pts2, tol_xyz)
    if isempty(Q1) || isempty(Q2)
        return 0.0
    end

    shared = 0
    if length(Q1) ≤ length(Q2)
        for p in Q1
            shared += Int(p in Q2)
        end
    else
        for p in Q2
            shared += Int(p in Q1)
        end
    end
    denom = min(length(Q1), length(Q2))
    return denom > 0 ? shared / denom : 0.0
end

function touches_by_nodes(
    face1::Face, face2::Face,
    block1::Block, block2::Block;
    tol_xyz::Real = 1e-8,
    min_shared_frac::Real = 0.02,
    min_shared_abs::Integer = 4,
    stride_u::Integer = 1,
    stride_v::Integer = 1,
)
    pts1 = _face_grid_points(block1, face1; stride_u=stride_u, stride_v=stride_v)
    pts2 = _face_grid_points(block2, face2; stride_u=stride_u, stride_v=stride_v)
    if isempty(pts1) || isempty(pts2)
        return false
    end

    Q1 = _quantize_points(pts1, tol_xyz)
    Q2 = _quantize_points(pts2, tol_xyz)
    if isempty(Q1) || isempty(Q2)
        return false
    end

    shared = 0
    if length(Q1) ≤ length(Q2)
        for p in Q1
            shared += Int(p in Q2)
        end
    else
        for p in Q2
            shared += Int(p in Q1)
        end
    end
    denom = min(length(Q1), length(Q2))
    frac = denom > 0 ? shared / denom : 0.0
    return (shared ≥ Int(min_shared_abs)) && (frac ≥ float(min_shared_frac))
end

function _select_seed_faces(
    outer_faces::Vector{Face},
    axis::AbstractString,
    side::AbstractString,
    plane::Float64,
    tol_abs::Real,
)
    side_l = lowercase(side)
    seeds = Face[]
    for f in outer_faces
        fmin, fmax = _face_axis_extreme(f, axis)
        face_ext = side_l == "min" ? fmin : fmax
        if abs(face_ext - plane) ≤ tol_abs
            push!(seeds, f)
        end
    end
    return seeds
end

function _bfs_collect_boundary(
    seed_faces::Vector{Face},
    all_outer_faces::Vector{Face},
    blocks::Vector{Block},
    axis::AbstractString,
    side::AbstractString,
    plane::Float64,
    tol_abs::Real,
    node_tol_xyz::Real,
    min_shared_abs::Integer,
    min_shared_frac::Real,
)
    side_l = lowercase(side)
    function key(f::Face)
        return (f.BlockIndex, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)
    end

    on_plane = Face[]
    for f in all_outer_faces
        fmin, fmax = _face_axis_extreme(f, axis)
        v = side_l == "min" ? fmin : fmax
        opp = side_l == "min" ? fmax : fmin
        touch_plane = abs(v - plane) ≤ tol_abs
        not_past = side_l == "min" ? (opp - plane ≤ tol_abs) : (plane - opp ≤ tol_abs)
        if touch_plane && not_past
            push!(on_plane, f)
        end
    end

    pool = Dict{NTuple{7,Int},Face}()
    for f in on_plane
        pool[key(f)] = f
    end

    queue = Face[]
    for f in seed_faces
        k = key(f)
        if haskey(pool, k)
            push!(queue, pool[k])
        end
    end

    visited = Set{NTuple{7,Int}}()
    result = Face[]

    while !isempty(queue)
        cur = pop!(queue)
        kcur = key(cur)
        if kcur in visited
            continue
        end
        push!(visited, kcur)
        push!(result, cur)

        bcur = blocks[cur.BlockIndex + 1]
        for cand in values(pool)
            k2 = key(cand)
            if (k2 in visited) || k2 == kcur
                continue
            end
            b2 = blocks[cand.BlockIndex + 1]
            if touches_by_nodes(cur, cand, bcur, b2;
                    tol_xyz=node_tol_xyz,
                    min_shared_abs=min_shared_abs,
                    min_shared_frac=min_shared_frac)
                push!(queue, cand)
            end
        end
    end

    return result
end

function _rescale_faces(faces::Vector{Face}, gcd_to_use::Int)
    uniq = Dict{NTuple{7,Int},Face}()
    for f in faces
        nf = deepcopy(f)
        if gcd_to_use ≠ 1
            nf.IMIN *= gcd_to_use; nf.IMAX *= gcd_to_use
            nf.JMIN *= gcd_to_use; nf.JMAX *= gcd_to_use
            nf.KMIN *= gcd_to_use; nf.KMAX *= gcd_to_use
            nf.I *= gcd_to_use; nf.J *= gcd_to_use; nf.K *= gcd_to_use
            nf.vertices = [(v[1], v[2], v[3],
                            v[4] * gcd_to_use,
                            v[5] * gcd_to_use,
                            v[6] * gcd_to_use) for v in nf.vertices]
        end
        nf.BlockIndex = f.BlockIndex
        nf.id = f.id
        uniq[(nf.BlockIndex, nf.IMIN, nf.JMIN, nf.KMIN, nf.IMAX, nf.JMAX, nf.KMAX)] = nf
    end
    return collect(values(uniq))
end

# -----------------------------------------------------------------------------
# find_connected_faces
# -----------------------------------------------------------------------------
function find_connected_faces(
    face_to_search::Face,
    outer_faces::Vector{Face},
    connectivity_matrix::AbstractMatrix{<:Integer},
    blocks::Vector{Block},
)
    conn = copy(connectivity_matrix)
    all_matching = Face[face_to_search]
    to_search = Face[face_to_search]

    while !isempty(to_search)
        matching_faces = Face[]
        for fsrc in to_search
            n1 = normal(fsrc, blocks[fsrc.BlockIndex + 1])
            row = view(conn, fsrc.BlockIndex + 1, :)
            @inbounds for f in outer_faces
                if row[f.BlockIndex + 1] == 1
                    if length(match_indices(fsrc, f)) == 2
                        n2 = normal(f, blocks[f.BlockIndex + 1])
                        denom = norm(n1) * norm(n2)
                        if denom > 0
                            c = (n1[1]*n2[1] + n1[2]*n2[2] + n1[3]*n2[3]) / denom
                            c = clamp(c, -1.0, 1.0)
                            ang = abs(acosd(c)); ang = ang > 90 ? 180 - ang : ang
                            if ang < 30
                                row[f.BlockIndex + 1] = 0
                                conn[f.BlockIndex + 1, fsrc.BlockIndex + 1] = 0
                                push!(matching_faces, f)
                            end
                        end
                    end
                end
            end
        end
        empty!(to_search)
        for m in matching_faces
            if all(n -> !isequal(n, m), all_matching)
                push!(to_search, m)
                push!(all_matching, m)
            end
        end
    end

    return unique(all_matching)
end

# -----------------------------------------------------------------------------
# find_closest_block
# -----------------------------------------------------------------------------
function find_closest_block(
    blocks::Vector{Block},
    x::AbstractVector, y::AbstractVector, z::AbstractVector,
    centroid::NTuple{3,Float64};
    translational_direction::AbstractString = "x",
    minvalue::Bool = true,
)
    cx, cy, cz = centroid
    target_x, target_y, target_z = cx, cy, cz

    if translational_direction == "x"
        xmins = map(b -> minimum(b.X), blocks)
        xmaxs = map(b -> maximum(b.X), blocks)
        xmin, xmax = minimum(xmins), maximum(xmaxs)
        dx = xmax - xmin
        target_x = minvalue ? xmin - 0.5dx : xmax + 0.5dx
        best = Inf; besti = 0
        @inbounds for i in eachindex(x)
            d = (target_x - x[i])^2 + (cy - y[i])^2 + (cz - z[i])^2
            if d < best; best = d; besti = i; end
        end
        target_y = cy; target_z = cz
        return (besti - 1), target_x, target_y, target_z

    elseif translational_direction == "y"
        ymins = map(b -> minimum(b.Y), blocks)
        ymaxs = map(b -> maximum(b.Y), blocks)
        ymin, ymax = minimum(ymins), maximum(ymaxs)
        dy = ymax - ymin
        target_y = minvalue ? ymin - 0.5dy : ymax + 0.5dy
        best = Inf; besti = 0
        @inbounds for i in eachindex(x)
            d = (cx - x[i])^2 + (target_y - y[i])^2 + (cz - z[i])^2
            if d < best; best = d; besti = i; end
        end
        target_x = cx; target_z = cz
        return (besti - 1), target_x, target_y, target_z

    else  # "z"
        zmins = map(b -> minimum(b.Z), blocks)
        zmaxs = map(b -> maximum(b.Z), blocks)
        zmin, zmax = minimum(zmins), maximum(zmaxs)
        dz = zmax - zmin
        target_z = minvalue ? zmin - 0.5dz : zmax + 0.5dz
        best = Inf; besti = 0
        @inbounds for i in eachindex(x)
            d = (cx - x[i])^2 + (cy - y[i])^2 + (target_z - z[i])^2
            if d < best; best = d; besti = i; end
        end
        target_x = cx; target_y = cy
        return (besti - 1), target_x, target_y, target_z
    end
end

# -----------------------------------------------------------------------------
# find_bounding_faces
# -----------------------------------------------------------------------------
function find_bounding_faces(
    blocks::AbstractVector{<:Block},
    outer_faces::AbstractVector{<:AbstractDict} = Dict{String,Int}[];
    direction::AbstractString = "z",
    side::AbstractString = "both",
    tol_rel::Real = 1e-8,
    node_tol_xyz::Real = 1e-6,
    min_shared_abs::Integer = 2,
    min_shared_frac::Real = 0.005,
)
    isempty(blocks) && return Dict{String,Int}[], Dict{String,Int}[], Face[], Face[]

    axis = lowercase(direction)
    axis in ("x", "y", "z") || throw(ArgumentError("direction must be \"x\", \"y\", or \"z\" (got $direction)"))
    side_l = lowercase(side)
    side_l in ("both", "min", "max") || throw(ArgumentError("side must be \"both\", \"min\", or \"max\" (got $side)"))

    gcd_vals = Int[]
    for b in blocks
        push!(gcd_vals, gcd(b.IMAX - 1, gcd(b.JMAX - 1, b.KMAX - 1)))
    end
    gcd_to_use = max(1, minimum(gcd_vals))

    blocks_copy = Block[]
    for b in blocks
        push!(blocks_copy, deepcopy(b))
    end
    blocks_red = reduce_blocks(blocks_copy, gcd_to_use)

    outer_faces_all::Vector{Face}
    if isempty(outer_faces)
        tmp = Face[]
        for (bi, b) in enumerate(blocks_red)
            outer_b, _ = get_outer_faces(b; block_index=bi-1)
            append!(tmp, outer_b)
        end
        outer_faces_all = tmp
    else
        outer_faces_all = outer_face_dict_to_list(blocks_red, outer_faces, gcd_to_use)
    end

    gmin, gmax = _global_axis_extreme(blocks_red, axis)
    tol_abs = max(1.0, abs(gmin) + abs(gmax)) * tol_rel

    lower_faces_red = Face[]
    upper_faces_red = Face[]

    if side_l == "min" || side_l == "both"
        seeds_min = _select_seed_faces(outer_faces_all, axis, "min", gmin, tol_abs)
        lower_faces_red = _bfs_collect_boundary(
            seeds_min, outer_faces_all, blocks_red,
            axis, "min", gmin, tol_abs,
            node_tol_xyz, min_shared_abs, min_shared_frac,
        )
    end

    if side_l == "max" || side_l == "both"
        seeds_max = _select_seed_faces(outer_faces_all, axis, "max", gmax, tol_abs)
        upper_faces_red = _bfs_collect_boundary(
            seeds_max, outer_faces_all, blocks_red,
            axis, "max", gmax, tol_abs,
            node_tol_xyz, min_shared_abs, min_shared_frac,
        )
    end

    lower_faces = _rescale_faces(lower_faces_red, gcd_to_use)
    upper_faces = _rescale_faces(upper_faces_red, gcd_to_use)

    lower_export = [to_dict(f) for f in lower_faces]
    upper_export = [to_dict(f) for f in upper_faces]

    return lower_export, upper_export, lower_faces, upper_faces
end

# -----------------------------------------------------------------------------
# split_face
# -----------------------------------------------------------------------------
function split_face(face_to_split::Face, block::Block,
    imin::Int, jmin::Int, kmin::Int, imax::Int, jmax::Int, kmax::Int)

    center_face = create_face_from_diagonals(block, imin, jmin, kmin, imax, jmax, kmax)

    if kmin == kmax
        left_face   = create_face_from_diagonals(block, face_to_split.IMIN, face_to_split.JMIN, kmin, imin,  face_to_split.JMAX, kmax)
        right_face  = create_face_from_diagonals(block, imax,  face_to_split.JMIN, kmin,  face_to_split.IMAX,  face_to_split.JMAX, kmax)
        top_face    = create_face_from_diagonals(block, imin,  jmax,  kmin, imax,  face_to_split.JMAX, kmax)
        bottom_face = create_face_from_diagonals(block, imin,  face_to_split.JMIN, kmin, imax,  jmin,  kmax)
    elseif imin == imax
        left_face   = create_face_from_diagonals(block, imin,  face_to_split.JMIN, face_to_split.KMIN, imax,  jmin,  face_to_split.KMAX)
        right_face  = create_face_from_diagonals(block, imin,  jmax,  face_to_split.KMIN, imax,  face_to_split.JMAX, face_to_split.KMAX)
        top_face    = create_face_from_diagonals(block, imin,  jmin,  kmax, imax,  jmax,  face_to_split.KMAX)
        bottom_face = create_face_from_diagonals(block, imin,  jmin,  face_to_split.KMIN, imax,  jmax,  kmin)
    elseif jmin == jmax
        left_face   = create_face_from_diagonals(block, face_to_split.IMIN, jmin,  face_to_split.KMIN, imin,  jmax,  face_to_split.KMAX)
        right_face  = create_face_from_diagonals(block, imax,  jmin,  face_to_split.KMIN, face_to_split.IMAX, jmax,  face_to_split.KMAX)
        top_face    = create_face_from_diagonals(block, imin,  jmin,  kmax, imax,  jmax,  face_to_split.KMAX)
        bottom_face = create_face_from_diagonals(block, imin,  jmin,  face_to_split.KMIN, imax,  jmax,  kmin)
    else
        error("split_face: expected one pair of indices to be equal (a constant coordinate plane).")
    end

    faces = Face[top_face, bottom_face, left_face, right_face]
    faces = [f for f in faces if !is_edge(f) && !index_equals(f, center_face)]
    for f in faces
        set_block_index(f, face_to_split.BlockIndex)
    end
    return faces
end

# -----------------------------------------------------------------------------
# find_face_nearest_point
# -----------------------------------------------------------------------------
function find_face_nearest_point(faces::Vector{Face}, x::Real, y::Real, z::Real)
    best = Inf; besti = 0
    @inbounds for i in eachindex(faces)
        f = faces[i]
        d = (x - f.cx)^2 + (y - f.cy)^2 + (z - f.cz)^2
        if d < best; best = d; besti = i; end
    end
    return besti
end

# -----------------------------------------------------------------------------
# outer_face_dict_to_list / match_faces_dict_to_list
# -----------------------------------------------------------------------------
function outer_face_dict_to_list(blocks::AbstractVector{<:Block}, outer_faces::AbstractVector{<:AbstractDict}, gcd::Int=1)
    out = Face[]
    @inbounds for o in outer_faces
        face = create_face_from_diagonals(
            blocks[o["block_index"] + 1],
            Int(o["IMIN"] ÷ gcd), Int(o["JMIN"] ÷ gcd), Int(o["KMIN"] ÷ gcd),
            Int(o["IMAX"] ÷ gcd), Int(o["JMAX"] ÷ gcd), Int(o["KMAX"] ÷ gcd),
        )
        if haskey(o, "id"); face.id = o["id"]; end
        set_block_index(face, o["block_index"])
        push!(out, face)
    end
    return out
end

function match_faces_dict_to_list(blocks::AbstractVector{<:Block}, matched_faces::AbstractVector{<:AbstractDict}, gcd::Int=1)
    out = Face[]
    @inbounds for m in matched_faces
        b1 = m["block1"]; b2 = m["block2"]
        face1 = create_face_from_diagonals(
            blocks[b1["block_index"] + 1],
            Int(b1["IMIN"] ÷ gcd), Int(b1["JMIN"] ÷ gcd), Int(b1["KMIN"] ÷ gcd),
            Int(b1["IMAX"] ÷ gcd), Int(b1["JMAX"] ÷ gcd), Int(b1["KMAX"] ÷ gcd),
        )
        face2 = create_face_from_diagonals(
            blocks[b2["block_index"] + 1],
            Int(b2["IMIN"] ÷ gcd), Int(b2["JMIN"] ÷ gcd), Int(b2["KMIN"] ÷ gcd),
            Int(b2["IMAX"] ÷ gcd), Int(b2["JMAX"] ÷ gcd), Int(b2["KMAX"] ÷ gcd),
        )
        set_block_index(face1, b1["block_index"]); haskey(b1, "id") && (face1.id = b1["id"])
        set_block_index(face2, b2["block_index"]); haskey(b2, "id") && (face2.id = b2["id"])
        push!(out, face1); push!(out, face2)
    end
    return out
end

# -----------------------------------------------------------------------------
# face_matches_to_dict (helper used by periodicity & connectivity)
# -----------------------------------------------------------------------------
function face_matches_to_dict(face1::Face, face2::Face, block1::Block, block2::Block)
    match = Dict{String,Any}(
        "block1"=>Dict{String,Int}("block_index"=>face1.BlockIndex,"IMIN"=>-1,"JMIN"=>-1,"KMIN"=>-1,"IMAX"=>-1,"JMAX"=>-1,"KMAX"=>-1,"id"=>face1.id),
        "block2"=>Dict{String,Int}("block_index"=>face2.BlockIndex,"IMIN"=>-1,"JMIN"=>-1,"KMIN"=>-1,"IMAX"=>-1,"JMAX"=>-1,"KMAX"=>-1,"id"=>face2.id)
    )
    I1=(face1.IMIN,face1.IMAX); J1=(face1.JMIN,face1.JMAX); K1=(face1.KMIN,face1.KMAX)
    I2=(face2.IMIN,face2.IMAX); J2=(face2.JMIN,face2.JMAX); K2=(face2.KMIN,face2.KMAX)

    x1l = block1.X[I1[1]+1,J1[1]+1,K1[1]+1]; y1l = block1.Y[I1[1]+1,J1[1]+1,K1[1]+1]; z1l = block1.Z[I1[1]+1,J1[1]+1,K1[1]+1]
    best = (+typemax(Float64), I2[1], J2[1], K2[1])
    for p in I2, q in J2, r in K2
        dx = block2.X[p+1,q+1,r+1] - x1l
        dy = block2.Y[p+1,q+1,r+1] - y1l
        dz = block2.Z[p+1,q+1,r+1] - z1l
        d = dx*dx+dy*dy+dz*dz
        if d < best[0]
            best = (d, p,q,r)
        end
    end
    match["block1"]["IMIN"]=I1[1]; match["block1"]["JMIN"]=J1[1]; match["block1"]["KMIN"]=K1[1]
    match["block2"]["IMIN"]=best[1]; match["block2"]["JMIN"]=best[2]; match["block2"]["KMIN"]=best[3]

    x1u = block1.X[I1[2]+1,J1[2]+1,K1[2]+1]; y1u = block1.Y[I1[2]+1,J1[2]+1,K1[2]+1]; z1u = block1.Z[I1[2]+1,J1[2]+1,K1[2]+1]
    best = (+typemax(Float64), I2[1], J2[1], K2[1])
    for p in I2, q in J2, r in K2
        dx = block2.X[p+1,q+1,r+1] - x1u
        dy = block2.Y[p+1,q+1,r+1] - y1u
        dz = block2.Z[p+1,q+1,r+1] - z1u
        d = dx*dx+dy*dy+dz*dz
        if d < best[0]
            best = (d, p,q,r)
        end
    end
    match["block1"]["IMAX"]=I1[2]; match["block1"]["JMAX"]=J1[2]; match["block1"]["KMAX"]=K1[2]
    match["block2"]["IMAX"]=best[1]; match["block2"]["JMAX"]=best[2]; match["block2"]["KMAX"]=best[3]
    return match
end

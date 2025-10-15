# block_face_functions.jl — merged Block+Face helpers (no circular imports)

import LinearAlgebra: norm, dot
import Statistics: mean

# Core types we depend on
import .Block3D: Block
import .Face3D: Face, add_vertex, vertices_equals, index_equals, normal, match_indices, to_dict, set_block_index, set_face_id, is_edge
import .Utils: unique_pairs

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
# get_outer_faces
# -----------------------------------------------------------------------------
function get_outer_faces(block::Block)
    I = (0, block.IMAX - 1)
    J = (0, block.JMAX - 1)
    K = (0, block.KMAX - 1)

    faces = Face[]

    # i = I[1]
    face = Face(4)
    i = I[1]
    @inbounds for j in J, k in K
        add_vertex(face, block.X[i+1, j+1, k+1], block.Y[i+1, j+1, k+1], block.Z[i+1, j+1, k+1], i, j, k)
    end
    push!(faces, face)

    # i = I[2]
    face = Face(4)
    i = I[2]
    @inbounds for j in J, k in K
        add_vertex(face, block.X[i+1, j+1, k+1], block.Y[i+1, j+1, k+1], block.Z[i+1, j+1, k+1], i, j, k)
    end
    push!(faces, face)

    # j = J[1]
    face = Face(4)
    j = J[1]
    @inbounds for i in I, k in K
        add_vertex(face, block.X[i+1, j+1, k+1], block.Y[i+1, j+1, k+1], block.Z[i+1, j+1, k+1], i, j, k)
    end
    push!(faces, face)

    # j = J[2]
    face = Face(4)
    j = J[2]
    @inbounds for i in I, k in K
        add_vertex(face, block.X[i+1, j+1, k+1], block.Y[i+1, j+1, k+1], block.Z[i+1, j+1, k+1], i, j, k)
    end
    push!(faces, face)

    if block.KMAX > 1
        # k = K[1]
        face = Face(4)
        k = K[1]
        @inbounds for i in I, j in J
            add_vertex(face, block.X[i+1, j+1, k+1], block.Y[i+1, j+1, k+1], block.Z[i+1, j+1, k+1], i, j, k)
        end
        push!(faces, face)

        # k = K[2]
        face = Face(4)
        k = K[2]
        @inbounds for i in I, j in J
            add_vertex(face, block.X[i+1, j+1, k+1], block.Y[i+1, j+1, k+1], block.Z[i+1, j+1, k+1], i, j, k)
        end
        push!(faces, face)
    end

    matching = Tuple{Int,Int}[]
    non_matching = Face[]

    @inbounds for a in eachindex(faces)
        matchFound = false
        for b in eachindex(faces)
            if a != b && vertices_equals(faces[a], faces[b])
                push!(matching, (a, b))
                matchFound = true
            end
        end
        !matchFound && push!(non_matching, faces[a])
    end

    matching = unique_pairs(matching)
    matching_pairs = [(faces[i], faces[j]) for (i, j) in matching]
    return non_matching, matching_pairs
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
    blocks::Vector{Block},
    connectivity_matrix::AbstractMatrix{<:Integer},
    outer_faces::Vector{Dict{String,Int}} = Dict{String,Int}[];
    direction::AbstractString = "z",
)
    gcd_array = Int[]
    @inbounds for b in blocks
        push!(gcd_array, gcd(b.IMAX - 1, gcd(b.JMAX - 1, b.KMAX - 1)))
    end
    gcd_to_use = minimum(gcd_array)
    blocks_red = reduce_blocks(deepcopy(blocks), gcd_to_use)

    xyz = [(b.cx, b.cy, b.cz) for b in blocks_red]
    xs = [p[1] for p in xyz]; ys = [p[2] for p in xyz]; zs = [p[3] for p in xyz]
    cx = mean(xs); cy = mean(ys); cz = mean(zs)
    x, y, z = xs, ys, zs

    outer_faces_all::Vector{Face}
    if isempty(outer_faces)
        tmp = Face[]
        for (i, b) in enumerate(blocks_red)
            outer, _ = get_outer_faces(b)
            for o in outer
                set_block_index(o, i - 1)
                push!(tmp, o)
            end
        end
        outer_faces_all = tmp
    else
        outer_faces_all = outer_face_dict_to_list(blocks_red, outer_faces, gcd_to_use)
    end

    sel_idx_min, tx, ty, tz = find_closest_block(blocks_red, x, y, z, (cx, cy, cz); translational_direction=direction, minvalue=true)
    faces_min = [f for f in outer_faces_all if f.BlockIndex == sel_idx_min]
    @inbounds begin
        mind = Inf; mini = 0
        for (i,f) in enumerate(faces_min)
            d = (f.cx - tx)^2 + (f.cy - ty)^2 + (f.cz - tz)^2
            if d < mind; mind = d; mini = i; end
        end
        min_face = faces_min[mini]
    end

    sel_idx_max, tx, ty, tz = find_closest_block(blocks_red, x, y, z, (cx, cy, cz); translational_direction=direction, minvalue=false)
    faces_max = [f for f in outer_faces_all if f.BlockIndex == sel_idx_max]
    @inbounds begin
        mind = Inf; mini = 0
        for (i,f) in enumerate(faces_max)
            d = (f.cx - tx)^2 + (f.cy - ty)^2 + (f.cz - tz)^2
            if d < mind; mind = d; mini = i; end
        end
        max_face = faces_max[mini]
    end

    conn = copy(connectivity_matrix)
    @inbounds for i in 1:size(conn, 1); conn[i, i] = 0; end

    outer_faces_all = [o for o in outer_faces_all if o.BlockIndex != min_face.BlockIndex]
    outer_faces_all = [o for o in outer_faces_all if o.BlockIndex != max_face.BlockIndex]

    lower_connected_faces = find_connected_faces(min_face, outer_faces_all, conn, blocks_red)
    upper_connected_faces = find_connected_faces(max_face, outer_faces_all, conn, blocks_red)

    push!(lower_connected_faces, min_face)
    push!(upper_connected_faces, max_face)

    lower_connected_faces = unique(lower_connected_faces)
    upper_connected_faces = unique(upper_connected_faces)

    @inbounds for l in lower_connected_faces
        l.I *= gcd_to_use; l.J *= gcd_to_use; l.K *= gcd_to_use
    end
    @inbounds for u in upper_connected_faces
        u.I *= gcd_to_use; u.J *= gcd_to_use; u.K *= gcd_to_use
    end

    lower_connected_faces_export = [to_dict(l) for l in lower_connected_faces]
    upper_connected_faces_export = [to_dict(u) for u in upper_connected_faces]
    return lower_connected_faces_export, upper_connected_faces_export, lower_connected_faces, upper_connected_faces
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
function outer_face_dict_to_list(blocks::Vector{Block}, outer_faces::Vector{Dict{String,Int}}, gcd::Int=1)
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

function match_faces_dict_to_list(blocks::Vector{Block}, matched_faces::Vector{Dict{String,Any}}, gcd::Int=1)
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

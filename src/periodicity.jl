# periodicity.jl — functions live directly in Plot3D (no module block)

import LinearAlgebra: norm, cross
import Statistics: mean

using .Block3D: Block, recompute_centroid!
using .Face3D: Face, add_vertex, set_block_index
using .Plot3D:  get_outer_faces, faces_match, create_face_from_diagonals,
                face_matches_to_dict, outer_face_dict_to_list, reduce_blocks

# -----------------------------
# Small utilities
# -----------------------------

# Apply a 3×3 rotation matrix R to all nodes of a Block (in-place if dest===src)
function _rotate_block!(dest::Block, src::Block, R::NTuple{9,Float64})
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = R
    I,J,K = size(src)
    @inbounds for k in 1:K, j in 1:J, i in 1:I
        x = src.X[i,j,k]; y = src.Y[i,j,k]; z = src.Z[i,j,k]
        dest.X[i,j,k] = r11*x + r12*y + r13*z
        dest.Y[i,j,k] = r21*x + r22*y + r23*z
        dest.Z[i,j,k] = r31*x + r32*y + r33*z
    end
    recompute_centroid!(dest)
    return dest
end

# Translate a block by (dx,dy,dz) (in-place if dest===src)
function _translate_block!(dest::Block, src::Block, dx::Float64, dy::Float64, dz::Float64)
    I,J,K = size(src)
    @inbounds for k in 1:K, j in 1:J, i in 1:I
        dest.X[i,j,k] = src.X[i,j,k] + dx
        dest.Y[i,j,k] = src.Y[i,j,k] + dy
        dest.Z[i,j,k] = src.Z[i,j,k] + dz
    end
    recompute_centroid!(dest)
    return dest
end

# Build a Face for a canonical face-name on a block, using full span
function _face_from_name(b::Block, face_name::String)
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
    elseif face_name === "kmax"
        return create_face_from_diagonals(b, 0, 0, b.KMAX-1, b.IMAX-1, b.JMAX-1, b.KMAX-1)
    else
        error("Unknown face name: $face_name")
    end
end

# List available canonical faces for a block
_all_face_names(b::Block) = (b.KMAX > 1) ? ["imin","imax","jmin","jmax","kmin","kmax"] : ["imin","imax","jmin","jmax"]

# -----------------------------
# Public periodicity helpers
# -----------------------------

"""
    create_rotation_matrix(rotation_angle::Real; rotation_axis::AbstractString=\"x\") -> NTuple{9,Float64}

Create a 3×3 rotation matrix (row-major in a tuple) for axis \"x\"|\"y\"|\"z\".
Angle in **radians**.
"""
function create_rotation_matrix(rotation_angle::Real; rotation_axis::AbstractString="x")
    c = cos(rotation_angle); s = sin(rotation_angle)
    if rotation_axis == "x"
        return (1.0, 0.0, 0.0,
                0.0,   c,  -s,
                0.0,   s,   c)
    elseif rotation_axis == "y"
        return (  c, 0.0,   s,
                 0.0, 1.0, 0.0,
                 -s, 0.0,   c)
    elseif rotation_axis == "z"
        return (  c,  -s, 0.0,
                   s,   c, 0.0,
                 0.0, 0.0, 1.0)
    else
        error("create_rotation_matrix: rotation_axis must be \"x\", \"y\", or \"z\"")
    end
end

"""
    linear_real_transform(face1::Face, face2::Face)
        -> (ang::Float64, R::NTuple{9,Float64})

Estimate the rotation that maps `face1` onto `face2`, assuming they lie on
cylindrical surfaces about the **x-axis** (typical annulus case). Returns the
angle (radians) and a rotation matrix for that axis.

This is a pragmatic approximation: it compares the mean polar angle in the
(y,z)-plane using face centroids and returns Δθ.
"""
function linear_real_transform(face1::Face, face2::Face)
    # Compare polar angle about x-axis using face centroids
    θ1 = atan(face1.cz, face1.cy)  # atan2(y,z) but we want angle in y–z plane; use atan(z,y) or atan(y,z)? Choose atan(y,z)
    θ2 = atan(face2.cz, face2.cy)
    Δ  = θ2 - θ1
    # wrap to [-π, π]
    if Δ >  π; Δ -= 2π; end
    if Δ < -π; Δ += 2π; end
    return (Δ, create_rotation_matrix(Δ; rotation_axis="x"))
end

"""
    periodicity(blocks, outer_faces, matched_faces;
                periodic_direction=\"k\", rotation_axis=\"x\", nblades::Integer=55)

Slow but straightforward periodicity scan:
- rotates candidate faces by `2π/nblades` around `rotation_axis`
- checks corner-consistent matches with `faces_match`
- returns:
    periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all
Where the *_export lists are dicts matching your Python schema.
"""
function periodicity(blocks::Vector{Block},
                     outer_faces::Vector{Dict{String,Int}},
                     matched_faces::Vector{Dict{String,Any}};
                     periodic_direction::AbstractString="k",
                     rotation_axis::AbstractString="x",
                     nblades::Integer=55,
                     tol::Real=1e-8)

    θ = 2π / nblades
    R = create_rotation_matrix(θ; rotation_axis=rotation_axis)

    # Build outer faces as Face objects
    faces_outer = outer_face_dict_to_list(blocks, outer_faces, 1)

    # Remove any faces that are already in `matched_faces` (convert both ends)
    matched_as_faces = let tmp = Face[]
        for m in matched_faces
            b1 = m["block1"]; b2 = m["block2"]
            push!(tmp, create_face_from_diagonals(blocks[b1["block_index"]+1], b1["IMIN"], b1["JMIN"], b1["KMIN"], b1["IMAX"], b1["JMAX"], b1["KMAX"]))
            push!(tmp, create_face_from_diagonals(blocks[b2["block_index"]+1], b2["IMIN"], b2["JMIN"], b2["KMIN"], b2["IMAX"], b2["JMAX"], b2["KMAX"]))
        end
        tmp
    end

    faces_outer = [f for f in faces_outer if all(g -> !((f.IMIN==g.IMIN)&&(f.IMAX==g.IMAX)&&(f.JMIN==g.JMIN)&&(f.JMAX==g.JMAX)&&(f.KMIN==g.KMIN)&&(f.KMAX==g.KMAX)&& (f.BlockIndex==g.BlockIndex)), matched_as_faces)]

    # Build rotated copy of all blocks
    rot_blocks = [Block(copy(b.X), copy(b.Y), copy(b.Z)) for b in blocks]
    for i in eachindex(blocks)
        _rotate_block!(rot_blocks[i], blocks[i], R)
    end

    periodic_pairs = Tuple{Face,Face}[]

    # Compare each outer face to the rotated geometry:
    for f in faces_outer
        # its counterpart should be on the **same** block after rotation if this is intra-block periodic;
        # however, in multi-block, allow any block.
        for (j, brot) in enumerate(rot_blocks)
            # try all canonical face names of rotated block
            for name in _all_face_names(brot)
                f2 = _face_from_name(brot, name)
                ok, _ = faces_match(
                    _arrays_from_face(blocks[f.BlockIndex+1], f),
                    _arrays_from_face(brot, f2);
                    tol=tol
                )
                if ok
                    # annotate block index fields
                    set_block_index(f, f.BlockIndex)
                    set_block_index(f2, j-1)
                    push!(periodic_pairs, (f, f2))
                end
            end
        end
    end

    # Unique by indices to avoid dupes
    function _key(face::Face)
        return (face.BlockIndex, face.IMIN, face.JMIN, face.KMIN, face.IMAX, face.JMAX, face.KMAX)
    end
    seen = Set{Tuple{Int,Int,Int,Int,Int,Int,Int}}()
    uniq_pairs = Tuple{Face,Face}[]
    for (a,b) in periodic_pairs
        k = (_key(a), _key(b))
        if !(k in seen)
            push!(seen, k)
            push!(uniq_pairs, (a,b))
        end
    end

    # Exports
    periodic_faces_export = [face_matches_to_dict(a,b, blocks[a.BlockIndex+1], blocks[b.BlockIndex+1]) for (a,b) in uniq_pairs]
    periodic_faces = uniq_pairs

    # faces not matched are still outer
    matched_set = Set{Tuple{Int,Int,Int,Int,Int,Int,Int}}((_key(a)) for (a,_) in uniq_pairs)
    outer_faces_all = [f for f in faces_outer if !((_key(f)) in matched_set)]
    outer_faces_export = [Dict(
        "block_index"=>f.BlockIndex,
        "IMIN"=>f.IMIN,"JMIN"=>f.JMIN,"KMIN"=>f.KMIN,
        "IMAX"=>f.IMAX,"JMAX"=>f.JMAX,"KMAX"=>f.KMAX,
        "id"=>f.id
    ) for f in outer_faces_all]

    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all
end

# helper: extract 2D arrays for a Face on a Block (copied from get_face_intersection style)
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

"""
    periodicity_fast(blocks, outer_faces, matched_faces; periodic_direction=\"k\", rotation_axis=\"x\", nblades=55)

Same return signature as `periodicity`, but first down-samples the mesh by the
minimum gcd across blocks for speed (then maps indices back).
"""
function periodicity_fast(blocks::Vector{Block},
                          outer_faces::Vector{Dict{String,Int}},
                          matched_faces::Vector{Dict{String,Any}};
                          periodic_direction::AbstractString="k",
                          rotation_axis::AbstractString="x",
                          nblades::Integer=55,
                          tol::Real=1e-8)

    gcds = [gcd(b.IMAX-1, gcd(b.JMAX-1, b.KMAX-1)) for b in blocks]
    g    = minimum(gcds) > 0 ? minimum(gcds) : 1
    blocks_red = reduce_blocks(deepcopy(blocks), g)

    pf_exp, of_exp, pf, of = periodicity(blocks_red, outer_faces, matched_faces;
                                         periodic_direction=periodic_direction,
                                         rotation_axis=rotation_axis,
                                         nblades=nblades, tol=tol)

    # map reduced-index faces back to original by multiplying index bounds by g
    function _scale!(d::Dict{String,Int})
        d["IMIN"] *= g; d["IMAX"] *= g
        d["JMIN"] *= g; d["JMAX"] *= g
        d["KMIN"] *= g; d["KMAX"] *= g
        return d
    end
    pf_exp = [_scale!(copy(d)) for d in pf_exp]
    of_exp = [_scale!(copy(d)) for d in of_exp]

    return pf_exp, of_exp, pf, of
end

"""
    rotated_periodicity(blocks, matched_faces, outer_faces;
                        rotation_angle::Real, rotation_axis::AbstractString=\"x\", ReduceMesh::Bool=true)

Rotate a copy of the geometry by `rotation_angle` (degrees) around the given axis.
Return same tuple as `periodicity`.
"""
function rotated_periodicity(blocks::Vector{Block},
                             matched_faces::Vector{Dict{String,Any}},
                             outer_faces::Vector{Dict{String,Int}};
                             rotation_angle::Real,
                             rotation_axis::AbstractString="x",
                             ReduceMesh::Bool=true,
                             tol::Real=1e-8)

    # optional reduction for speed
    if ReduceMesh
        gcds = [gcd(b.IMAX-1, gcd(b.JMAX-1, b.KMAX-1)) for b in blocks]
        g    = minimum(gcds) > 0 ? minimum(gcds) : 1
        blocks = reduce_blocks(deepcopy(blocks), g)
    end

    R = create_rotation_matrix(deg2rad(rotation_angle); rotation_axis=rotation_axis)
    rot_blocks = [Block(copy(b.X), copy(b.Y), copy(b.Z)) for b in blocks]
    for i in eachindex(blocks)
        _rotate_block!(rot_blocks[i], blocks[i], R)
    end

    faces_outer = outer_face_dict_to_list(blocks, outer_faces, 1)
    periodic_pairs = Tuple{Face,Face}[]

    # Compare outer faces vs rotated geometry
    for f in faces_outer
        for (j, brot) in enumerate(rot_blocks)
            for name in _all_face_names(brot)
                f2 = _face_from_name(brot, name)
                ok, _ = faces_match(_arrays_from_face(blocks[f.BlockIndex+1], f), _arrays_from_face(brot, f2); tol=tol)
                if ok
                    set_block_index(f, f.BlockIndex)
                    set_block_index(f2, j-1)
                    push!(periodic_pairs, (f, f2))
                end
            end
        end
    end

    periodic_faces_export = [face_matches_to_dict(a,b, blocks[a.BlockIndex+1], blocks[b.BlockIndex+1]) for (a,b) in periodic_pairs]
    matched_set = Set{Tuple{Int,Int,Int,Int,Int,Int,Int}}((f.BlockIndex,f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX) for (f,_) in periodic_pairs)
    outer_faces_all = [f for f in faces_outer if !((f.BlockIndex,f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX) in matched_set)]
    outer_faces_export = [Dict("block_index"=>f.BlockIndex,"IMIN"=>f.IMIN,"JMIN"=>f.JMIN,"KMIN"=>f.KMIN,"IMAX"=>f.IMAX,"JMAX"=>f.JMAX,"KMAX"=>f.KMAX,"id"=>f.id) for f in outer_faces_all]

    return periodic_faces_export, outer_faces_export, periodic_pairs, outer_faces_all
end
# periodicity.jl — replace your current translational_periodicity with this one

"""
    translational_periodicity(
        blocks::Vector{Block},
        outer_faces::Vector{Dict{String,Int}};
        delta::Union{Nothing,Real}=nothing,
        translational_direction::AbstractString="z",
        node_tol_xyz::Union{Nothing,Real}=nothing,
        min_shared_frac::Real = 0.02,
        min_shared_abs::Int  = 4,
        stride_u::Int = 1,
        stride_v::Int = 1,
    ) -> Tuple{
            Vector{Dict{String,Any}},                  # periodic_faces_export
            Vector{Tuple{Face,Face,Dict{String,String}}}, # periodic_pairs (Face,Face,mapping)
            Vector{Dict{String,Int}}                   # outer_faces_remaining
       }

    Detect translational periodicity between block faces along a given axis (`"x"`, `"y"`, or `"z"`).

    This mirrors the Python implementation:

    1. **Bounded faces:** Calls `find_bounding_faces(blocks, outer_faces, axis, "both")` to get
    **lower** and **upper** candidate faces (as dicts) for the requested axis.
    2. **GCD reduction:** Reduces all blocks to a common index scale via the minimum
    `gcd(IMAX-1, JMAX-1, KMAX-1)` so face index ranges become compatible across blocks.
    3. **Shifted copies:** Creates *translated* copies of the reduced blocks by ±Δ along the axis,
    where Δ is inferred from global min/max of coordinates on that axis if not provided.
    4. **Adaptive tolerance:** For each candidate face pair, computes a per-pair absolute node tolerance
    from the face’s in-plane median edge length (≈ 3% of the larger in-plane spacing), unless
    `node_tol_xyz` is supplied.
    5. **Orthogonal-plane precheck:** Before doing a node-by-node match, projects the two faces into the
    plane orthogonal to the periodic axis (after shifting one face by +Δ), snaps points to a tolerance
    grid, and requires sufficient overlap by **count** and **fraction** (`min_shared_abs` and
    `min_shared_frac`). This avoids expensive matches when faces are clearly incompatible.
    6. **Node match:** Uses the existing `faces_match` (array-based) routine to check node sharing with
    the adaptive tolerance, trying both “lower shifted up vs upper original” and the symmetric guard
    (“lower original vs upper shifted down”).
    7. **Pairing:** Greedily pairs each lower face with the first upper face that matches; records a
    per-axis **min→min / min→max** index mapping for I/J/K similar to Python’s `mapping_minmax`.
    8. **Scale back:** Scales all reported face index bounds back up by the GCD factor to match the
    original mesh resolution.
    9. **Prune outer faces:** Removes matched periodic faces from the `outer_faces` list (preserving
    any custom `id` fields left on the remaining ones).
    10. **Return:**
        - `periodic_faces_export`: JSON-like dicts (block indices, IJK extents, mapping, match mode)
        - `periodic_pairs`: `(Face, Face, Dict{String,String})` at original index scale
        - `outer_faces_remaining`: updated `outer_faces` (dicts) with periodic ones removed

    Arguments
    ---------
    - `blocks`: vector of `Block`s.
    - `outer_faces`: outer faces as dictionaries (`"block_index"`, `"IMIN"`, `"JMIN"`, `"KMIN"`, `"IMAX"`, `"JMAX"`, `"KMAX"`).
    - `delta`: periodic spacing along the chosen axis; if `nothing`, inferred from global bounds.
    - `translational_direction`: `"x"`, `"y"`, or `"z"` (default `"z"`).
    - `node_tol_xyz`: absolute coordinate tolerance for node matching. If not set, an adaptive tolerance
    is computed per pair.
    - `min_shared_frac`: minimum fraction of shared nodes in the precheck / match.
    - `min_shared_abs`: minimum absolute number of shared nodes.
    - `stride_u`, `stride_v`: optional subsampling strides passed down to face matching.

    Notes
    -----
    - Requires `find_bounding_faces`, `outer_face_dict_to_list`, `reduce_blocks`, and `faces_match` to be available.
    - Assumes `faces_match` accepts `(X1,Y1,Z1)` and `(X2,Y2,Z2)` array views and supports `tol`, `stride_u`,
    `stride_v`, `min_shared_frac`, and `min_shared_abs` keywords (align with your existing signature).
"""
function translational_periodicity(
    blocks::Vector{Block},
    outer_faces::Vector{Dict{String,Int}};
    delta::Union{Nothing,Real}=nothing,
    translational_direction::AbstractString="z",
    node_tol_xyz::Union{Nothing,Real}=nothing,
    min_shared_frac::Real = 0.02,
    min_shared_abs::Int  = 4,
    stride_u::Int = 1,
    stride_v::Int = 1,
)
    axis = lowercase(translational_direction)
    @assert axis == "x" || axis == "y" || axis == "z"

    # --- 0) find lower/upper connected faces for this axis (Python parity) ----
    lower_connected_faces, upper_connected_faces, _, _ =
        find_bounding_faces(blocks, outer_faces, translational_direction, "both")

    # --- 1) GCD reduce mesh (indexes become compatible across blocks) ---------
    gcds = [gcd(b.IMAX-1, gcd(b.JMAX-1, b.KMAX-1)) for b in blocks]
    g    = maximum([1, minimum(gcds)])  # guard against 0
    blocks_r = reduce_blocks(deepcopy(blocks), g)

    lower_faces_r = outer_face_dict_to_list(blocks, lower_connected_faces, g)
    upper_faces_r = outer_face_dict_to_list(blocks, upper_connected_faces, g)

    # --- 2) infer Δ along axis if not provided --------------------------------
    if delta === nothing
        if axis == "x"
            a_min = minimum(b -> minimum(b.X), blocks_r)
            a_max = maximum(b -> maximum(b.X), blocks_r)
        elseif axis == "y"
            a_min = minimum(b -> minimum(b.Y), blocks_r)
            a_max = maximum(b -> maximum(b.Y), blocks_r)
        else # "z"
            a_min = minimum(b -> minimum(b.Z), blocks_r)
            a_max = maximum(b -> maximum(b.Z), blocks_r)
        end
        delta = a_max - a_min
    end
    Δ = float(delta)

    # --- 3) build shifted copies (up and down) --------------------------------
    function _shifted(bb::Vector{Block}, amount::Float64)
        cp = [Block(copy(b.X), copy(b.Y), copy(b.Z)) for b in bb]
        if axis == "x"
            for i in eachindex(bb); _translate_block!(cp[i], bb[i], amount, 0.0, 0.0); end
        elseif axis == "y"
            for i in eachindex(bb); _translate_block!(cp[i], bb[i], 0.0, amount, 0.0); end
        else
            for i in eachindex(bb); _translate_block!(cp[i], bb[i], 0.0, 0.0, amount); end
        end
        return cp
    end
    blocks_up = _shifted(blocks_r, +Δ)
    blocks_dn = _shifted(blocks_r, -Δ)

    # quick selectors
    @inline B(which::Symbol, i::Int) = which === :orig ? blocks_r[i] :
                                       which === :up   ? blocks_up[i] :
                                                         blocks_dn[i]

    # --- helpers for adaptive tolerance & precheck ----------------------------
    # median in-plane spacing of a face
    function _median_inplane_spacing(f::Face, b::Block)
        I0,I1,J0,J1,K0,K1 = f.IMIN, f.IMAX, f.JMIN, f.JMAX, f.KMIN, f.KMAX
        xs = Float64[]  # collect edge lengths
        if I0 == I1        # vary (J,K)
            i = I0 + 1
            X = view(b.X, i, J0+1:J1+1, K0+1:K1+1)
            Y = view(b.Y, i, J0+1:J1+1, K0+1:K1+1)
            Z = view(b.Z, i, J0+1:J1+1, K0+1:K1+1)
        elseif J0 == J1    # vary (I,K)
            j = J0 + 1
            X = view(b.X, I0+1:I1+1, j, K0+1:K1+1)
            Y = view(b.Y, I0+1:I1+1, j, K0+1:K1+1)
            Z = view(b.Z, I0+1:I1+1, j, K0+1:K1+1)
        else               # vary (I,J)
            k = K0 + 1
            X = view(b.X, I0+1:I1+1, J0+1:J1+1, k)
            Y = view(b.Y, I0+1:I1+1, J0+1:J1+1, k)
            Z = view(b.Z, I0+1:I1+1, J0+1:J1+1, k)
        end
        # along dim1
        if size(X,1) > 1
            for a in axes(X,2), c in axes(X,3)
                for u in 1:size(X,1)-1
                    dx = X[u+1,a,c]-X[u,a,c]; dy = Y[u+1,a,c]-Y[u,a,c]; dz = Z[u+1,a,c]-Z[u,a,c]
                    push!(xs, sqrt(dx*dx+dy*dy+dz*dz))
                end
            end
        end
        # along dim2
        if size(X,2) > 1
            for u in axes(X,1), c in axes(X,3)
                for a in 1:size(X,2)-1
                    dx = X[u,a+1,c]-X[u,a,c]; dy = Y[u,a+1,c]-Y[u,a,c]; dz = Z[u,a+1,c]-Z[u,a,c]
                    push!(xs, sqrt(dx*dx+dy*dy+dz*dz))
                end
            end
        end
        isempty(xs) ? 1.0 : Statistics.median(xs)
    end

    @inline function _pair_tol(fA::Face, fB::Face)
        node_tol_xyz !== nothing && return float(node_tol_xyz)
        sA = _median_inplane_spacing(fA, B(:orig, fA.BlockIndex+1))
        sB = _median_inplane_spacing(fB, B(:orig, fB.BlockIndex+1))
        return max(0.03 * max(sA, sB), 1e-4)  # ~3% of in-plane spacing
    end

    # very fast orthogonal-plane precheck (projects away the periodic axis)
    function _orthogonal_precheck(fA::Face, fB::Face, d::Float64, tol::Float64)
        XA,YA,ZA = _arrays_from_face(B(:orig, fA.BlockIndex+1), fA)
        XB,YB,ZB = _arrays_from_face(B(:orig, fB.BlockIndex+1), fB)
        # shift A by +Δ along axis
        if axis == "x"
            XA = XA .+ d
            # compare (y,z)
            A1, A2 = YA, ZA; B1, B2 = YB, ZB
        elseif axis == "y"
            YA = YA .+ d
            A1, A2 = XA, ZA; B1, B2 = XB, ZB
        else
            ZA = ZA .+ d
            A1, A2 = XA, YA; B1, B2 = XB, YB
        end
        # hash points by tol grid
        function _codes(U::AbstractArray, V::AbstractArray)
            n1,n2 = size(U)
            out = Tuple{Int,Int}[]
            invt = 1.0 / tol
            @inbounds for j in 1:n2, i in 1:n1
                push!(out, (floor(Int, U[i,j]*invt), floor(Int, V[i,j]*invt)))
            end
            out
        end
        SA = Set(_codes(A1, A2))
        SB = Set(_codes(B1, B2))
        inter = intersect(SA, SB)
        return length(inter) >= max(min_shared_abs, Int(floor(min_shared_frac * min(length(SA), length(SB)))))
    end

    # node-sharing style matcher built on faces_match (array-based)
    function _faces_match_pair(fL::Face, fU::Face)
        tol_pair = _pair_tol(fL, fU)
        # orthogonal-plane guards
        if _orthogonal_precheck(fL, fU, Δ, tol_pair); return true; end
        ok, _ = faces_match(_arrays_from_face(B(:orig, fL.BlockIndex+1), fL),
                            _arrays_from_face(B(:up,   fU.BlockIndex+1), fU);
                            tol = tol_pair, stride_u=stride_u, stride_v=stride_v,
                            min_shared_frac=min_shared_frac, min_shared_abs=min_shared_abs)
        ok && return true
        ok2, _ = faces_match(_arrays_from_face(B(:dn,  fL.BlockIndex+1), fL),
                             _arrays_from_face(B(:orig, fU.BlockIndex+1), fU);
                             tol = tol_pair, stride_u=stride_u, stride_v=stride_v,
                             min_shared_frac=min_shared_frac, min_shared_abs=min_shared_abs)
        return ok2
    end

    # index mapping (min->min or min->max, like Python)
    function _mapping_minmax(fA::Face, fB::Face)
        out = Dict{String,String}()
        for (Ax,Amin,Amax,Bmin,Bmax) in (("I",fA.IMIN,fA.IMAX,fB.IMIN,fB.IMAX),
                                         ("J",fA.JMIN,fA.JMAX,fB.JMIN,fB.JMAX),
                                         ("K",fA.KMIN,fA.KMAX,fB.KMIN,fB.KMAX))
            if (Amin == Bmin) && (Amax == Bmax)
                out[Ax] = "min->min"
            elseif (Amin == Bmax) && (Amax == Bmin)
                out[Ax] = "min->max"
            else
                d_mm = abs(Amin-Bmin)+abs(Amax-Bmax)
                d_mM = abs(Amin-Bmax)+abs(Amax-Bmin)
                out[Ax] = d_mm <= d_mM ? "min->min" : "min->max"
            end
        end
        out
    end

    # pools to pair greedily
    lower_pool = collect(Dict{Tuple,Face}(( (f.BlockIndex,f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX) => f ) for f in lower_faces_r) |> values)
    upper_pool = collect(Dict{Tuple,Face}(( (f.BlockIndex,f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX) => f ) for f in upper_faces_r) |> values)

    periodic_pairs_r = Tuple{Face,Face,Dict{String,String}}[]
    periodic_export  = Vector{Dict{String,Any}}()

    for fL in lower_pool
        matched = false
        for (j,fU) in enumerate(upper_pool)
            if _faces_match_pair(fL, fU)
                m = _mapping_minmax(fL, fU)
                push!(periodic_pairs_r, (fL, fU, m))
                push!(periodic_export, Dict(
                    "block1" => Dict(
                        "block_index"=>fL.BlockIndex,
                        "IMIN"=>fL.IMIN, "JMIN"=>fL.JMIN, "KMIN"=>fL.KMIN,
                        "IMAX"=>fL.IMAX, "JMAX"=>fL.JMAX, "KMAX"=>fL.KMAX),
                    "block2" => Dict(
                        "block_index"=>fU.BlockIndex,
                        "IMIN"=>fU.IMIN, "JMIN"=>fU.JMIN, "KMIN"=>fU.KMIN,
                        "IMAX"=>fU.IMAX, "JMAX"=>fU.JMAX, "KMAX"=>fU.KMAX),
                    "mapping" => m,
                    "mode"    => "$(axis)_precheck/match"
                ))
                deleteat!(upper_pool, j)
                matched = true
                break
            end
        end
        matched || nothing
    end

    # --- 9) scale indices back up by g ----------------------------------------
    for rec in periodic_export
        for side in ("block1","block2")
            d = rec[side]::Dict{String,Int}
            d["IMIN"] *= g; d["IMAX"] *= g
            d["JMIN"] *= g; d["JMAX"] *= g
            d["KMIN"] *= g; d["KMAX"] *= g
        end
    end
    periodic_pairs = Tuple{Face,Face,Dict{String,String}}[]
    for (fL,fU,m) in periodic_pairs_r
        gL = deepcopy(fL); gU = deepcopy(fU)
        gL.I *= g; gL.J *= g; gL.K *= g
        gU.I *= g; gU.J *= g; gU.K *= g
        push!(periodic_pairs, (gL, gU, m))
    end

    # --- 10) remove periodic faces from outer_faces ---------------------------
    periodic_keys = Set{NTuple{7,Int}}()
    for rec in periodic_export
        for side in ("block1","block2")
            d = rec[side]::Dict{String,Int}
            push!(periodic_keys, (d["block_index"], d["IMIN"], d["JMIN"], d["KMIN"],
                                              d["IMAX"], d["JMAX"], d["KMAX"]))
        end
    end
    outer_faces_remaining = [o for o in outer_faces if
        (o["block_index"], o["IMIN"], o["JMIN"], o["KMIN"], o["IMAX"], o["JMAX"], o["KMAX"]) ∉ periodic_keys]

    return periodic_export, periodic_pairs, outer_faces_remaining
end

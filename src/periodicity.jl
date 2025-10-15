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

"""
    translational_periodicity(blocks, lower_connected_faces, upper_connected_faces;
                              delta::Union{Nothing,Real}=nothing, translational_direction::AbstractString=\"z\")

Copy the mesh and translate it along the specified axis by `delta`
(estimated from bounds if not provided), then match faces.

Returns: periodic_faces_export, outer_faces_export, periodic_pairs, outer_faces_all
"""
function translational_periodicity(blocks::Vector{Block},
                                   lower_connected_faces::Vector{Dict{String,Int}},
                                   upper_connected_faces::Vector{Dict{String,Int}};
                                   delta::Union{Nothing,Real}=nothing,
                                   translational_direction::AbstractString="z",
                                   tol::Real=1e-8)

    # Estimate delta from block bounds if not supplied
    if delta === nothing
        if translational_direction == "x"
            xmin = minimum(b->minimum(b.X), blocks); xmax = maximum(b->maximum(b.X), blocks)
            delta = xmax - xmin
        elseif translational_direction == "y"
            ymin = minimum(b->minimum(b.Y), blocks); ymax = maximum(b->maximum(b.Y), blocks)
            delta = ymax - ymin
        else
            zmin = minimum(b->minimum(b.Z), blocks); zmax = maximum(b->maximum(b.Z), blocks)
            delta = zmax - zmin
        end
    end
    δ = Float64(delta)

    # Make a translated copy of the whole mesh
    tblocks = [Block(copy(b.X), copy(b.Y), copy(b.Z)) for b in blocks]
    if translational_direction == "x"
        for i in eachindex(blocks); _translate_block!(tblocks[i], blocks[i], δ, 0.0, 0.0); end
    elseif translational_direction == "y"
        for i in eachindex(blocks); _translate_block!(tblocks[i], blocks[i], 0.0, δ, 0.0); end
    else
        for i in eachindex(blocks); _translate_block!(tblocks[i], blocks[i], 0.0, 0.0, δ); end
    end

    # Build face lists from dicts
    lower_faces = outer_face_dict_to_list(blocks, lower_connected_faces, 1)
    upper_faces = outer_face_dict_to_list(blocks, upper_connected_faces, 1)

    periodic_pairs = Tuple{Face,Face}[]

    # For each lower face on original mesh, find matching upper face on translated mesh
    for f in lower_faces
        for (j, bt) in enumerate(tblocks)
            for name in _all_face_names(bt)
                f2 = _face_from_name(bt, name)
                ok, _ = faces_match(_arrays_from_face(blocks[f.BlockIndex+1], f), _arrays_from_face(bt, f2); tol=tol)
                if ok
                    set_block_index(f, f.BlockIndex); set_block_index(f2, j-1)
                    push!(periodic_pairs, (f, f2))
                end
            end
        end
    end

    periodic_faces_export = [face_matches_to_dict(a,b, blocks[a.BlockIndex+1], tblocks[b.BlockIndex+1]) for (a,b) in periodic_pairs]

    matched_set = Set{Tuple{Int,Int,Int,Int,Int,Int,Int}}((f.BlockIndex,f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX) for (f,_) in periodic_pairs)
    outer_faces_all = [f for f in vcat(lower_faces, upper_faces) if !((f.BlockIndex,f.IMIN,f.JMIN,f.KMIN,f.IMAX,f.JMAX,f.KMAX) in matched_set)]
    outer_faces_export = [Dict("block_index"=>f.BlockIndex,"IMIN"=>f.IMIN,"JMIN"=>f.JMIN,"KMIN"=>f.KMIN,"IMAX"=>f.IMAX,"JMAX"=>f.JMAX,"KMAX"=>f.KMAX,"id"=>f.id) for f in outer_faces_all]

    return periodic_faces_export, outer_faces_export, periodic_pairs, outer_faces_all
end

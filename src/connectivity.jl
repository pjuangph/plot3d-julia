# connectivity.jl — explicit imports, no DataFrames

import LinearAlgebra: norm
import Statistics: mean

import .Block3D: Block
import .Face3D: Face
# Do NOT import a module named "facefunctions" — those functions are in this same Plot3D module.
# Just call them directly, but ensure facefunctions.jl is included before this file.
import .Utils: unique_pairs


# -----------------------------------------------------------------------------
# Small container for match results (no DataFrames)
# -----------------------------------------------------------------------------
struct FaceMatchSet
    i1::Vector{Int}; j1::Vector{Int}; k1::Vector{Int}
    i2::Vector{Int}; j2::Vector{Int}; k2::Vector{Int}
end
FaceMatchSet() = FaceMatchSet(Int[], Int[], Int[], Int[], Int[], Int[])
npoints(m::FaceMatchSet) = length(m.i1)
isempty(m::FaceMatchSet) = npoints(m) == 0

function push_match!(m::FaceMatchSet, i1::Int,j1::Int,k1::Int,i2::Int,j2::Int,k2::Int)
    push!(m.i1, i1); push!(m.j1, j1); push!(m.k1, k1)
    push!(m.i2, i2); push!(m.j2, j2); push!(m.k2, k2)
    return m
end

function _mask!(m::FaceMatchSet, keep::Vector{Bool})
    m.i1 = m.i1[keep]; m.j1 = m.j1[keep]; m.k1 = m.k1[keep]
    m.i2 = m.i2[keep]; m.j2 = m.j2[keep]; m.k2 = m.k2[keep]
    return m
end

# -----------------------------------------------------------------------------
# point_match — search (i,j) location of closest X2,Y2,Z2 to (x,y,z)
# -----------------------------------------------------------------------------
function point_match(x::Real, y::Real, z::Real,
                     X2::AbstractMatrix, Y2::AbstractMatrix, Z2::AbstractMatrix; tol::Real=1e-6)
    best = Inf; best_i = 0; best_j = 0
    @inbounds for j in axes(X2, 2)
        @simd for i in axes(X2, 1)
            dx = x - X2[i, j]; dy = y - Y2[i, j]; dz = z - Z2[i, j]
            d = sqrt(dx*dx + dy*dy + dz*dz)
            if d < best
                best = d; best_i = i; best_j = j
            end
        end
    end
    return best < tol ? (best_i - 1, best_j - 1) : (-1, -1)  # return 0-based like Python
end

# -----------------------------------------------------------------------------
# select_multi_dimensional — slice a face (constant i or j or k)
# inputs are 0-based (Python-style); views are 1-based
# -----------------------------------------------------------------------------
function select_multi_dimensional(T::AbstractArray,
                                  dim1::Tuple{Int,Int},
                                  dim2::Tuple{Int,Int},
                                  dim3::Tuple{Int,Int})
    i1,i2 = dim1; j1,j2 = dim2; k1,k2 = dim3
    if i1 == i2
        return @view T[i1+1, j1+1:j2+1, k1+1:k2+1]
    elseif j1 == j2
        return @view T[i1+1:i2+1, j1+1, k1+1:k2+1]
    elseif k1 == k2
        return @view T[i1+1:i2+1, j1+1:j2+1, k1+1]
    else
        return @view T[i1+1:i2+1, j1+1:j2+1, k1+1:k2+1]
    end
end

# -----------------------------------------------------------------------------
# Edge check and increasing filters (ported logic, vectorized)
# -----------------------------------------------------------------------------
function __check_edge(m::FaceMatchSet)
    imin = minimum(m.i1); jmin = minimum(m.j1); kmin = minimum(m.k1)
    imax = maximum(m.i1); jmax = maximum(m.j1); kmax = maximum(m.k1)
    edge_matches = 0
    edge_matches += (imin == imax) ? 1 : 0
    edge_matches += (jmin == jmax) ? 1 : 0
    edge_matches += (kmin == kmax) ? 1 : 0
    return edge_matches >= 2
end

function __filter_block_increasing!(m::FaceMatchSet, key::Symbol)
    vals = unique(getfield(m, key))
    sort!(vals)
    if length(vals) <= 1
        # This indicates a degenerate (edge-like) set on this axis; clear results
        m.i1 = Int[]; m.j1 = Int[]; m.k1 = Int[]
        m.i2 = Int[]; m.j2 = Int[]; m.k2 = Int[]
        return m
    end
    keep_set = Set{Int}()
    for idx in 1:length(vals)-1
        if (vals[idx+1] - vals[idx]) == 1
            push!(keep_set, vals[idx])
        end
    end
    if (vals[end] - vals[end-1]) == 1
        push!(keep_set, vals[end])
    end
    keep = [getfield(m, key)[t] in keep_set for t in eachindex(getfield(m,key))]
    return _mask!(m, keep)
end

# -----------------------------------------------------------------------------
# get_face_intersection — core matching of two faces (from different blocks)
# -----------------------------------------------------------------------------
function get_face_intersection(face1::Face, face2::Face, block1::Block, block2::Block; tol::Real=1e-6)
    matches = FaceMatchSet()
    split_faces1 = Face[]; split_faces2 = Face[]

    I1 = (face1.IMIN, face1.IMAX); J1 = (face1.JMIN, face1.JMAX); K1 = (face1.KMIN, face1.KMAX)
    I2 = (face2.IMIN, face2.IMAX); J2 = (face2.JMIN, face2.JMAX); K2 = (face2.KMIN, face2.KMAX)

    X1 = select_multi_dimensional(block1.X, I1, J1, K1)
    Y1 = select_multi_dimensional(block1.Y, I1, J1, K1)
    Z1 = select_multi_dimensional(block1.Z, I1, J1, K1)

    X2 = select_multi_dimensional(block2.X, I2, J2, K2)
    Y2 = select_multi_dimensional(block2.Y, I2, J2, K2)
    Z2 = select_multi_dimensional(block2.Z, I2, J2, K2)

    if I1[1] == I1[2]           # Face 1 constant-i
        @inbounds for p in axes(X1,1), q in axes(X1,2)
            pm, qm = point_match(X1[p,q], Y1[p,q], Z1[p,q], X2, Y2, Z2; tol=tol)
            if pm != -1
                (I2[1]==I2[2]) && push_match!(matches, I1[1], (p-1)+J1[1], (q-1)+K1[1], I2[1], pm+J2[1], qm+K2[1])
                (J2[1]==J2[2]) && push_match!(matches, I1[1], (p-1)+J1[1], (q-1)+K1[1], pm+I2[1], J2[1], qm+K2[1])
                (K2[1]==K2[2]) && push_match!(matches, I1[1], (p-1)+J1[1], (q-1)+K1[1], pm+I2[1], qm+J2[1], K2[1])
            end
        end

    elseif J1[1] == J1[2]       # Face 1 constant-j
        @inbounds for p in axes(X1,1), q in axes(X1,2)
            pm, qm = point_match(X1[p,q], Y1[p,q], Z1[p,q], X2, Y2, Z2; tol=tol)
            if pm != -1
                (I2[1]==I2[2]) && push_match!(matches, (p-1)+I1[1], J1[1], (q-1)+K1[1], I2[1], pm+J2[1], qm+K2[1])
                (J2[1]==J2[2]) && push_match!(matches, (p-1)+I1[1], J1[1], (q-1)+K1[1], pm+I2[1], J2[1], qm+K2[1])
                (K2[1]==K2[2]) && push_match!(matches, (p-1)+I1[1], J1[1], (q-1)+K1[1], pm+I2[1], qm+J2[1], K2[1])
            end
        end

    elseif K1[1] == K1[2]       # Face 1 constant-k
        @inbounds for p in axes(X1,1), q in axes(X1,2)
            pm, qm = point_match(X1[p,q], Y1[p,q], Z1[p,q], X2, Y2, Z2; tol=tol)
            if pm != -1
                (I2[1]==I2[2]) && push_match!(matches, (p-1)+I1[1], (q-1)+J1[1], K1[1], I2[1], pm+J2[1], qm+K2[1])
                (J2[1]==J2[2]) && push_match!(matches, (p-1)+I1[1], (q-1)+J1[1], K1[1], pm+I2[1], J2[1], qm+K2[1])
                (K2[1]==K2[2]) && push_match!(matches, (p-1)+I1[1], (q-1)+J1[1], K1[1], pm+I2[1], qm+J2[1], K2[1])
            end
        end
    end

    # Final validation & optional splits
    if npoints(matches) < 4 || __check_edge(matches)
        return FaceMatchSet(), split_faces1, split_faces2
    end

    # Filter-increasing (uniqueness)
    if I1[1]==I1[2]; __filter_block_increasing!(matches, :j1); __filter_block_increasing!(matches, :k1); end
    if J1[1]==J1[2]; __filter_block_increasing!(matches, :i1); __filter_block_increasing!(matches, :k1); end
    if K1[1]==K1[2]; __filter_block_increasing!(matches, :i1); __filter_block_increasing!(matches, :j1); end
    if I2[1]==I2[2]; __filter_block_increasing!(matches, :j2); __filter_block_increasing!(matches, :k2); end
    if J2[1]==J2[2]; __filter_block_increasing!(matches, :i2); __filter_block_increasing!(matches, :k2); end
    if K2[1]==K2[2]; __filter_block_increasing!(matches, :i2); __filter_block_increasing!(matches, :j2); end

    if npoints(matches) < 4
        return FaceMatchSet(), split_faces1, split_faces2
    end

    # Split faces if intersection defines a proper sub-face
    main_face = create_face_from_diagonals(block1, I1[1], J1[1], K1[1], I1[2], J1[2], K1[2])
    imin, jmin, kmin = minimum(matches.i1), minimum(matches.j1), minimum(matches.k1)
    imax, jmax, kmax = maximum(matches.i1), maximum(matches.j1), maximum(matches.k1)
    if Int(imin==imax) + Int(jmin==jmax) + Int(kmin==kmax) == 1
        split_faces1 = split_face(main_face, block1, imin, jmin, kmin, imax, jmax, kmax)
        for s in split_faces1; s.BlockIndex = face1.BlockIndex; end
    end

    main_face2 = create_face_from_diagonals(block2, I2[1], J2[1], K2[1], I2[2], J2[2], K2[2])
    imin2, jmin2, kmin2 = minimum(matches.i2), minimum(matches.j2), minimum(matches.k2)
    imax2, jmax2, kmax2 = maximum(matches.i2), maximum(matches.j2), maximum(matches.k2)
    if Int(imin2==imax2) + Int(jmin2==jmax2) + Int(kmin2==kmax2) == 1
        split_faces2 = split_face(main_face2, block2, imin2, jmin2, kmin2, imax2, jmax2, kmax2)
        for s in split_faces2; s.BlockIndex = face2.BlockIndex; end
    end

    return matches, split_faces1, split_faces2
end

# -----------------------------------------------------------------------------
# find_matching_blocks — iteratively split & match outer faces
# -----------------------------------------------------------------------------
function find_matching_blocks(block1::Block,block2::Block,block1_outer::Vector{Face}, block2_outer::Vector{Face}, tol::Real=1e-6)
    block_match_indices = Vector{FaceMatchSet}()
    block1_split_faces = Face[]
    block2_split_faces = Face[]

    match = true
    while match
        match = false
        local p_idx = 0; local q_idx = 0
        for p in eachindex(block1_outer)
            block1_face = block1_outer[p]
            for q in eachindex(block2_outer)
                block2_face = block2_outer[q]
                df, split_faces1, split_faces2 = get_face_intersection(block1_face, block2_face, block1, block2; tol=tol)
                if npoints(df) > 0
                    push!(block_match_indices, df)
                    append!(block1_split_faces, split_faces1)
                    append!(block2_split_faces, split_faces2)
                    match = true; p_idx = p; q_idx = q
                    break
                end
            end
            if match; break; end
        end
        if match
            deleteat!(block1_outer, p_idx)
            deleteat!(block2_outer, q_idx)
            append!(block1_outer, block1_split_faces)
            append!(block2_outer, block2_split_faces)
            empty!(block1_split_faces)
            empty!(block2_split_faces)
        end
    end
    return block_match_indices, block1_outer, block2_outer
end

# -----------------------------------------------------------------------------
# combinations_of_nearest_blocks — neighbor pairs by centroid distance
# -----------------------------------------------------------------------------
function combinations_of_nearest_blocks(blocks::Vector{Block}; nearest_nblocks::Int=4)
    centroids = [(mean(b.X), mean(b.Y), mean(b.Z)) for b in blocks]
    distance_matrix = fill(1.0e10, length(blocks), length(blocks))
    @inbounds for i in 1:length(blocks), j in 1:length(blocks)
        if i != j
            dx = centroids[i][1]-centroids[j][1]
            dy = centroids[i][2]-centroids[j][2]
            dz = centroids[i][3]-centroids[j][3]
            distance_matrix[i,j] = sqrt(dx*dx+dy*dy+dz*dz)
        end
    end
    new_combos = Tuple{Int,Int}[]  # 0-based pairs to match the rest of the library
    for i in 1:length(blocks)
        idx = sortperm(view(distance_matrix, i, :))
        for j in idx[1:min(nearest_nblocks, length(idx))]
            if distance_matrix[i,j] < 1.0e10
                push!(new_combos, (i-1, j-1))
            end
        end
    end
    return new_combos
end

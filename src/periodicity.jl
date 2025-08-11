
# periodicity.jl
# (Only Base trig/linear algebra; no external imports)

"""
    create_rotation_matrix(axis::Symbol, θ::Real) -> Matrix{Float64}

Right-handed rotation by θ radians about axis ∈ (:x,:y,:z).
"""
function create_rotation_matrix(axis::Symbol, θ::Real)
    c = cos(θ); s = sin(θ)
    if axis === :x
        return [1 0 0; 0 c -s; 0 s c]
    elseif axis === :y
        return [c 0 s; 0 1 0; -s 0 c]
    elseif axis === :z
        return [c -s 0; s c 0; 0 0 1]
    else
        error("axis must be :x, :y, or :z")
    end
end

"""
    linear_real_transform!(X, Y, Z, A, t)

Apply affine transform: [x';y';z'] = A*[x;y;z] + t, in-place.
"""
function linear_real_transform!(X::AbstractArray, Y::AbstractArray, Z::AbstractArray,
                                A::AbstractMatrix, t::NTuple{3,Real})
    nx, ny, nk = size(X)
    tx, ty, tz = t
    @inbounds for k in 1:nk, j in 1:ny, i in 1:nx
        x = X[i,j,k]; y = Y[i,j,k]; z = Z[i,j,k]
        X[i,j,k] = A[1,1]*x + A[1,2]*y + A[1,3]*z + tx
        Y[i,j,k] = A[2,1]*x + A[2,2]*y + A[2,3]*z + ty
        Z[i,j,k] = A[3,1]*x + A[3,2]*y + A[3,3]*z + tz
    end
    return nothing
end

"""
    shift_blocks!(blocks, dx, dy, dz=0)

Translate all blocks by (dx,dy,dz). Works for 2D (K=1).
"""
function shift_blocks!(blocks::Vector{Block}, dx::Real, dy::Real, dz::Real=0)
    @inbounds for b in blocks
        b.X .+= dx; b.Y .+= dy; b.Z .+= dz
    end
    return blocks
end

"""
    translational_periodicity(blocks, direction::Symbol, L; tol=1e-8)

Check opposing boundaries along `direction` ∈ (:x,:y,:z) match under translation of length `L`.
Returns `(matches, lower_export, upper_export)`.
"""
function translational_periodicity(blocks::Vector{Block}, direction::Symbol, L::Real; tol::Real=1e-8)
    n = length(blocks)
    conn = ones(Int8, n, n); @inbounds for i in 1:n; conn[i,i] = 0; end

    lower_exp, upper_exp, lower_faces, upper_faces =
        find_bounding_faces(blocks, conn, Dict{String,Int}[]; direction=String(direction))

    blocks_shifted = deepcopy(blocks)
    if direction === :x
        shift_blocks!(blocks_shifted,  L, 0, 0)
    elseif direction === :y
        shift_blocks!(blocks_shifted,  0, L, 0)
    else
        shift_blocks!(blocks_shifted,  0, 0, L)
    end

    matches = Vector{Dict{String,Any}}()
    for lf in lower_faces
        for uf in upper_faces
            d = face_matches_to_dict(lf, uf, blocks_shifted[lf.BlockIndex+1], blocks[uf.BlockIndex+1])
            push!(matches, d)
        end
    end
    return matches, lower_exp, upper_exp
end

"""
    rotated_periodicity(blocks, axis::Symbol, θ; origin=(0,0,0))

Rotate a copy about `axis` by `θ`, then match bounding faces. Returns vector of match dicts.
"""
function rotated_periodicity(blocks::Vector{Block}, axis::Symbol, θ::Real; origin::NTuple{3,Real}=(0,0,0))
    A = create_rotation_matrix(axis, θ)
    ox, oy, oz = origin
    blocks_rot = deepcopy(blocks)
    @inbounds for b in blocks_rot
        b.X .-= ox; b.Y .-= oy; b.Z .-= oz
        linear_real_transform!(b.X, b.Y, b.Z, A, (ox, oy, oz))
    end
    n = length(blocks_rot)
    conn = ones(Int8, n, n); @inbounds for i in 1:n; conn[i,i] = 0; end
    lower_exp, upper_exp, lower_faces, upper_faces =
        find_bounding_faces(blocks_rot, conn, Dict{String,Int}[]; direction="z")
    matches = Vector{Dict{String,Any}}()
    for lf in lower_faces, uf in upper_faces
        push!(matches, face_matches_to_dict(lf, uf, blocks_rot[lf.BlockIndex+1], blocks_rot[uf.BlockIndex+1]))
    end
    return matches
end

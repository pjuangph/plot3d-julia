module Face3D

import LinearAlgebra: norm, dot

export Face, add_vertex, vertices_equals, index_equals, match_indices,
       normal, to_dict, set_block_index, set_face_id, is_edge

"""
    Face(nverts)

Face container recording vertex xyz and (i,j,k) index triplets.
Indices are stored 0-based (to match Python heritage).
"""
mutable struct Face
    # geometry
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    # index triplets (0-based)
    I::Vector{Int}
    J::Vector{Int}
    K::Vector{Int}
    # metadata
    BlockIndex::Int
    id::Int
    cx::Float64
    cy::Float64
    cz::Float64
    IMIN::Int
    JMIN::Int
    KMIN::Int
    IMAX::Int
    JMAX::Int
    KMAX::Int
    function Face(n::Int)
        new(Float64[], Float64[], Float64[],
            Int[], Int[], Int[],
            -1, 0, 0.0, 0.0, 0.0,
            typemax(Int), typemax(Int), typemax(Int),
            typemin(Int), typemin(Int), typemin(Int))
    end
end

# --- basic mutators ---
set_block_index(f::Face, idx::Int) = (f.BlockIndex = idx)
set_face_id(f::Face, id::Int)      = (f.id = id)

# --- building & stats ---
function _update_bounds!(f::Face, i::Int, j::Int, k::Int)
    f.IMIN = min(f.IMIN, i)
    f.JMIN = min(f.JMIN, j)
    f.KMIN = min(f.KMIN, k)
    f.IMAX = max(f.IMAX, i)
    f.JMAX = max(f.JMAX, j)
    f.KMAX = max(f.KMAX, k)
end

function add_vertex(f::Face, x::Real, y::Real, z::Real, i::Int, j::Int, k::Int)
    push!(f.X, float(x)); push!(f.Y, float(y)); push!(f.Z, float(z))
    push!(f.I, i);        push!(f.J, j);        push!(f.K, k)
    _update_bounds!(f, i, j, k)
    n = length(f.X)                    # update centroid incrementally
    f.cx += (x - f.cx)/n
    f.cy += (y - f.cy)/n
    f.cz += (z - f.cz)/n
    return f
end

# --- comparisons ---
index_equals(f1::Face, f2::Face) =
    length(f1.I) == length(f2.I) &&
    f1.I == f2.I && f1.J == f2.J && f1.K == f2.K

function match_indices(f1::Face, f2::Face)
    s2 = Set(zip(f2.I, f2.J, f2.K))
    [ (f1.I[t], f1.J[t], f1.K[t]) for t in eachindex(f1.I) if (f1.I[t], f1.J[t], f1.K[t]) in s2 ]
end

vertices_equals(f1::Face, f2::Face) =
    Set(zip(f1.I, f1.J, f1.K)) == Set(zip(f2.I, f2.J, f2.K))

# --- normal & edge-ness ---
function normal(f::Face, block_unused=nothing)
    if length(f.X) < 3
        return [0.0, 0.0, 0.0]
    end
    v1 = (f.X[2]-f.X[1], f.Y[2]-f.Y[1], f.Z[2]-f.Z[1])
    v2 = (f.X[3]-f.X[1], f.Y[3]-f.Y[1], f.Z[3]-f.Z[1])
    n  = (v1[2]*v2[3]-v1[3]*v2[2],
          v1[3]*v2[1]-v1[1]*v2[3],
          v1[1]*v2[2]-v1[2]*v2[1])
    nv = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
    nv == 0 ? [0.0,0.0,0.0] : [n[1]/nv, n[2]/nv, n[3]/nv]
end

# an "edge" has two constant directions among I/J/K
function is_edge(f::Face)
    m = 0
    m += (f.IMIN == f.IMAX) ? 1 : 0
    m += (f.JMIN == f.JMAX) ? 1 : 0
    m += (f.KMIN == f.KMAX) ? 1 : 0
    return m >= 2
end

# --- export as dict ---
function to_dict(f::Face)
    Dict{String,Int}(
        "IMIN"=>f.IMIN, "JMIN"=>f.JMIN, "KMIN"=>f.KMIN,
        "IMAX"=>f.IMAX, "JMAX"=>f.JMAX, "KMAX"=>f.KMAX,
        "block_index"=>f.BlockIndex, "id"=>f.id
    )
end

end # module

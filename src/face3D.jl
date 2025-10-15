module Face3D

using LinearAlgebra

# -----------------------------------------------------------------------------
# Face type
# -----------------------------------------------------------------------------
"""
Face for a structured Plot3D block.

Notes
- Indices IMIN/IMAX/JMIN/JMAX/KMIN/KMAX are **zero-based** node indices
  (consistent with how your facefunctions build them).
- `BlockIndex` is zero-based to match your usage (you often set it with `i-1`).
- `cx, cy, cz` are running centroids updated as vertices are added.
- `I, J, K` are generic integer attributes that some of your code scales
  (e.g., `l.I *= gcd_to_use`). They start at 1 by default.
"""
mutable struct Face
    # indexing/identity
    BlockIndex::Int
    id::Int

    # zero-based index bounds on each axis
    IMIN::Int; IMAX::Int
    JMIN::Int; JMAX::Int
    KMIN::Int; KMAX::Int

    # centroid (running)
    cx::Float64; cy::Float64; cz::Float64

    # freeform scalar slots (used by callers for scaling)
    I::Int; J::Int; K::Int

    # storage of vertices: (x,y,z,i,j,k) zero-based idx
    vertices::Vector{NTuple{6,Float64}}

    function Face(BlockIndex::Int, id::Int,
                  IMIN::Int, IMAX::Int, JMIN::Int, JMAX::Int, KMIN::Int, KMAX::Int,
                  cx::Float64, cy::Float64, cz::Float64, I::Int, J::Int, K::Int,
                  vertices::Vector{NTuple{6,Float64}})
        new(BlockIndex, id, IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, cx, cy, cz, I, J, K, vertices)
    end
end

# Convenience constructors
Face() = Face(-1, -1,  typemax(Int), typemin(Int), typemax(Int), typemin(Int), typemax(Int), typemin(Int),
              0.0, 0.0, 0.0, 1, 1, 1, NTuple{6,Float64}[])

# Keep compatibility with existing calls like `Face(4)`
Face(::Int) = Face()

# -----------------------------------------------------------------------------
# Basic mutators
# -----------------------------------------------------------------------------
set_block_index(f::Face, i::Int) = (f.BlockIndex = i; f)
set_face_id(f::Face, id::Int)    = (f.id = id; f)

# -----------------------------------------------------------------------------
# add_vertex
# -----------------------------------------------------------------------------
"""
    add_vertex(f::Face, x, y, z, i, j, k)

Append a vertex to the face. `(i,j,k)` must be **zero-based** node indices.
Updates:
- vertex list
- index bounds (IMIN/IMAX/...),
- running centroid `(cx,cy,cz)`.
"""
function add_vertex(f::Face, x::Real, y::Real, z::Real, i::Int, j::Int, k::Int)
    # push vertex
    push!(f.vertices, (Float64(x), Float64(y), Float64(z), Float64(i), Float64(j), Float64(k)))

    # update bounds
    if i < f.IMIN; f.IMIN = i; end
    if i > f.IMAX; f.IMAX = i; end
    if j < f.JMIN; f.JMIN = j; end
    if j > f.JMAX; f.JMAX = j; end
    if k < f.KMIN; f.KMIN = k; end
    if k > f.KMAX; f.KMAX = k; end

    # update centroid (running mean)
    n = length(f.vertices)
    if n == 1
        f.cx = x; f.cy = y; f.cz = z
    else
        α = 1.0 / n
        f.cx = (1-α)*f.cx + α*Float64(x)
        f.cy = (1-α)*f.cy + α*Float64(y)
        f.cz = (1-α)*f.cz + α*Float64(z)
    end
    return f
end

# -----------------------------------------------------------------------------
# vertices_equals (order-insensitive; compares coordinates only)
# -----------------------------------------------------------------------------
"""
    vertices_equals(a::Face, b::Face; tol=1e-8)

Return true if faces have the same set of vertex coordinates (order-insensitive).
"""
function vertices_equals(a::Face, b::Face; tol::Real=1e-8)
    la = length(a.vertices); lb = length(b.vertices)
    la == lb || return false
    used = falses(lb)
    @inbounds for va in a.vertices
        matched = false
        for (j, vb) in enumerate(b.vertices)
            used[j] && continue
            # compare xyz only
            if isapprox(va[1], vb[1]; atol=tol) &&
               isapprox(va[2], vb[2]; atol=tol) &&
               isapprox(va[3], vb[3]; atol=tol)
                used[j] = true
                matched = true
                break
            end
        end
        matched || return false
    end
    return true
end

# -----------------------------------------------------------------------------
# index_equals (exact equality of index ranges)
# -----------------------------------------------------------------------------
"""
    index_equals(a::Face, b::Face)

Return true if both faces have identical (IMIN/IMAX/JMIN/JMAX/KMIN/KMAX).
"""
index_equals(a::Face, b::Face) =
    a.IMIN==b.IMIN && a.IMAX==b.IMAX &&
    a.JMIN==b.JMIN && a.JMAX==b.JMAX &&
    a.KMIN==b.KMIN && a.KMAX==b.KMAX

# -----------------------------------------------------------------------------
# is_edge
# -----------------------------------------------------------------------------
"""
    is_edge(f::Face)

Return true if the face degenerates into an edge (i.e., one in-face span is 0).
"""
function is_edge(f::Face)
    # Determine which axis is constant for the face and check the two varying spans
    if f.IMIN == f.IMAX
        return (f.JMIN == f.JMAX) || (f.KMIN == f.KMAX)
    elseif f.JMIN == f.JMAX
        return (f.IMIN == f.IMAX) || (f.KMIN == f.KMAX)
    elseif f.KMIN == f.KMAX
        return (f.IMIN == f.IMAX) || (f.JMIN == f.JMAX)
    else
        # Not a valid face (shouldn't happen): no axis is constant
        return true
    end
end

# -----------------------------------------------------------------------------
# normal
# -----------------------------------------------------------------------------
"""
    normal(f::Face, block) -> (nx, ny, nz)

Estimate a unit normal using the four corners inferred from the face bounds on `block`.
For warped faces, two triangle normals are averaged.
"""
function normal(f::Face, block)
    # infer the 4 corners from min/max bounds (zero-based -> add 1 for array access)
    i0,i1 = f.IMIN+1, f.IMAX+1
    j0,j1 = f.JMIN+1, f.JMAX+1
    k0,k1 = f.KMIN+1, f.KMAX+1

    corners = Tuple{Float64,Float64,Float64}[]
    if f.IMIN == f.IMAX
        idxs = [(i0,j0,k0), (i0,j1,k0), (i0,j1,k1), (i0,j0,k1)]
    elseif f.JMIN == f.JMAX
        idxs = [(i0,j0,k0), (i1,j0,k0), (i1,j0,k1), (i0,j0,k1)]
    else
        idxs = [(i0,j0,k0), (i1,j0,k0), (i1,j1,k0), (i0,j1,k0)]
    end
    for (i,j,k) in idxs
        push!(corners, (block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k]))
    end

    p1,p2,p3,p4 = corners
    v12 = [p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3]]
    v13 = [p3[1]-p1[1], p3[2]-p1[2], p3[3]-p1[3]]
    v23 = [p3[1]-p2[1], p3[2]-p2[2], p3[3]-p2[3]]

    n1 = cross(v12, v13)
    n2 = cross(v23, [-v12[1], -v12[2], -v12[3]])
    n  = n1 .+ n2
    mag = norm(n)
    return mag == 0 ? (0.0,0.0,0.0) : (n[1]/mag, n[2]/mag, n[3]/mag)
end

# -----------------------------------------------------------------------------
# match_indices
# -----------------------------------------------------------------------------
"""
    match_indices(a::Face, b::Face) -> Vector{Symbol}

Return the set of **in-face axes** if both faces are parallel (same constant axis):
- For i-constant faces, returns `[:j, :k]`
- For j-constant faces, returns `[:i, :k]`
- For k-constant faces, returns `[:i, :j]`
Return `Symbol[]` if faces are not parallel.

This matches how your connectivity/search uses it (checking for length==2
to confirm faces are coplanar/parallel before comparing normals).
"""
function match_indices(a::Face, b::Face)
    if (a.IMIN==a.IMAX) && (b.IMIN==b.IMAX)
        return Symbol[:j, :k]
    elseif (a.JMIN==a.JMAX) && (b.JMIN==b.JMAX)
        return Symbol[:i, :k]
    elseif (a.KMIN==a.KMAX) && (b.KMIN==b.KMAX)
        return Symbol[:i, :j]
    else
        return Symbol[]
    end
end

# -----------------------------------------------------------------------------
# to_dict
# -----------------------------------------------------------------------------
"""
    to_dict(f::Face) -> Dict{String,Any}

Export a minimal dictionary describing the face index window, id, and block index.
"""
function to_dict(f::Face)
    return Dict{String,Any}(
        "block_index" => f.BlockIndex,
        "IMIN" => f.IMIN, "IMAX" => f.IMAX,
        "JMIN" => f.JMIN, "JMAX" => f.JMAX,
        "KMIN" => f.KMIN, "KMAX" => f.KMAX,
        "id"   => f.id,
        "cx"   => f.cx, "cy" => f.cy, "cz" => f.cz,
    )
end

end # module

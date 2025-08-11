module Utils

# Small, generic helpers used across modules.
# Keep explicit import style minimal here (no external deps).

export unique_pairs, ensure3d

"""
    unique_pairs(pairs::Vector{<:Tuple{Int,Int}}) -> Vector{Tuple{Int,Int}}

Return a list of unique unordered pairs.  
E.g. both `(3,7)` and `(7,3)` collapse to `(3,7)`, and duplicates are removed.
"""
function unique_pairs(pairs::Vector{<:Tuple{Int,Int}})
    seen = Set{Tuple{Int,Int}}()
    out  = Tuple{Int,Int}[]
    for (a,b) in pairs
        u = a <= b ? (a,b) : (b,a)
        if !(u in seen)
            push!(seen, u)
            push!(out, u)
        end
    end
    return out
end

"""
    ensure3d(A::AbstractArray) -> Array

Ensure array has 3 dimensions by adding a trailing singleton dimension if needed.
- If `A` is 2D (IMAX×JMAX), returns a view with size (IMAX, JMAX, 1).
- If `A` is 3D already, returns `A` as-is.
"""
function ensure3d(A::AbstractArray)
    nd = ndims(A)
    if nd == 3
        return A
    elseif nd == 2
        return reshape(A, size(A,1), size(A,2), 1)
    else
        error("ensure3d: expected 2D or 3D array, got ndims=$(ndims(A))")
    end
end

end # module

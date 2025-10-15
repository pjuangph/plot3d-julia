module Utils

    # -----------------------------------------------------------------------------
    # unique_pairs
    # -----------------------------------------------------------------------------
    """
        unique_pairs(n::Integer)

        Iterator of unordered index pairs `(i, j)` with `1 ≤ i < j ≤ n`.

        Example
        -------
        ```julia
        collect(unique_pairs(4))  # => [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]

        unique_pairs(n::Integer) = ((i,j) for i in 1:n for j in i+1:n)
        
        unique_pairs(ps::Vector{<:Tuple{Int,Int}}) -> Vector{Tuple{Int,Int}}
        Deduplicate an existing list of index pairs as unordered pairs.
        Ensures i < j for all returned pairs and removes duplicates.
    """
    function unique_pairs(ps::Vector{<:Tuple{Int,Int}})
        s = Set{Tuple{Int,Int}}()
        @inbounds for (i,j) in ps
            i == j && continue
            a,b = i < j ? (i,j) : (j,i)
            push!(s, (a,b))
        end
        return collect(s)
    end

    """
    -----------------------------------------------------------------------------
    ensure3d
    -----------------------------------------------------------------------------
    ensure3d(A)
    Normalize coordinate containers to a 3-component representation.
    Accepted inputs
    (X, Y, Z) where each is a 3D array (I×J×K) → returned as-is.
    A::Array{T,3} with size (I, J, 3) → returned as-is.
    A::Array{T,2} with size (N, 3) → reshaped to (1, N, 3).
    Throws an ArgumentError otherwise.
    """
    function ensure3d(A)
        if A isa Tuple && length(A) == 3
            return A
        elseif A isa AbstractArray{<:Real,3} && size(A,3) == 3
            return A
        elseif A isa AbstractArray{<:Real,2} && size(A,2) == 3
            return reshape(A, 1, size(A,1), 3)
        else
            throw(ArgumentError("ensure3d: unsupported coordinate container (got $(summary(A)))"))
        end
    end

end # module
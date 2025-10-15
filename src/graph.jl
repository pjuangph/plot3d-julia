# graph.jl — functions live directly in Plot3D; returns edge lists (no deps)

# Helpers to linearize (i,j,k) -> 1-based global vertex index
@inline _lin(i,j,k, IMAX,JMAX,KMAX) = i + (j-1)*IMAX + (k-1)*IMAX*JMAX

"""
    get_starting_vertex(blockIndex, block_sizes) -> Int

Return the **1-based** starting vertex id for a block, where `block_sizes` is
a `Vector{Tuple{Int,Int,Int}}` of (IMAX,JMAX,KMAX) for each prior block.
"""
function get_starting_vertex(blockIndex::Int, block_sizes::Vector{Tuple{Int,Int,Int}})
    blockIndex ≥ 0 || error("blockIndex must be ≥ 0 (0-based)")
    offset = 0
    for b in 1:blockIndex
        I,J,K = block_sizes[b]
        offset += I*J*K
    end
    return offset + 1
end

"""
    get_face_vertex_indices(IMIN,JMIN,KMIN, IMAX,JMAX,KMAX, block_size)
        -> Vector{Int}

Return the **1-based** linear vertex ids for the face bounded by the
zero-based ranges and the provided block size `(IMAX,JMAX,KMAX)`.
"""
function get_face_vertex_indices(IMIN::Int,JMIN::Int,KMIN::Int,
                                 IMAX::Int,JMAX::Int,KMAX::Int,
                                 block_size::Tuple{Int,Int,Int})
    I,J,K = block_size
    ids = Int[]
    if IMIN == IMAX
        i = IMIN+1
        for k in KMIN+1:KMAX+1, j in JMIN+1:JMAX+1
            push!(ids, _lin(i,j,k, I,J,K))
        end
    elseif JMIN == JMAX
        j = JMIN+1
        for k in KMIN+1:KMAX+1, i in IMIN+1:IMAX+1
            push!(ids, _lin(i,j,k, I,J,K))
        end
    else
        k = KMIN+1
        for j in JMIN+1:JMAX+1, i in IMIN+1:IMAX+1
            push!(ids, _lin(i,j,k, I,J,K))
        end
    end
    return ids
end

"""
    block_to_graph(IMAX,JMAX,KMAX; offset=0) -> Vector{Tuple{Int,Int}}

Return the undirected edge list (1-based ids) for the structured grid of size
`IMAX×JMAX×KMAX`. `offset` shifts all ids (use
`offset = get_starting_vertex(blockIndex, block_sizes)-1`).
"""
function block_to_graph(IMAX::Int, JMAX::Int, KMAX::Int; offset::Int=0)
    edges = Tuple{Int,Int}[]
    @inbounds for k in 1:KMAX, j in 1:JMAX, i in 1:IMAX
        v  = offset + _lin(i,j,k, IMAX,JMAX,KMAX)
        if i < IMAX; push!(edges, (v, offset + _lin(i+1,j,k, IMAX,JMAX,KMAX))); end
        if j < JMAX; push!(edges, (v, offset + _lin(i,j+1,k, IMAX,JMAX,KMAX))); end
        if k < KMAX; push!(edges, (v, offset + _lin(i,j,k+1, IMAX,JMAX,KMAX))); end
    end
    return edges
end

"""
    add_connectivity_to_graph(G, block_sizes, connectivities) -> edges

Given:
- `G`: any edge container supporting `push!` of `(Int,Int)` **or** `nothing`.
- `block_sizes::Vector{Tuple{Int,Int,Int}}` for each block.
- `connectivities::Vector{Dict{String,Int}}` with keys
  "block1"/"block2" subdicts (IMIN..KMAX, block_index).

Returns a **Vector{Tuple{Int,Int}}** of edges that connect matching faces
vertex-to-vertex (you can merge into your graph structure as needed).
"""
function add_connectivity_to_graph(G,
                                   block_sizes::Vector{Tuple{Int,Int,Int}},
                                   connectivities::Vector{Dict{String,Any}})
    new_edges = Tuple{Int,Int}[]
    for m in connectivities
        b1 = m["block1"]; b2 = m["block2"]
        b1_idx = b1["block_index"]; b2_idx = b2["block_index"]
        off1 = get_starting_vertex(b1_idx, block_sizes) - 1
        off2 = get_starting_vertex(b2_idx, block_sizes) - 1
        I1 = (b1["IMIN"], b1["IMAX"]); J1=(b1["JMIN"], b1["JMAX"]); K1=(b1["KMIN"], b1["KMAX"])
        I2 = (b2["IMIN"], b2["IMAX"]); J2=(b2["JMIN"], b2["JMAX"]); K2=(b2["KMIN"], b2["KMAX"])

        ids1 = get_face_vertex_indices(I1[1],J1[1],K1[1], I1[2],J1[2],K1[2], block_sizes[b1_idx+1])
        ids2 = get_face_vertex_indices(I2[1],J2[1],K2[1], I2[2],J2[2],K2[2], block_sizes[b2_idx+1])

        length(ids1) == length(ids2) || error("Mismatched face cardinality between matched faces")

        @inbounds for n in eachindex(ids1)
            e = (off1 + ids1[n], off2 + ids2[n])
            push!(new_edges, e)
            try
                G === nothing || push!(G, e)
            catch
                # ignore if G isn't a simple edge container
            end
        end
    end
    return new_edges
end


# graph.jl
# (No external imports; uses your existing connectivity code)

"""
    block_to_graph(blocks) -> Dict

Return a simple graph:
- "nodes" => Vector of block indices (0-based)
- "edges" => Vector of unique (i,j) pairs (0-based) where blocks share at least one matching face.
"""
function block_to_graph(blocks::Vector{Block})
    n = length(blocks)
    nodes = collect(0:n-1)
    edges = Tuple{Int,Int}[]

    # cache outer faces
    ofaces = [get_outer_faces(b)[1] for b in blocks]

    # scan nearest pairs to keep it cheap
    combos = combinations_of_nearest_blocks(blocks; nearest_nblocks=6)
    for (i,j) in combos
        msets, _, _ = find_matching_blocks(blocks[i+1], blocks[j+1],
                                           copy(ofaces[i+1]), copy(ofaces[j+1]))
        if !isempty(msets)
            push!(edges, (i,j))
        end
    end

    return Dict("nodes"=>nodes, "edges"=>unique_pairs(edges))
end

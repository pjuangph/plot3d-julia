module GlennHTImport

export read_ght_conn

"""
    read_ght_conn(filename) -> (max_block_index, face_matches, outer_faces, num_connections_foreach_block, nZones, zone_types, Zones, GIFs)

Julia port of `glennht/import_functions.py`. All indices converted to **0-based** like the Python.
"""
function read_ght_conn(filename::AbstractString)
    open(filename, "r") do io
        # matches
        line = readline(io)
        pairs = parse(Int, split(strip(line))[1])
        idx = 0
        blk_id1 = 0
        face_matches = Vector{Dict{String,Any}}()
        IMIN=0;JMIN=0;KMIN=0;IMAX=0;JMAX=0;KMAX=0
        for _ in 1:(pairs*2)
            line = readline(io)
            idx += 1
            toks = filter(!isempty, split(strip(line), ' '))
            nums = [parse(Int, t)-1 for t in toks]  # to 0-based
            if length(nums) == 7
                if iseven(idx)
                    push!(face_matches, Dict(
                        "block1"=>Dict("block_index"=>blk_id1, "IMIN"=>IMIN,"JMIN"=>JMIN,"KMIN"=>KMIN,"IMAX"=>IMAX,"JMAX"=>JMAX,"KMAX"=>KMAX),
                        "block2"=>Dict("block_index"=>nums[1], "IMIN"=>nums[2],"JMIN"=>nums[3],"KMIN"=>nums[4],"IMAX"=>nums[5],"JMAX"=>nums[6],"KMAX"=>nums[7])
                    ))
                else
                    blk_id1 = nums[1]
                    IMIN=nums[2]; JMIN=nums[3]; KMIN=nums[4]; IMAX=nums[5]; JMAX=nums[6]; KMAX=nums[7]
                end
            end
        end
        # outer faces
        line = readline(io)
        nFaces = parse(Int, split(strip(line))[1])
        outer_faces = Vector{Dict{String,Int}}()
        for _ in 1:nFaces
            line = readline(io)
            toks = filter(!isempty, split(strip(line), ' '))
            vals = [parse(Int, t)-1 for t in toks]
            push!(outer_faces, Dict(
                "block_index"=>vals[1],
                "IMIN"=>vals[2],"JMIN"=>vals[3],"KMIN"=>vals[4],
                "IMAX"=>vals[5],"JMAX"=>vals[6],"KMAX"=>vals[7],
                "id"=>vals[8]
            ))
        end
        # GIFs
        line = readline(io)
        nGIF = parse(Int, split(strip(line))[1])
        GIFs = Vector{Dict{String,Int}}()
        for _ in 1:nGIF
            line = readline(io)
            toks = filter(!isempty, split(strip(line), ' '))
            vals = [parse(Int, t) for t in toks]
            push!(GIFs, Dict("S1"=>vals[1], "S2"=>vals[2], "GIF_TYPE"=>vals[3], "GIF_ORDER"=>vals[4]))
        end
        # Zones
        line = readline(io)
        nZones = parse(Int, split(strip(line))[1])
        temp = split(strip(readline(io)), ' ')
        zone_types = [parse(Int, t) for t in temp if !isempty(t)]
        Zones = Int[]
        while !eof(io)
            line = readline(io)
            isempty(strip(line)) && continue
            toks = filter(!isempty, split(strip(line), ' '))
            append!(Zones, [parse(Int, t) for t in toks])
        end

        # connections per block (like Python)
        l1 = [(row["block1"]["block_index"],
               (row["block1"]["IMAX"]-row["block1"]["IMIN"])*
               (row["block1"]["JMAX"]-row["block1"]["IMIN"])*
               (row["block1"]["KMAX"]-row["block1"]["KMIN"]))
               for row in face_matches]
        l2 = [(row["block2"]["block_index"],
               (row["block2"]["IMAX"]-row["block2"]["IMIN"])*
               (row["block2"]["JMAX"]-row["block2"]["IMIN"])*
               (row["block2"]["KMAX"]-row["block2"]["KMIN"]))
               for row in face_matches]
        both = vcat(l1, l2)
        if isempty(both)
            max_block_index = -1
            num_connections_foreach_block = Array{Int}(undef, 0, 2)
        else
            max_block_index = maximum(first.(both))
            # unique by block id, keep first, then sort by first column
            seen = Dict{Int,Int}()
            for (i,(bid,_)) in enumerate(both)
                haskey(seen, bid) || (seen[bid] = i)
            end
            rows = collect(values(seen))
            arr = [both[i] for i in rows]
            sort!(arr; by=x->x[1])
            num_connections_foreach_block = reshape(collect(Iterators.flatten(arr)), :, 2)
        end

        return max_block_index, face_matches, outer_faces, num_connections_foreach_block, nZones, zone_types, Zones, GIFs
    end
end

end # module

module block_search

using Parameters, StaticArrays
export read_plot3D, identify_faces, match_faces

struct Block
    IMAX::Int64
    JMAX::Int64
    KMAX::Int64
    X::Array{Float64,3}
    Y::Array{Float64,3}
    Z::Array{Float64,3}
end

function read_plot3D(file::String)
    block_dims = Tuple[]
    reader = String[]
    NBLOCKS = 0
    if isfile(file)
        reader = readlines(file)
        NBLOCKS = parse(Int64, reader[1])
        #println("NBLOCKS: ", NBLOCKS)
        starting_line_num = 2 + NBLOCKS
        #println("starting_line_num: ", starting_line_num)
        for a in 1:NBLOCKS
            block_size = Tuple(parse.(Int64, split(reader[a + 1])))
            push!(block_dims, block_size)
        end
    else
        println("File not found!")
    end

    blocks = Vector{Block}(undef, 0)

    for c in 1:NBLOCKS
        X_Flat = Float64[]; Y_Flat = Float64[]; Z_Flat = Float64[]
        (IMAX, JMAX, KMAX) = block_dims[c]
        Max_Value = max(IMAX, JMAX, KMAX)

        for i in 1:JMAX
            xblock_rows = parse.(Float64, split(reader[i + starting_line_num - 1])) 
            for a in 1:IMAX
                try
                    append!(X_Flat, xblock_rows[a])
                catch e
                    if e isa BoundsError
                    else
                        rethrow(e)
                    end
                end
            end
        end
        X = reshape(X_Flat, IMAX, JMAX, 1)

        for j in 1:JMAX
            yblock_rows = parse.(Float64, split(reader[j + JMAX + starting_line_num - 1]))
            for a in 1:IMAX
                try
                    append!(Y_Flat, yblock_rows[a])
                catch e
                    if e isa BoundsError
                    else
                        rethrow(e)
                    end
                end           
            end
        end
        Y = reshape(Y_Flat, IMAX, JMAX, 1)

        for k in 1:JMAX
            zblock_rows = parse.(Float64, split(reader[k + (2*JMAX) + starting_line_num - 1]))
            for a in 1:IMAX
                try
                    append!(Z_Flat, zblock_rows[a])
                catch e
                    if e isa BoundsError
                    else
                        rethrow(e)
                    end
                end
            end
        end
        Z = reshape(Z_Flat, IMAX, JMAX, 1)

        starting_line_num += JMAX * 3
        b = Block(IMAX, JMAX, KMAX, X, Y, Z)
        push!(blocks, b)
    end
    return blocks
end

function identify_faces(blocks)
    faces = []
    for block in blocks
        IMAX, JMAX, KMAX = block.IMAX, block.JMAX, block.KMAX
        X, Y, Z = block.X, block.Y, block.Z

        # Extract faces using X[begin, :, end] format
        block_faces = [
            [(X[begin, j, 1], Y[begin, j, 1], Z[begin, j, 1]) for j in 1:JMAX],  # k = 1 (bottom face) = klo
            [(X[end, j, KMAX], Y[end, j, KMAX], Z[end, j, KMAX]) for j in 1:JMAX],  # k = KMAX (top face) = khi
            [(X[i, begin, k], Y[i, begin, k], Z[i, begin, k]) for i in 1:IMAX, k in 1:KMAX],  # j = begin (front face) = jlo
            [(X[i, end, k], Y[i, end, k], Z[i, end, k]) for i in 1:IMAX, k in 1:KMAX],  # j = end (back face) = jhi
            [(X[begin, j, k], Y[begin, j, k], Z[begin, j, k]) for j in 1:JMAX, k in 1:KMAX],  # i = begin (left face) = ilo
            [(X[end, j, k], Y[end, j, k], Z[end, j, k]) for j in 1:JMAX, k in 1:KMAX]  # i = end (right face) = ihi
        ]
        # Convert each face to a set of unique vertices
        block_faces = [collect(Set(reshape(face, :))) for face in block_faces]
        push!(faces, block_faces)
    end
    return faces
end

function match_faces(faces)
    face_labels = [:klo, :khi, :jlo, :jhi, :ilo, :ihi]  # Labels for the faces
    block_connections = []  # Store connections as tuples

    for (i, block_faces) in enumerate(faces)
        #println("Checking faces of block $i:")
        for (face_index, face) in enumerate(block_faces)
            matched = false
            for (j, other_block_faces) in enumerate(faces)
                if i != j  # Avoid comparing the block with itself
                    for (other_face_index, other_face) in enumerate(other_block_faces)
                        if Set(face) == Set(other_face)  # Check if the faces match
                            #println("    Matches with face $other_face_index of block $j")
                            matched = true
                            # Add the connection to the list
                            push!(block_connections, 
                                ((Symbol("blk$i"), face_labels[face_index]) => (Symbol("blk$j"), face_labels[other_face_index]))
                            )
                            break  # Stop checking once a match is found
                        end
                    end
                    if matched
                        break
                    end
                end
            end
            if !matched
                #println("    No match found for face $face_index of block $i")
            end
        end
    end

    # Convert the list of connections to a tuple
    return tuple(block_connections...)
end

# Example usage
blocks = read_plot3D("/home/maiaross/JuliaGMSH/Plot3D/test/Copied_C_grid.p3d")
faces = identify_faces(blocks)
block_connections = match_faces(faces)
#println("Block Connections:")
println(block_connections)

end

# split_block.jl
# (No external imports required; all Base)

# Simple Direction enum; avoid K in 2D (KMAX==1)
@enum Direction::UInt8 begin
    DirI = 0
    DirJ = 1
    DirK = 2
end

# Internal step search (ported from Python)
function _step_search(total_cells::Int, gcd_cells::Int, ncells_per_block::Int, denominator::Int; forward::Bool)
    initial_guess = max(1, Int(cld(ncells_per_block, denominator)))
    step_size = initial_guess
    rem_cells = total_cells % (step_size * denominator)
    ijkmax_rem = rem_cells ÷ denominator
    inc = forward ? 1 : -1

    while ((step_size % gcd_cells != 0) ||
           ((ijkmax_rem - 1) % gcd_cells != 0)) &&
          step_size > initial_guess ÷ 2 &&
          step_size < Int(round(initial_guess * 1.5))
        if (step_size % gcd_cells == 0) && ((ijkmax_rem - 1) % gcd_cells == 0)
            break
        end
        step_size += inc
        rem_cells = total_cells % (step_size * denominator)
        ijkmax_rem = rem_cells ÷ denominator
    end

    if (step_size % gcd_cells != 0) || ((ijkmax_rem - 1) % gcd_cells != 0)
        return -1
    end
    return step_size
end

"""
    split_blocks(blocks, ncells_per_block; prefer::Union{Nothing,Direction}=nothing)

Split each block into sub-blocks targeting about `ncells_per_block` cells, along the
dominant dimension unless `prefer` is specified. Honors the grid GCD so faces line up.

2D blocks are supported (KMAX==1): DirK is never chosen for those blocks.
"""
function split_blocks(blocks::Vector{Block}, ncells_per_block::Int; prefer::Union{Nothing,Direction}=nothing)
    new_blocks = Block[]
    for block in blocks
        IMAX, JMAX, KMAX = block.IMAX, block.JMAX, block.KMAX
        ci = max(1, IMAX - 1); cj = max(1, JMAX - 1); ck = max(1, KMAX - 1)
        total_cells = ci * cj * ck
        if total_cells <= ncells_per_block
            push!(new_blocks, block); continue
        end

        # choose split direction
        lengths = [(DirI, ci), (DirJ, cj)]
        (KMAX > 1) && push!(lengths, (DirK, ck))
        direction_to_use = isnothing(prefer) ? (reduce((a,b)->a[2]≥b[2] ? a : b, lengths)[1]) : prefer

        gcd_cells = gcd(IMAX-1, gcd(JMAX-1, max(1, KMAX-1)))

        if direction_to_use == DirI
            denom = JMAX * KMAX
            ss = _step_search(total_cells, gcd_cells, ncells_per_block, denom; forward=false)
            ss == -1 && (ss = _step_search(total_cells, gcd_cells, ncells_per_block, denom; forward=true))
            ss == -1 && error("No valid step size found for I-splits.")

            iprev = 1
            i = ss
            while i < IMAX
                X = @view block.X[iprev:i+1, :, :]
                Y = @view block.Y[iprev:i+1, :, :]
                Z = @view block.Z[iprev:i+1, :, :]
                push!(new_blocks, Block(X, Y, Z))
                iprev = i + 1
                i += ss
            end
            if iprev < IMAX
                X = @view block.X[iprev:end, :, :]
                Y = @view block.Y[iprev:end, :, :]
                Z = @view block.Z[iprev:end, :, :]
                push!(new_blocks, Block(X, Y, Z))
            end

        elseif direction_to_use == DirJ
            denom = IMAX * KMAX
            ss = _step_search(total_cells, gcd_cells, ncells_per_block, denom; forward=false)
            ss == -1 && (ss = _step_search(total_cells, gcd_cells, ncells_per_block, denom; forward=true))
            ss == -1 && error("No valid step size found for J-splits.")

            jprev = 1
            j = ss
            while j < JMAX
                X = @view block.X[:, jprev:j+1, :]
                Y = @view block.Y[:, jprev:j+1, :]
                Z = @view block.Z[:, jprev:j+1, :]
                push!(new_blocks, Block(X, Y, Z))
                jprev = j + 1
                j += ss
            end
            if jprev < JMAX
                X = @view block.X[:, jprev:end, :]
                Y = @view block.Y[:, jprev:end, :]
                Z = @view block.Z[:, jprev:end, :]
                push!(new_blocks, Block(X, Y, Z))
            end

        else # DirK
            @assert KMAX > 1 "DirK chosen but this is a 2D block (KMAX==1)."
            denom = IMAX * JMAX
            ss = _step_search(total_cells, gcd_cells, ncells_per_block, denom; forward=false)
            ss == -1 && (ss = _step_search(total_cells, gcd_cells, ncells_per_block, denom; forward=true))
            ss == -1 && error("No valid step size found for K-splits.")

            kprev = 1
            k = ss
            while k < KMAX
                X = @view block.X[:, :, kprev:k+1]
                Y = @view block.Y[:, :, kprev:k+1]
                Z = @view block.Z[:, :, kprev:k+1]
                push!(new_blocks, Block(X, Y, Z))
                kprev = k + 1
                k += ss
            end
            if kprev < KMAX
                X = @view block.X[:, :, kprev:end]
                Y = @view block.Y[:, :, kprev:end]
                Z = @view block.Z[:, :, kprev:end]
                push!(new_blocks, Block(X, Y, Z))
            end
        end
    end
    return new_blocks
end

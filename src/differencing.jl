
# differencing.jl
# (No external imports needed; uses Base math)

"""
    ddx(F, X) -> Array

Centered difference in i-direction for scalar field `F` on (i,j,k).
One-sided at i=1 and i=end. Works when KMAX==1 (2D grids).
"""
function ddx(F::AbstractArray, X::AbstractArray)
    @assert size(F) == size(X)
    ni, nj, nk = size(F)
    G = similar(F)
    @inbounds for k in 1:nk, j in 1:nj
        # forward difference at i=1
        G[1, j, k] = (F[2, j, k] - F[1, j, k]) / (X[2, j, k] - X[1, j, k])
        # centered interior
        for i in 2:ni-1
            dx = X[i+1, j, k] - X[i-1, j, k]
            G[i, j, k] = (F[i+1, j, k] - F[i-1, j, k]) / dx
        end
        # backward difference at i=end
        G[ni, j, k] = (F[ni, j, k] - F[ni-1, j, k]) / (X[ni, j, k] - X[ni-1, j, k])
    end
    return G
end

"""
    ddy(F, Y) -> Array

Centered difference in j-direction; one-sided at ends. Works when KMAX==1.
"""
function ddy(F::AbstractArray, Y::AbstractArray)
    @assert size(F) == size(Y)
    ni, nj, nk = size(F)
    G = similar(F)
    @inbounds for k in 1:nk, i in 1:ni
        G[i, 1, k] = (F[i, 2, k] - F[i, 1, k]) / (Y[i, 2, k] - Y[i, 1, k])
        for j in 2:nj-1
            dy = Y[i, j+1, k] - Y[i, j-1, k]
            G[i, j, k] = (F[i, j+1, k] - F[i, j-1, k]) / dy
        end
        G[i, nj, k] = (F[i, nj, k] - F[i, nj-1, k]) / (Y[i, nj, k] - Y[i, nj-1, k])
    end
    return G
end

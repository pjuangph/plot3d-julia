# differencing.jl — functions live directly in Plot3D

# Return types:
# - find_edges: Vector of NamedTuples with fields (i,j,k, di, dj, dk)
#   where di, dj, dk are tuples: ((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f))
# - find_face_edges: same idea for 2D face arrays (p,q indices; dp,dq)

"""
    find_edges(X, Y, Z)

Compute forward/backward edge vectors along each axis for a 3D block.
Returns a vector of NamedTuples:
    (i,j,k, di=((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)),
            dj=((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)),
            dk=((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)))
"""
function find_edges(X::AbstractArray{<:Real,3},
                    Y::AbstractArray{<:Real,3},
                    Z::AbstractArray{<:Real,3})
    size(X) == size(Y) == size(Z) || error("X,Y,Z must have identical sizes")
    I,J,K = size(X)
    out = NamedTuple[]
    @inbounds for k in 1:K, j in 1:J, i in 1:I
        # neighbors with clamped indices for backward/forward differences
        ib = max(i-1, 1); ifw = min(i+1, I)
        jb = max(j-1, 1); jfw = min(j+1, J)
        kb = max(k-1, 1); kfw = min(k+1, K)

        di_b = (X[i,j,k]-X[ib,j,k], Y[i,j,k]-Y[ib,j,k], Z[i,j,k]-Z[ib,j,k])
        di_f = (X[ifw,j,k]-X[i,j,k], Y[ifw,j,k]-Y[i,j,k], Z[ifw,j,k]-Z[i,j,k])
        dj_b = (X[i,j,k]-X[i,jb,k], Y[i,j,k]-Y[i,jb,k], Z[i,j,k]-Z[i,jb,k])
        dj_f = (X[i,jfw,k]-X[i,j,k], Y[i,jfw,k]-Y[i,j,k], Z[i,jfw,k]-Z[i,j,k])
        dk_b = (X[i,j,k]-X[i,j,kb], Y[i,j,k]-Y[i,j,kb], Z[i,j,k]-Z[i,j,kb])
        dk_f = (X[i,j,kfw]-X[i,j,k], Y[i,j,kfw]-Y[i,j,k], Z[i,j,kfw]-Z[i,j,k])

        push!(out, (i=i-1, j=j-1, k=k-1,
                    di=(di_b, di_f), dj=(dj_b, dj_f), dk=(dk_b, dk_f)))
    end
    return out
end

"""
    find_face_edges(X, Y, Z)

2D version for a face array (PMAX×QMAX). Returns a vector of NamedTuples:
    (p,q, dp=((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)),
          dq=((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)))
"""
function find_face_edges(X::AbstractArray{<:Real,2},
                         Y::AbstractArray{<:Real,2},
                         Z::AbstractArray{<:Real,2})
    size(X) == size(Y) == size(Z) || error("X,Y,Z must have identical sizes")
    P,Q = size(X)
    out = NamedTuple[]
    @inbounds for q in 1:Q, p in 1:P
        pb = max(p-1, 1); pf = min(p+1, P)
        qb = max(q-1, 1); qf = min(q+1, Q)

        dp_b = (X[p,q]-X[pb,q], Y[p,q]-Y[pb,q], Z[p,q]-Z[pb,q])
        dp_f = (X[pf,q]-X[p,q], Y[pf,q]-Y[p,q], Z[pf,q]-Z[p,q])
        dq_b = (X[p,q]-X[p,qb], Y[p,q]-Y[p,qb], Z[p,q]-Z[p,qb])
        dq_f = (X[p,qf]-X[p,q], Y[p,qf]-Y[p,q], Z[p,qf]-Z[p,q])

        push!(out, (p=p-1, q=q-1, dp=(dp_b, dp_f), dq=(dq_b, dq_f)))
    end
    return out
end

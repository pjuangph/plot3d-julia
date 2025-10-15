using Test
using Downloads
using LinearAlgebra
using Plot3D

# ---- helpers (local to the test) ---------------------------------------------

"Build a 3×3 rotation matrix about axis ∈ ('x','y','z') by angle θ (radians)."
function _rotation_matrix(axis::AbstractString, θ::Real)
    c = cos(θ); s = sin(θ)
    if axis == "x"
        return SMatrix{3,3,Float64}((1,0,0,  0,c,-s,  0,s,c))
    elseif axis == "y"
        return SMatrix{3,3,Float64}((c,0,s,  0,1,0,  -s,0,c))
    elseif axis == "z"
        return SMatrix{3,3,Float64}((c,-s,0,  s,c,0,  0,0,1))
    else
        error("axis must be 'x','y', or 'z'")
    end
end

"Rotate a Block’s coordinates in-place by 3×3 matrix R."
function _rotate_block!(B::Plot3D.Block, R::AbstractMatrix{<:Real})
    @inbounds for k in 1:B.KMAX, j in 1:B.JMAX, i in 1:B.IMAX
        x = B.X[i,j,k]; y = B.Y[i,j,k]; z = B.Z[i,j,k]
        rx = R[1,1]*x + R[1,2]*y + R[1,3]*z
        ry = R[2,1]*x + R[2,2]*y + R[2,3]*z
        rz = R[3,1]*x + R[3,2]*y + R[3,3]*z
        B.X[i,j,k] = rx; B.Y[i,j,k] = ry; B.Z[i,j,k] = rz
    end
    return B
end

"Deep-copy a Block and return a rotated copy (non-mutating convenience)."
function _rotated_copy(B::Plot3D.Block, R::AbstractMatrix{<:Real})
    C = Plot3D.Block(copy(B.X), copy(B.Y), copy(B.Z))
    return _rotate_block!(C, R)
end

# ---- test --------------------------------------------------------------------

@testset "Axial duplication: rotate, connectivity, rotated periodicity" begin
    mktempdir() do tmp
        cd(tmp) do
            # --- download original mesh (ASCII) ---
            url  = "https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz"
            path = "VSPT_ASCII.xyz"
            if !isfile(path)
                @info "Downloading $url ..."
                Downloads.download(url, path)
            end
            @test isfile(path) && filesize(path) > 0

            # --- read base blocks (ASCII) ---
            base_blocks = read_plot3D_ascii(path)
            @test !isempty(base_blocks)

            # --- geometry/rotation settings (mirror notebook) ---
            number_of_blades = 55
            rotation_angle_deg = 360.0 / number_of_blades
            θ = deg2rad(rotation_angle_deg)
            copies = 3
            R = _rotation_matrix("x", θ)

            # --- build rotated set (original + 2 rotated copies) ---
            rotated_blocks = Plot3D.Block[]
            append!(rotated_blocks, base_blocks)
            for c in 2:copies
                Rc = _rotation_matrix("x", (c-1)*θ)
                for b in base_blocks
                    push!(rotated_blocks, _rotated_copy(b, Rc))
                end
            end
            @test length(rotated_blocks) == copies * length(base_blocks)

            # --- connectivity on ORIGINAL (unrotated) blocks ---
            face_matches, outer_faces = find_matching_blocks(base_blocks; tol=1e-8)
            @test isa(face_matches, Vector)
            @test isa(outer_faces, Vector)

            # --- rotated periodicity using original blocks/outer faces ---
            # API: rotated_periodicity(blocks, face_matches, outer_faces; rotation_angle, rotation_axis)
            periodic_faces, outer_keep, lower_faces, upper_faces =
                rotated_periodicity(base_blocks, face_matches, outer_faces;
                                    rotation_angle = rotation_angle_deg,
                                    rotation_axis = "x")

            @test isa(periodic_faces, Vector)
            @test isa(outer_keep, Vector)
            @test isa(lower_faces, Vector)
            @test isa(upper_faces, Vector)

            # structural checks for a pair dict
            function _ok_face(d::Dict)
                req = ("block_index","IMIN","JMIN","KMIN","IMAX","JMAX","KMAX")
                all(haskey(d,k) for k in req)
            end
            function _ok_pair(d::Dict)
                haskey(d,"block1") && haskey(d,"block2") && _ok_face(d["block1"]) && _ok_face(d["block2"])
            end
            @test all(_ok_pair, periodic_faces)

            @info "rotated periodic matches: $(length(periodic_faces)) outer_keep: $(length(outer_keep))"

            # --- inner periodicities (between duplicated rings) ---
            # next ring (c=2) faces map to previous ring (c=1), etc.
            inner_period = Dict{String,Any}[]
            nb = length(base_blocks)
            for ring in 2:copies
                offset_left  = (ring-1)*nb
                offset_right = (ring-2)*nb
                for p in periodic_faces
                    q = deepcopy(p)
                    q["block1"]["block_index"] += offset_left
                    q["block2"]["block_index"] += offset_right
                    push!(inner_period, q)
                end
            end
            @test !isempty(inner_period)
            @test all(_ok_pair, inner_period)

            # --- outer periodicities (wrap last ring to first ring) ---
            outer_period = deepcopy(periodic_faces)
            for p in outer_period
                p["block2"]["block_index"] += (copies-1)*nb
            end
            @test all(_ok_pair, outer_period)

            # very light sanity: no OOB indices vs rotated_blocks
            maxidx = length(rotated_blocks) - 1
            function _inbounds_pair(d)
                0 ≤ d["block1"]["block_index"] ≤ maxidx &&
                0 ≤ d["block2"]["block_index"] ≤ maxidx
            end
            @test all(_inbounds_pair, vcat(inner_period, outer_period))
        end
    end
end

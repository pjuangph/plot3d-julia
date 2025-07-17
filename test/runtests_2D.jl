using Test
using Plot3D
using StaticArrays

include("../src/Reformat_Match.jl")
using.Faces


@testset "Faces match" begin #CORRECT :)
    f1 = Face2D([SVector(0.0, 0.0), SVector(1.0, 0.0)])
    f2 = Face2D([SVector(0.0, 0.0), SVector(1.0, 0.0)])
    forward_match_check = forward_match(f1, f2)
        @test forward_match_check == true
    end

@testset "Faces match reverse" begin #CORRECT :)
    f1 = Face2D([SVector(0.0, 0.0), SVector(1.0, 0.0)])
    f2 = Face2D([SVector(1.0, 0.0), SVector(0.0, 0.0)])
    reverse_match_check = reverse_match(f1, f2)
        @test reverse_match_check == true
    end

@testset "Faces do not match" begin #CORRECT :)
    f1 = Face2D([SVector(0.0, 0.0), SVector(1.0, 0.0)])
    f2 = Face2D([SVector(2.0, 0.0), SVector(2.0, 1.0)])
    forward_match_check = forward_match(f1, f2)
    reverse_match_check = reverse_match(f1, f2)
        @test final_match_check(f1, f2) == false
    end 
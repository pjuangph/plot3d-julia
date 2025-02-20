using Plot3D
using Test

@testitem "test_read.jl" begin
    include("test_read.jl")
    @test read_VSPT() == "completed"
end

@testitem "test_write.jl" begin
    include("test_write.jl")
    @test WritePlot3DMultiBlock() == "completed"
end

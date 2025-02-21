using Plot3D
using Test
import Downloads
using Plot3D

@testitem "test_read.jl" begin
    mkpath("test_files")
    if !isfile("test_files/VSPT_ASCII.xyz")
        Downloads.download(raw"https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz","test_files/VSPT_ASCII.xyz")
    end
    blocks = ReadMultiBlock("test_files/VSPT_ASCII.xyz";binary=false)   
    if (size(blocks[1].X)[1] == 257 && size(blocks[1].X)[2] == 101 && size(blocks[1].X)[3] == 33) && 
        (size(blocks[2].X)[1] == 269 && size(blocks[2].X)[2] == 101 && size(blocks[2].X)[3] == 53)
        @test true
    else
        @test false
    end
end

@testitem "test_write.jl" begin
    mkpath("test_files")
    if !isfile("test_files/VSPT_ASCII.xyz")
        Downloads.download(raw"https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz","test_files/VSPT_ASCII.xyz")
    end
    blocks = ReadMultiBlock("test_files/VSPT_ASCII.xyz",binary=false)
    WriteMultiBlock("test_files/VSPT_Binary.xyz",blocks,binary=true,double_precision=false)
    blocks_written = ReadMultiBlock("test_files/VSPT_Binary.xyz",binary=true,double_precision=false)
    if (size(blocks_written[1].X)[1] == size(blocks[1].X)[1] && size(blocks_written[1].X)[2] == size(blocks[1].X)[2] && size(blocks_written[1].X)[3] == size(blocks[1].X)[3]) && 
        (size(blocks_written[2].X)[1] == size(blocks[2].X)[1] && size(blocks_written[2].X)[2] == size(blocks[2].X)[2] && size(blocks_written[2].X)[3] == size(blocks[2].X)[3])
        @test true
    else
        @test false
    end    
end

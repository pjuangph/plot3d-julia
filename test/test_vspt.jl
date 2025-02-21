import Downloads
using Plot3D

function read_VSPT()
    mkpath("test_files")
    if !isfile("test_files/VSPT_ASCII.xyz")
        Downloads.download(raw"https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz","test_files/VSPT_ASCII.xyz")
    end
    blocks = ReadMultiBlock("test_files/VSPT_ASCII.xyz";binary=false)   
    if (size(blocks[1].X)[1] == 257 && size(blocks[1].X)[2] == 101 && size(blocks[1].X)[3] == 33) && 
        (size(blocks[2].X)[1] == 269 && size(blocks[2].X)[2] == 101 && size(blocks[2].X)[3] == 53)
        return "completed"
    else
        return "fail"
    end
end



function write_VSPT()
    mkpath("test_files")
    if !isfile("test_files/VSPT_ASCII.xyz")
        Downloads.download(raw"https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz","test_files/VSPT_ASCII.xyz")
    end
    blocks = ReadMultiBlock("test_files/VSPT_ASCII.xyz";binary=false)
    WriteMultiBlock("test_files/VSPT_Binary.xyz",blocks;binary=true)
    blocks_written = ReadMultiBlock("test_files/VSPT_Binary.xyz",binary=true)
    if (size(blocks_written[1].X)[1] == size(blocks[1].X)[1] && size(blocks_written[1].X)[2] == size(blocks[1].X)[2] && size(blocks_written[1].X)[3] == size(blocks[1].X)[3]) && 
        (size(blocks_written[2].X)[1] == size(blocks[2].X)[1] && size(blocks_written[2].X)[2] == size(blocks[2].X)[2] && size(blocks_written[2].X)[3] == size(blocks[2].X)[3])
        return "completed"
    else
        return "fail"
    end
end
using Test
using Plot3D
using StaticArrays


@testset "Faces match" begin 
    f1 = FaceND([SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)])
    f2 = FaceND([SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)])
    forward_match_check = forward_match(f1, f2)
        @test forward_match_check == true
    end
@testset "Faces match reverse" begin
    f1 = FaceND([SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)])
    f2 = FaceND([SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)])
    reverse_match_check = reverse_match(f1, f2)
        @test reverse_match_check == true
    end
@testset "Faces do not match" begin
    f1 = FaceND([SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)])
    f2 = FaceND([SVector(2.0, 0.0, 0.0), SVector(2.0, 1.0, 0.0)])
    forward_match_check = forward_match(f1, f2)
    reverse_match_check = reverse_match(f1, f2)
        @test final_match_check(f1, f2) == false
    end

@testset "Block3DXYZ struct" begin
    X = reshape(collect(1.0:8.0), 2, 2, 2)
    Y = reshape(collect(11.0:18.0), 2, 2, 2)
    Z = reshape(collect(21.0:28.0), 2, 2, 2)
    b = Block3DXYZ(X, Y, Z)
    @test b.x == X
    @test b.y == Y
    @test b.z == Z
end

@testset "read_xyz_block" begin
    # Fake lines for a 2x2x2 block
    lines = [
        "1.0 2.0", "3.0 4.0",  # x, k=1, j=1:2
        "5.0 6.0", "7.0 8.0",  # x, k=2, j=1:2
        "11.0 12.0", "13.0 14.0",
        "15.0 16.0", "17.0 18.0",
        "21.0 22.0", "23.0 24.0",
        "25.0 26.0", "27.0 28.0"
    ]
    dims = (2, 2, 2)
    block, nextidx = read_xyz_block(1, dims, lines)
    @test block.x[1,1,1] == 1.0
    @test block.x[2,2,2] == 8.0
    @test block.y[2,2,2] == 18.0
    @test block.z[2,2,2] == 28.0
    @test nextidx == 13
end

@testset "read_structured_xyz" begin
    # Write a temporary .xyz file for testing
    fname = tempname()
    open(fname, "w") do io
        write(io, "1\n2 2 2\n")
        # x
        write(io, "1.0 2.0\n3.0 4.0\n5.0 6.0\n7.0 8.0\n")
        # y
        write(io, "11.0 12.0\n13.0 14.0\n15.0 16.0\n17.0 18.0\n")
        # z
        write(io, "21.0 22.0\n23.0 24.0\n25.0 26.0\n27.0 28.0\n")
    end
    blocks = read_structured_xyz(fname)
    @test length(blocks) == 1
    b = blocks[1]
    @test b.x[1,1,1] == 1.0
    @test b.x[2,2,2] == 8.0
    @test b.y[1,1,1] == 11.0
    @test b.z[2,2,2] == 28.0
end

@testset "Block Connections" begin

    X1 = [0.0 1.0; 0.0 1.0]
    Y1 = [0.0 0.0; 1.0 1.0]
    X2 = [1.0 2.0; 1.0 2.0]
    Y2 = [0.0 0.0; 1.0 1.0]

    # Create Block2D objects
    b1 = Plot3D.BlockType.Block2D(X1, Y1)
    b2 = Plot3D.BlockType.Block2D(X2, Y2)
    blocks = [b1, b2]

    # Generate faces and connections
    block_faces = [Plot3D.Faces.Block2DFaces(b) for b in blocks]
    connections = print_connections(block_faces)

    # Define expected connections
    expected = [
        ((:blk1, :jhi) => (:blk2, :jlo))
    ]

    # Test that all expected connections are present
    for conn in expected
        @test conn in connections
    end

    # Optionally, test that there are no extra connections
    @test length(connections) == length(expected)
end


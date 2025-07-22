# Plot3D.jl

A Julia package designed to take .p3d and .xyz files modeling multi-block meshes and output the 
block connections in the form ((:blk1, :ihi) => (:blk2, :ilo), (:blk1, :ilo) => (:blk2, :ihi)).

## Installation and Usage

```
git clone https://github.com/ReadingRocks2973/plot3d-julia.git
```

##### From the Plot3D folder, for just the multiblock connections: 

```
julia ./test/Run_File3D.jl ./path/to/.xyz/or/.p3d/file
```

### API Documentation

```@autodocs
Modules = [Plot3D]
```

```@docs
Plot3D.Faces.reverse_match
Plot3D.Faces.compare_blocks
Plot3D.xyzreader.read_structured_xyz
Plot3D.BlockType.read_block
Plot3D.Faces.FaceND
Plot3D.Faces.Block2DFaces
Plot3D.Faces.Block3DFaces
Plot3D.Faces.final_match_check
Plot3D.Faces.print_connections
Plot3D.Faces.forward_match
Plot3D.BlockType.read_blocks
Plot3D.Faces.block_xyz_faces
Plot3D.xyzreader.read_xyz_block 
Plot3D.BlockType.Block2D
Plot3D.BlockType.Block3D
Plot3D.BlockType.read_2d_blocks 
Plot3D.BlockType.read_3d_blocks 

```
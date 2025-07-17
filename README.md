#Plot3D.jl

A Julia package designed to take .p3d and .xyz files modeling multi-block meshes and output the 
block connections in the form ((:blk1, :ihi) => (:blk2, :ilo), (:blk1, :ilo) => (:blk2, :ihi)). This
can then be input into Cygnus.jl to simulate complex multiphysics problems such as fluid dynamics
and heat transfer. 

## Installation

FILL IN HERE - git clone? or registry like Cygnus?

## Usage

From Plot3D folder:

julia ./test/Run_File3D.jl ./data/path_to_.xyz_or.p3d_file
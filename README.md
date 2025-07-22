# Plot3D.jl

A Julia package designed to take .p3d and .xyz files modeling multi-block meshes and output the 
block connections in the form ((:blk1, :ihi) => (:blk2, :ilo), (:blk1, :ilo) => (:blk2, :ihi)).
## Installation

'''
git clone https://github.com/ReadingRocks2973/plot3d-julia.git
'''

## Usage

From Plot3D folder:

julia ./test/Run_File3D.jl ./data/path_to_.xyz_or.p3d_file
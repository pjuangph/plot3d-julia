module XYZReader

import ..ReadPlot3D: read_plot3D_ascii
import ..Block3D: Block

export read_xyz_ascii

"""
    read_xyz_ascii(path) -> Vector{Block}

Thin alias for `read_plot3D_ascii(path)`. Keeps API parity with prior code.
"""
read_xyz_ascii(path::AbstractString) = read_plot3D_ascii(path)

end # module

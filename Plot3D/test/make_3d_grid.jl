using CurvilinearGrids
using Printf

#     |---------|
#     |    m2   |
#     |---------|--------|
#     |    m1   |   m3   |
#     |---------|--------|

m1 = rectlinear_grid(
  0:0.1:5,   # x
  0:0.2:2,  # y
  0:0.25:5, # z
  :meg6
);

m2 = rectlinear_grid(
  0:0.1:5,   # x
  2:0.2:7, # y
  0:0.25:5, # z
  :meg6
);

m3 = rectlinear_grid(
  5:0.1:8,   # x
  0:0.2:2,  # y
  0:0.25:5, # z
  :meg6
);

function write_plot3d(blocks)
  max_col_width = 10
  dims = ("x", "y", "z")
  open("mesh.xyz", "w") do file
    write(file, "$(length(blocks))\n")

    for (blk_id, block) in enumerate(blocks)
      blk_size = join(block.nnodes, " ")
      write(file, "$blk_size\n")
    end

    for (blk_id, block) in enumerate(blocks)
      @info "Block $blk_id"
      node_coordinates = coords(block)
      for coordinate_dim in node_coordinates # x, y, or z

        ni, nj, nk = size(coordinate_dim)
        for k in 1:nk
          for j in 1:nj
            for i in 1:ni
              # d = @sprintf("%.4e", coordinate_dim[i, j, k])
              # write(file, "$d ")
              write(file, "$(coordinate_dim[i, j, k]) ")
            end
            write(file, "\n")
          end
        end
      end
    end
  end

  return nothing
end

write_plot3d([m1, m2, m3])

save_vtk(m1, "m1")
save_vtk(m2, "m2")
save_vtk(m3, "m3")
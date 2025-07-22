using Documenter
using Plot3D

example_dict = [
    "O-grid Example" => ("examples/ogrid_fewerblocks.md"),
    "C-grid Example" => ("examples/cgrid_copiedcgrid.md"),
]

examples_md_dir = joinpath(@__DIR__,"src/examples")


makedocs(
    sitename = "Plot3D.jl",
    authors = "Maia Ross",
    modules = [Plot3D],
    pages = [
        "Home" => "index.md",
        # "Installation" => "installation.md",
        # "Usage" => "usage.md",
        "Examples" => "examples.md",
        # "Contributing" => "contributing.md",
        # "Contact" => "contact.md"
    ],
    format = Documenter.HTML(),
)From Plot3D folder:

julia ./test/Run_File3D.jl ./data/path_to_.xyz_or.p3d_file
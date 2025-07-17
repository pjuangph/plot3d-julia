using Documenter
using Plot3D

makedocs(
    sitename = "Plot3D.jl",
    authors = "Maia Ross",
    modules = [Plot3D],
    pages = [
        "Home" => "index.md",
        # "Installation" => "installation.md",
        # "Usage" => "usage.md",
        # "Examples" => "examples.md",
        # "Contributing" => "contributing.md",
        # "Contact" => "contact.md"
    ],
    format = Documenter.HTML(),
)
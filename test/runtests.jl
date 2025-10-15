using Test

# Keep each major area in its own file so you can run it standalone during debugging
include("test_read_write.jl")                    # your existing file (kept as-is)
include("test_translated_periodicity.jl")        # translated periodicity workflow (x,y,z)
include("test_glennht_io.jl")                    # GlennHT import/export smoke tests
include("test_gridpro_io.jl")                    # GridPro connectivity parsing smoke tests
include("test_axial_duplication.jl")

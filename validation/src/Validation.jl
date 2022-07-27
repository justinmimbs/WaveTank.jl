module Validation

using GLMakie

include("../../extra/plotting.jl")
include("data_extraction.jl")

to_path(file) = relpath(normpath(@__DIR__, "..", file), pwd())

include("SolitaryWave.jl")
include("SimpleBeach.jl")
include("ConicalIsland.jl")
include("MonaiValley.jl")

export SolitaryWave, SimpleBeach, ConicalIsland, MonaiValley

end # module

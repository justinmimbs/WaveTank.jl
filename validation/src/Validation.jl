module Validation

using GLMakie

include("../../extra/basics.jl")
include("../../extra/plotting.jl")
include("data_extraction.jl")

include("SolitaryWave.jl")
include("SimpleBeach.jl")
include("ConicalIsland.jl")
include("MonaiValley.jl")

export SolitaryWave, SimpleBeach, ConicalIsland, MonaiValley

end # module

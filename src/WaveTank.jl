module WaveTank

using Printf

include("core.jl")
include("model.jl")
include("run.jl")

export Grid, Model, step!, run!, Results

end # module

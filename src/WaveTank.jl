module WaveTank

using Printf

include("core.jl")
include("model.jl")
include("basics.jl")
include("run.jl")

export Grid, Model, step!, run!, Results
export solitary_wave, planar_beach, truncated_cone

end # module

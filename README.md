# WaveTank.jl

A depth-averaged, nonhydrostatic fluid model capable of simulating water waves as they run up on shore.

For a brief introduction, see this [JuliaCon 2022 talk](https://www.youtube.com/watch?v=QM-32FzMOcQ&list=PLP8iPy9hna6TRg6qJaBLJ-FRMi9Cp7gSX).

## Example

Below is a complete example that models a solitary wave running up on a planar beach. Output from the simulation is saved to a file, and the results are animated using Makie with Observables.

```julia
using WaveTank: Grid, Model, Results, run!, solitary_wave, planar_beach
using GLMakie
using Printf

# 1. define model

grid = Grid(
    (-20.0, 10.0),  # x bounds (m)
    (0.0, 5.0),     # y bounds (m)
    (300, 50)       # resolution
)
depth = 0.5         # maximum water depth (m)
wave = solitary_wave(0.15; h=depth, x0=-15.0) # (amplitude; h, x0)
basin = planar_beach(depth, 1 / 20) # (depth, slope) waterline is at x = 0
model = Model(;
    grid,
    bcx=:open,                  # boundary conditions at the x bounds
    h=(x, y) -> basin(x),       # depth of the basin
    eta=(x, y) -> wave.eta(x),  # initial surface elevation
    u=(x, y) -> wave.u(x),      # initial particle velocity in the x direction
)

# 2. run simulation

run!(model, "out/example.jld2";
    seconds=40.0,       # duration
    frequency=5,        # output frequency (per second)
    output=(; eta=:eta) # output fields
)

# 3. plot results

results = Results("out/example.jld2")
(; grid, h, dt) = results # results contain static attributes of the model

snapshot = Observable(results[0]) # results contain model snapshots, indexed by timestep
eta = @lift $snapshot.eta[:, end ÷ 2]
title = @lift @sprintf("time = %.1f s", $snapshot.t * dt)

fig = Figure(; resolution=(900, 300))
ax = Axis(fig[1, 1]; xlabel="x (m)", ylabel="surface elevation (m)", title)
lines!(ax, grid.xc, -h[:, 1]; color=:gray)
lines!(ax, grid.xc, eta; color=:dodgerblue)
display(fig)

# 4. animate figure

for s in results # results are iterable
    snapshot[] = s
    sleep(0.016)
end
```

For more examples, see the _validation_ directory.

## Installation

This package is unregistered and it currently lacks documentation, but if you're interested, you can add it by URL.

```
pkg> add https://github.com/justinmimbs/WaveTank.jl
```

----

## References

The model used in WaveTank.jl is based on the one formulated in:

Yamazaki, Y., Kowalik, Z., & Cheung, K.F. (2009). Depth-integrated, non-hydrostatic model for wave breaking and run-up. _International Journal for Numerical Methods in Fluids_, 61, 473-497. https://doi.org/10.1002/fld.1952

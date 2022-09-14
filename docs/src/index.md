# WaveTank.jl

_A depth-averaged, nonhydrostatic fluid model capable of simulating water waves as they run up on shore._

This page presents an overview of how to
- define a model by its domain and initial conditions,
- define wave generator boundary conditions,
- run a simulation while saving output to a file, and
- access results in the output file.

## Defining a model

The model uses meters for all distances and seconds for all durations.

```@docs
Grid
Model
```

The still-water level, `h`, defines the bathymetry over the domain. Positive and negative heights indicate submerged and emerged surfaces, respectively.

The surface elevation, `eta`, is the distance from still-water level. Here, positive and negative values indicate above and below, respectively.

The depth of the water, from the bottom to the surface, is given by `h + eta`.

## Defining a wave generator

Waves can be generated from the lower x boundary only. Currently a wave generator must be specified as a timeseries of either one or both of
- particle velocity in the x direction, or
- surface elevation.

A timeseries must be a `Vector{Float64}`; each value is applied to the entire lower x boundary for a single timestep.

## Running a simulation

```@docs
step!
run!
```

### Defining output fields

The `output` parameter of `run!` defines fields to write to the file. It can be given as a `NamedTuple` or any iterator of key-value pairs.

Output fields will be written at the interval implied by the given `frequency`, and this interval may be greater than the timestep interval. Because multiple timesteps may pass between writing output, we may want to output a _snapshot_ from the current timestep, or we may want to output an _aggregate_ of all timesteps since the last output time.

To take snapshots from the model at output time, define a field as one of the following.
- `name::Symbol`: the name of a model field to copy (one of `[:eta, :u, :v]`)
- `f::Function`: a function to map from the model to an arbitrary value

To make aggregates using models at all timesteps between each output time, define a field as such.
- `(; f, combine)::NamedTuple`: a function to map from each model to an arbitrary value, and a function to combine these values

For example, to aggregate the maximum water depth reached at each grid cell over the duration of each interval and write it to field named `max_depth`, you could define `output` as follows.

```julia
output = (; max_depth=(f=m -> m.h .+ m.eta, combine=(a, b) -> max.(a, b)))
```

## Working with output

To make it easier to access the output from a simulation run, we can use the `Results` type. This type represents the sequence of output as it was written during the simulation; it is iterable and it can be indexed by timestep number.

```julia
results = Results("path/to/example.jld2")
for x in results
    # `x` is a named tuple containing the output fields
    # for convenience, `x` also contains static fields from the model
end
results[end] # get the output at the last output time
```

```@docs
Results(filepath)
Results
```

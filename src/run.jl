using JLD2: jldsave, jldopen

"""
    step!(model[, iterations])

Step a model forward in time.
"""
function step!(m, n::Int; wavemaker=[], waveinput=[])
    scale = m.bcx == :open ? 2.0 : 1.0
    for _ in 1:n
        t = m.t + 1
        if t < length(wavemaker)
            m.u[1, :] .+= scale * (wavemaker[t + 1] - wavemaker[t])
        end
        step!(m)
        if t <= length(waveinput)
            m.eta[1, :] .= waveinput[t]
        end
    end
    m
end

# simulation output

GridSerialized = NamedTuple{
    (:xmin, :xmax, :ymin, :ymax, :nx, :ny),
    Tuple{Float64, Float64, Float64, Float64, Int64, Int64}
}
function Base.convert(::Type{Grid}, a::GridSerialized)
    Grid((a.xmin, a.xmax), (a.ymin, a.ymax), (a.nx, a.ny))
end
function Base.convert(::Type{GridSerialized}, grid::Grid)
    (
        xmin=grid.xf[1],
        xmax=grid.xf[end],
        ymin=grid.yf[1],
        ymax=grid.yf[end],
        nx=grid.nx,
        ny=grid.ny,
    )
end

"""
    run!(model, filepath; seconds, frequency, output, wavemaker=[], waveinput=[])

Run a model forward in time while writing output to a file.

## Keyword arguments
- `seconds`: duration to run the simulation for
- `frequency`: number of times per simulation-second write output
- `output`: key-value pairs indicating the name of the output field and how to compute it
- `wavemaker=[]`: timeseries of particle velocity to apply to the lower x boundary
- `waveinput=[]`: timeseries of surface elevation to apply to the lower x boundary
"""
function run!(m, filepath; seconds, frequency, output,
        wavemaker=[], waveinput=[],
    )
    @assert frequency <= 1.0 / m.dt
    # wavemaker from waveinput if needed
    if isempty(wavemaker) && !isempty(waveinput)
        wavemaker = particle_velocity(waveinput, maximum(m.h[1, :]))
    end
    # split output fields into snapshots and aggregates
    snapshots = [] # computed at the end of each output interval
    aggregates = [] # computed at every timestep and combined over the interval
    for (key, v) in pairs(output)
        if isa(v, Symbol) && hasproperty(m, v)
            push!(snapshots, key => m -> getproperty(m, v))
        elseif isa(v, Function)
            push!(snapshots, key => v)
        elseif isa(v, NamedTuple) && haskey(v, :f) && haskey(v, :combine)
            push!(aggregates, key => v)
        else
            throw(ArgumentError("invalid value for output.$key'"))
        end
    end
    # write base data to file
    base = Dict(
        :grid => convert(GridSerialized, m.grid),
        :dt => m.dt,
        :h => m.h,
    )
    mkpath(dirname(filepath))
    jldsave(filepath; base...)
    # write initial output to file (t = 0)
    jldopen(filepath, "r+") do file
        for (key, f) in snapshots
            write(file, "output/$key/$(m.t)", f(m))
        end
        for (key, (; f)) in aggregates
            write(file, "output/$key/$(m.t)", f(m))
        end
    end
    # for n output intervals
    n = ceil(Int, seconds * frequency)
    for i in 1:n
        ti = round(Int, (i / frequency) / m.dt)
        aggregated = Dict()
        for t in (m.t + 1):ti
            # step
            isempty(wavemaker) ? step!(m) : step!(m, 1; wavemaker, waveinput)
            # check for NaN, Inf
            if !all(isfinite, m.eta)
                @error "encountered non-finite float: stopping at t = $(m.t)"
                return nothing
            end
            # aggregate
            for (key, (; f, combine)) in aggregates
                r = haskey(aggregated, key) ? combine(aggregated[key], f(m)) : f(m)
                aggregated[key] = r
            end
        end
        # write output to file (t = ti)
        jldopen(filepath, "r+") do file
            for (key, f) in snapshots
                write(file, "output/$key/$(m.t)", f(m))
            end
            for (key, r) in aggregated
                write(file, "output/$key/$(m.t)", r)
            end
        end
        filename = basename(filepath)
        @info "saved output to '$filename' at t = $(m.t) (interval $i of $n)"
    end
    m
end

# Results

"""
Represent the output from a simulation run.

Static fields from the model (`grid`, `dt`, `h`) are available directly.

## Fields
- `filepath::String`
- `grid::Grid`
- `dt::Float64`
- `h::Matrix{Float64}`
- `ts::Vector{Int}`: timestep numbers for all output times
- `output::Vector{Symbol}`: keys of the fields available at each output time
"""
struct Results
    filepath::String
    grid::Grid
    dt::Float64
    h::Matrix{Float64}
    output::Vector{Symbol}
    ts::Vector{Int}
    #
    @doc """
        Results(filepath)

    Load the results from a file.
    """
    function Results(filepath)
        jldopen(filepath, "r") do file
            grid = convert(Grid, file["grid"])
            dt = file["dt"]
            h = file["h"]
            output = Symbol.(keys(file["output"]))
            ts = sort!(parse.(Int, keys(file["output/$(output[1])"])))
            new(filepath, grid, dt, h, output, ts)
        end
    end
end

function Base.show(io::IO, r::Results)
    (; filepath, dt, ts) = r
    @printf(io, "WaveTank.Results \"%s\" (%g s)", filepath, dt * ts[end])
end
function Base.show(io::IO, ::MIME"text/plain", r::Results)
    (; filepath, dt, ts, output) = r
    @printf(io, "WaveTank.Results \"%s\" (%g s)", filepath, dt * ts[end])
    @printf(io, "\n   %s", format_ts(ts))
    @printf(io, "\n   output = [%s]", join(output, ", "))
end

function format_ts(ts)
    if 2 < length(ts)
        "$(length(ts)) snapshots @ t = [0, $(ts[2]), ..., $(ts[end])]"
    elseif length(ts) == 2
        "2 snapshots @ t =[0, $(ts[end])]"
    else
        "1 snapshot @ t = 0"
    end
end

# interface: iterator (over timeseries data)

Base.length(r::Results) = length(r.ts)
Base.eltype(::Type{Results}) = NamedTuple

function Base.iterate(r::Results)
    iterate(r, (jldopen(r.filepath, "r"), 0))
end
function Base.iterate(r::Results, (file, i))
    if i < length(r.ts)
        t = r.ts[i + 1]
        v = (; r.grid, r.dt, r.h, t, ( key => file["output/$key/$t"] for key in r.output )...)
        v, (file, i + 1)
    else
        close(file)
        nothing
    end
end

# interface: indexing (by timestep)

function Base.getindex(r::Results, t::Int)
    insorted(t, r.ts) || throw(BoundsError())
    jldopen(r.filepath, "r") do file
        (; r.grid, r.dt, r.h, t, ( key => file["output/$key/$t"] for key in r.output )...)
    end
end

Base.firstindex(r::Results) = 0
Base.lastindex(r::Results) = r.ts[end]

module WaveTank

using JLD2: jldsave, jldopen

include("core.jl")
include("model.jl")

# step! multiple iterations, optionally with wavemaker/waveinput
# 	wavemaker = velocity timeseries to apply to lower x boundary
# 	waveinput = eta timeseries to apply to lower x boundary
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

function velocity_from_eta(eta, h)
    # particle speed / phase speed = amplitude / depth
    [ a / (h + a) * sqrt(g * (h + a)) for a in eta ]
end

function run!(m, filepath; seconds, frequency, output,
        wavemaker=[], waveinput=[],
    )
    @assert frequency <= 1.0 / m.dt
    # wavemaker from waveinput if needed
    if isempty(wavemaker) && !isempty(waveinput)
        wavemaker = velocity_from_eta(waveinput, maximum(m.h[1, :]))
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
                @warn "encountered non-finite float: stopping at t = $(m.t)"
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
        # filename = basename(filepath)
        # @info "run! saved output to '$filename' at t = $(m.t) (interval $i of $n)"
    end
    m
end

# TODO remove foldtime, maptime; use map/foldl(Result(filepath)) instead
# TODO remove load; implement indexing for Result instead: Result(filepath)[t]

# extract a named tuple from a jld2 file at time step t
function load(filepath::AbstractString, t::Int)
    jldopen(filepath, "r") do file
        base = (grid=convert(Grid, file["grid"]), dt=file["dt"], h=file["h"], t=t)
        output = Base.map(Symbol, keys(file["output"]))
        (; base..., ( key => file["output/$key/$t"] for key in output )...)
    end
end
# return all timesteps available in output
function load(filepath::AbstractString)
    jldopen(filepath, "r") do file
        output = keys(file["output"])
        sort!(parse.(Int, keys(file["output/$(output[1])"])))
    end
end

function foldtime(f::Function, filepath; init)
    jldopen(filepath, "r") do file
        base = (grid=convert(Grid, file["grid"]), dt=file["dt"], h=file["h"])
        output = Base.map(Symbol, keys(file["output"]))
        ts = sort!(parse.(Int, keys(file["output/$(output[1])"])))
        foldl(ts, init=init) do r, t
            m = (; base..., t, ( key => file["output/$key/$t"] for key in output )...)
            f(r, m)
        end
    end
end

function maptime(f::Function, filepath)
    foldtime((a, m) -> push!(a, f(m)), filepath, init=[])
end

# Results

# respresent results from a saved simulation run

struct Results
    filepath::AbstractString
    grid::Grid
    dt::Float64
    h::Matrix{Float64}
    output::Vector{Symbol}
    ts::Vector{Int}
    #
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

# iterator interface to timeseries data

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

end # module

import Pkg
Pkg.activate(Base.current_project(@__DIR__))

using BenchmarkTools
using ProfileView
using WaveTank: g, Grid, Model, step!

include("../extra/basics.jl")

function model()
    h = 1.0
    a = 0.3 * h
    k = sqrt(3 * a / (4 * h^3))
    l = (2.0 / k) * acosh(sqrt(0.05^-1))
    beta = 1 / 19.85 # beach slope
    # beach (waterline at x = 0)
    bmin = -h / beta
    bh = 2 * h # beach height
    bmax = bmin + (bh / beta)
    # x domain
    xmin = bmin - 1.5 * l
    dx = h / 8
    nx = ceil(Int, (bmax - xmin) / dx)
    xmax = xmin + nx * dx
    #
    wave = sech2wave(a=a, x0=bmin - 0.5 * l, k=k)
    speed = particle_velocity(wave, h, a)
    bathymetry = piecewiselinear([bmin, bmax, bmax + l], h .- [0, bh, bh], h)
    Model(
        grid=Grid((xmin, xmax), (0, 80dx), (nx, 80)),
        h=(x, y) -> bathymetry(x),
        bcx=:wall,
        bcy=:wall,
        eta=(x, y) -> wave(x),
        u=(x, y) -> speed(x),
        n=0.01,
        dt=0.01,
    )
end

let
    m = step!(model())
    @info typeof(m.solver), hasproperty(m, :lu) ? typeof(m.lu) : nothing
end

function benchmark(samples=10)
    b1 = @benchmarkable step!(m, 200) setup=(m = model()) #teardown=(println(m.t))
    run(b1; evals=1, samples, seconds=300.0, gcsample=true)
end

function profile(steps=1)
    ProfileView.@profview step!(model(), steps)
end
